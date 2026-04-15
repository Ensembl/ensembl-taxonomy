# Copyright [2009-2026] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Module for working with NCBI Taxonomy taxdumps."""

__all__ = [
    "TaxonCache",
]

from collections import defaultdict
from collections.abc import MutableMapping
from contextlib import ExitStack
import datetime
from ftplib import FTP
import itertools
import os
from pathlib import Path
import re
import shutil
from tempfile import TemporaryDirectory
from typing import Optional, Union
from warnings import warn
from zipfile import ZipFile, ZIP_DEFLATED

import duckdb
import numpy as np
import pandas as pd
import requests

from ensembl.taxonomy.db import _NCBI_TAXA_TABLE_NAMES, TaxonDB
from ensembl.taxonomy.utils import group_by_array, inclusive_range_between_values, NcbiTaxdumpDialect


_NCBI_TAXDUMP_METAINFO: dict = {
    "merged.dmp": {
        "canonic_columns": ["old_taxon_id", "new_taxon_id"],
        "canonic_sort_keys": ["old_taxon_id", "new_taxon_id"],
        "dtypes": {
            "old_taxon_id": np.int32,
            "new_taxon_id": np.int32,
        },
        # We treat dump file delimiter as "\t", so we need to double
        # column indices w.r.t the original pipe-delimited columns.
        "dump_col_map": {
            0: "old_taxon_id",
            2: "new_taxon_id",
        },
    },
    "names.dmp": {
        "canonic_columns": ["taxon_id", "name", "name_class"],
        "canonic_sort_keys": ["taxon_id", "name_class", "name"],
        "dtypes": {
            "taxon_id": np.int32,
            "name": np.str_,
            "name_class": "category",
        },
        # We treat dump file delimiter as "\t", so we need to double
        # column indices w.r.t the original pipe-delimited columns.
        "dump_col_map": {
            0: "taxon_id",
            2: "name",
            6: "name_class",
        },
    },
    "nodes.dmp": {
        "canonic_columns": ["taxon_id", "parent_id", "rank", "genbank_hidden_flag"],
        "canonic_sort_keys": ["taxon_id"],
        "dtypes": {
            "taxon_id": np.int32,
            "parent_id": np.int32,
            "rank": "category",
            "genbank_hidden_flag": np.int8,
        },
        # We treat dump file delimiter as "\t", so we need to double
        # column indices w.r.t the original pipe-delimited columns.
        "dump_col_map": {
            0: "taxon_id",
            2: "parent_id",
            4: "rank",
            20: "genbank_hidden_flag",
        },
    },
}

_NH_RANKS = {"clade", "no rank"}

# Most of the rank hierarchy information presented here is as described in
# Schoch et al. (2020) NCBI Taxonomy: a comprehensive update on curation, resources and tools.
# <https://doi.org/10.1093/database/baaa062>,
# with some updates based on NCBI Insights (2024-06-04) Upcoming changes to NCBI Taxonomy classifications.
# <https://ncbiinsights.ncbi.nlm.nih.gov/2024/06/04/changes-ncbi-taxonomy-classifications/>
# and NCBI Insights (2025-04-25) NCBI Taxonomy updates to virus classification.
# <https://ncbiinsights.ncbi.nlm.nih.gov/2025/04/25/ncbi-taxonomy-updates-virus-classification-april-2025/>.
_RANK_HIERARCHY = [
    ["isolate"],
    ["strain"],
    ["serotype", "biotype", "genotype"],
    ["serogroup", "pathogroup"],
    ["forma"],
    ["subvariety"],
    ["varietas", "form", "morph"],
    ["subspecies"],
    ["forma specialis", "special form"],
    ["species"],
    ["species subgroup"],
    ["species group"],
    ["subseries"],
    ["series"],
    ["subsection"],
    ["section"],
    ["subgenus"],
    ["genus"],
    ["subtribe"],
    ["tribe"],
    ["subfamily"],
    ["family"],
    ["superfamily"],
    ["parvorder"],
    ["infraorder"],
    ["suborder"],
    ["order"],
    ["superorder"],
    ["subcohort"],
    ["cohort"],
    ["infraclass"],
    ["subclass"],
    ["class"],
    ["superclass"],
    ["infraphylum", "infradivision"],
    ["subphylum", "subdivision"],
    ["phylum", "division"],
    ["superphylum", "superdivision"],
    ["subkingdom"],
    ["kingdom"],
    ["domain", "realm"],
    ["acellular root", "cellular root"],
]


def _build_aux_node_data(node_df: pd.DataFrame, root_id: int) -> tuple[pd.DataFrame, set[str]]:
    leftright_df, unorderable_ranks = _build_leftright_indices(node_df, root_id)
    leftright_df.set_index("taxon_id", drop=True, inplace=True)
    node_df = node_df.join(leftright_df, on="taxon_id", how="inner", validate="one_to_one")
    node_df.sort_values("left_index", inplace=True)
    node_df.reset_index(drop=True, inplace=True)
    return node_df, unorderable_ranks


def _build_leftright_indices(node_df: pd.DataFrame, root_taxon_id: int) -> tuple[pd.DataFrame, set[str]]:
    node_ranks = dict(zip(node_df["taxon_id"].tolist(), node_df["rank"].tolist()))
    rank_data: dict = {"pairs": defaultdict(int), "ranks_by_node": node_ranks, "stack": []}

    parent_child_map = group_by_array(
        node_df["taxon_id"].to_numpy(),
        node_df["parent_id"].to_numpy(),
        sort_values=True,
    )

    left_right_recs: list[tuple[int, int, int]] = []
    _recursive_leftright_indexing(parent_child_map, rank_data, left_right_recs, root_taxon_id)
    leftright_col_names = ["taxon_id", "left_index", "right_index"]
    leftright_df = pd.DataFrame(left_right_recs, columns=leftright_col_names, dtype=np.int32)

    unorderable_ranks = _find_unorderable_ranks(rank_data["pairs"])
    return leftright_df, unorderable_ranks


def _current_timestamp_utc() -> str:
    curr_time_utc = datetime.datetime.now(tz=datetime.timezone.utc)
    return curr_time_utc.strftime("%Y-%m-%dT%H:%M:%SZ")


def _find_unorderable_ranks(rank_pair_counts: dict[tuple[str, str], int]) -> set[str]:

    exp_rank_pairs = set()
    for desc_ranks, anc_ranks in itertools.combinations(_RANK_HIERARCHY, 2):
        for desc_rank, anc_rank in itertools.product(desc_ranks, anc_ranks):
            exp_rank_pairs.add((desc_rank, anc_rank))

    # From the observed set of rank pairs, we flag pairs
    # of ranks inconsistent with the expected hierarchy.
    obs_rank_pairs = set(rank_pair_counts.keys())
    flagged_rank_pairs = obs_rank_pairs - exp_rank_pairs

    unorderable_ranks = set()
    while flagged_rank_pairs:
        flagged_ranks = {rank for rank_pair in flagged_rank_pairs for rank in rank_pair}

        rank_to_exp_pair_count: MutableMapping[str, int] = defaultdict(int)
        rank_to_odd_pair_count: MutableMapping[str, int] = defaultdict(int)
        for rank_pair in obs_rank_pairs:
            if rank_pair in flagged_rank_pairs:
                for rank in rank_pair:
                    if rank in flagged_ranks:
                        rank_to_odd_pair_count[rank] += rank_pair_counts[rank_pair]
            else:
                for rank in rank_pair:
                    if rank in flagged_ranks:
                        rank_to_exp_pair_count[rank] += rank_pair_counts[rank_pair]

        # From the set of ranks involved in a flagged rank pair, we try to select as unorderable
        # the unique rank associated with the maximum number of unexpected rank pairs.
        max_num_odd_rank_pairs = max(rank_to_odd_pair_count.values())
        maximal_odd_ranks = [
            rank for rank in flagged_ranks if rank_to_odd_pair_count[rank] == max_num_odd_rank_pairs
        ]

        if len(maximal_odd_ranks) == 1:
            unorderable_rank = maximal_odd_ranks[0]
        else:
            # If multiple flagged ranks are involved in the maximum number of unexpected rank
            # pairs, we take from these the rank with the minimum number of expected rank pairs.
            min_num_exp_rank_pairs = min(rank_to_exp_pair_count.values())
            minimal_exp_ranks = [
                rank for rank in maximal_odd_ranks if rank_to_exp_pair_count[rank] == min_num_exp_rank_pairs
            ]

            if len(minimal_exp_ranks) == 1:
                unorderable_rank = minimal_exp_ranks[0]
            else:
                # If all else fails, we arbitrarily take the first-sorting rank name.
                unorderable_rank = min(minimal_exp_ranks)

        unorderable_ranks.add(unorderable_rank)

        obs_rank_pairs = {rank_pair for rank_pair in obs_rank_pairs if unorderable_rank not in rank_pair}
        flagged_rank_pairs = obs_rank_pairs - exp_rank_pairs

    return unorderable_ranks


def _recursive_leftright_indexing(
    parent_child_map: dict[int, list[int]],
    rank_data: dict,
    left_right_recs: list[tuple[int, int, int]],
    taxon_id: int,
    counter: int = 1,
) -> int:

    taxon_rank = rank_data["ranks_by_node"][taxon_id]
    taxon_rank_is_hierarchical = taxon_rank not in _NH_RANKS
    if taxon_rank_is_hierarchical:
        rank_data["stack"].append(taxon_rank)

    left_index = counter
    counter += 1

    if taxon_id in parent_child_map:
        for child_id in parent_child_map[taxon_id]:
            counter = _recursive_leftright_indexing(
                parent_child_map,
                rank_data,
                left_right_recs,
                child_id,
                counter,
            )

    right_index = counter
    counter += 1
    left_right_recs.append((taxon_id, left_index, right_index))

    if taxon_rank_is_hierarchical:
        taxon_rank = rank_data["stack"].pop()
        for anc_rank in rank_data["stack"]:
            if taxon_rank != anc_rank:
                rank_data["pairs"][(taxon_rank, anc_rank)] += 1

    return counter


def _taxdump_id_from_last_modified(last_modified: str) -> str:
    # We assume the last-modified header has
    # English-language month abbreviations.
    month_abbr_to_num = {
        "Jan": "01",
        "Feb": "02",
        "Mar": "03",
        "Apr": "04",
        "May": "05",
        "Jun": "06",
        "Jul": "07",
        "Aug": "08",
        "Sep": "09",
        "Oct": "10",
        "Nov": "11",
        "Dec": "12",
    }

    timestamp_re = re.compile(
        r"[^,]+, (?P<day>[0-9]{2}) (?P<month_abbr>\S+) (?P<year>[0-9]{4,})"
        r" (?P<hour>[0-9]{2}):(?P<min>[0-9]{2}):(?P<sec>[0-9]{2}) GMT"
    )

    if m := timestamp_re.fullmatch(last_modified):
        month = month_abbr_to_num[m["month_abbr"]]
        taxdump_id = f"{m['year']}-{month}-{m['day']}T{m['hour']}:{m['min']}:{m['sec']}Z"
    else:
        raise ValueError(f"failed to generate taxdump ID from 'Last-Modified' header: '{last_modified}'")

    return taxdump_id


class TaxonCache:
    """A class for caching NCBI Taxonomy taxdumps."""

    def __init__(self, path: Union[os.PathLike, str], threads: int = 1) -> None:
        self._threads = threads

        self._cache_root_path = Path(path)
        self._taxdump_root_path = self._cache_root_path / "taxdumps"
        self._taxdump_archive_meta_file = self._cache_root_path / "archive_meta.parquet"

        try:
            self._archive_meta = pd.read_parquet(self._taxdump_archive_meta_file)
        except FileNotFoundError:
            self._archive_meta = self._update_archive_metadata()

    def _get_archive_metadata(self) -> pd.DataFrame:
        latest_archive_date = self._archive_meta["archive_id"].iat[-1]

        # We expect the archive list needs updating if there has been
        # no update in the 4 weeks prior to the first of this month.
        utc_time_now = datetime.datetime.now(datetime.timezone.utc)
        first_of_this_month = np.datetime64(utc_time_now.date().replace(day=1))
        if first_of_this_month - np.datetime64(latest_archive_date) >= np.timedelta64(28, "D"):
            self._archive_meta = self._update_archive_metadata()

        return self._archive_meta

    def _get_taxdump_dir_path(self, taxdump_id: str) -> Path:
        taxdump_path_part = self._taxdump_id_to_path_part(taxdump_id)
        taxdump_dt = self._taxdump_id_to_datetime(taxdump_id)
        year_path_part = str(taxdump_dt.year)
        taxdump_dir_path = self._taxdump_root_path / year_path_part / taxdump_path_part
        return taxdump_dir_path

    def _load_data_from_taxdump(
        self, taxdump_path: Path, table_name: str, taxdump_id: str, archive_id: Optional[str] = None
    ) -> pd.DataFrame:

        if table_name == "ncbi_taxa_name":
            dump_file_names = ["names.dmp", "merged.dmp"]
        elif table_name == "ncbi_taxa_node":
            dump_file_names = ["nodes.dmp"]
        else:
            raise ValueError(f"unknown NCBI Taxonomy table name: '{table_name}'")

        for dump_file_name in dump_file_names:

            metainfo = _NCBI_TAXDUMP_METAINFO[dump_file_name]
            col_idxs = sorted(metainfo["dump_col_map"])

            col_names = []
            dump_dtypes = {}
            for col_idx in col_idxs:
                col_name = metainfo["dump_col_map"][col_idx]
                dump_dtypes[col_idx] = metainfo["dtypes"][col_name]
                col_names.append(col_name)

            with ZipFile(taxdump_path, compression=ZIP_DEFLATED) as zip_file_obj:
                with zip_file_obj.open(dump_file_name, mode="r") as in_file_obj:
                    dump_file_df = pd.read_csv(  # type: ignore
                        in_file_obj,
                        header=None,
                        index_col=None,
                        usecols=col_idxs,
                        names=col_names,
                        dtype=dump_dtypes,
                        keep_default_na=False,
                        na_filter=False,
                        compression=None,
                        lineterminator=NcbiTaxdumpDialect.lineterminator,
                        encoding="utf-8",
                        dialect=NcbiTaxdumpDialect,
                    )

            if dump_file_name == "names.dmp":
                df = dump_file_df

                root_mask = (df["name_class"] == "scientific name") & (df["name"] == "root")
                root_id, *surplus_root_ids = df.loc[root_mask, "taxon_id"].tolist()

                if surplus_root_ids:
                    raise ValueError(f"{len(surplus_root_ids)} surplus root IDs found in 'nodes.dmp'")

                provenance_recs = [
                    [root_id, _current_timestamp_utc(), "import date"],
                    [root_id, taxdump_id, "taxdump_id"],
                ]
                if archive_id:
                    provenance_recs.append([root_id, archive_id, "archive_id"])
                provenance_df = pd.DataFrame(provenance_recs, columns=col_names)

                df = pd.concat([provenance_df, df])

            elif dump_file_name == "merged.dmp":

                df_part = pd.concat(
                    [dump_file_df["new_taxon_id"], dump_file_df["old_taxon_id"]],
                    axis="columns",
                    keys=["taxon_id", "name"],
                )
                df_part["name_class"] = "merged_taxon_id"
                df = pd.concat([df, df_part], axis="index")
                df.sort_values(["taxon_id", "name_class"], inplace=True)
                df.reset_index(drop=True, inplace=True)

            elif dump_file_name == "nodes.dmp":
                df = dump_file_df

                root_mask = df["taxon_id"] == df["parent_id"]
                root_id, *surplus_root_ids = df.loc[root_mask, "taxon_id"].tolist()

                if surplus_root_ids:
                    raise ValueError(f"{len(surplus_root_ids)} surplus root IDs found in 'nodes.dmp'")

                df.loc[root_mask, "parent_id"] = 0
                df, unorderable_ranks = _build_aux_node_data(df, root_id)
                df["root_id"] = root_id

                for rank in unorderable_ranks:
                    warn(f"rank '{rank}' cannot be ordered consistently")
            else:
                raise ValueError(f"unknown taxdump file name '{dump_file_name}'")

        return df

    def _taxdump_id_from_path_part(self, path_part: str) -> str:

        path_part_re = re.compile(
            r"(?P<year>[0-9]{4,})(?P<month>[0-9]{2})(?P<day>[0-9]{2})"
            r"T(?P<hour>[0-9]{2})(?P<min>[0-9]{2})(?P<sec>[0-9]{2})Z"
        )

        if match := path_part_re.fullmatch(path_part):
            year, month, day, hour, minute, second = match.groups()
            taxdump_id = f"{year}-{month}-{day}T{hour}:{minute}:{second}Z"
        else:
            raise ValueError(f"failed to generate taxdump ID from path part: '{path_part}'")

        return taxdump_id

    def _taxdump_id_to_datetime(self, taxdump_id: str) -> datetime.datetime:
        try:
            taxdump_dt = datetime.datetime.strptime(taxdump_id, "%Y-%m-%dT%H:%M:%SZ")
        except ValueError as exc:
            raise ValueError(f"failed to obtain datetime for taxdump '{taxdump_id}'") from exc
        return taxdump_dt.replace(tzinfo=datetime.timezone.utc)

    def _taxdump_id_to_path_part(self, taxdump_id: str) -> str:

        taxdump_id_re = re.compile(
            r"(?P<year>[0-9]{4})-(?P<month>[0-9]{2})-(?P<day>[0-9]{2})"
            r"T(?P<hour>[0-9]{2}):(?P<min>[0-9]{2}):(?P<sec>[0-9]{2})Z"
        )

        if match := taxdump_id_re.fullmatch(taxdump_id):
            year, month, day, hour, minute, second = match.groups()
            path_part = f"{year}{month}{day}T{hour}{minute}{second}Z"
        else:
            raise ValueError(f"failed to generate path part from taxdump ID: '{taxdump_id}'")

        return path_part

    def _update_archive_metadata(self):
        archive_file_re = re.compile(
            r"(?P<prefix>taxdmp|new_taxdump)_(?P<date>[0-9]{4,}-[0-9]{2}-[0-9]{2})\.zip"
        )

        with FTP("ftp.ncbi.nlm.nih.gov") as ftp:
            ftp.login()
            ftp.cwd("pub/taxonomy/taxdump_archive")
            archive_file_recs = list(ftp.mlsd(facts=["modify"]))

        archive_meta_recs = []
        for file_name, file_rec in archive_file_recs:
            if match := archive_file_re.fullmatch(file_name):
                if match["prefix"] == "taxdmp":
                    archive_id = match["date"]
                    taxdump_dt = datetime.datetime.strptime(file_rec["modify"], "%Y%m%d%H%M%S").replace(
                        tzinfo=datetime.timezone.utc
                    )
                    taxdump_id = taxdump_dt.strftime("%Y-%m-%dT%H:%M:%SZ")
                    archive_meta_recs.append([archive_id, taxdump_id])
            elif not file_name.startswith("."):
                warn(
                    f"unrecognised taxdump file name: '{file_name}';"
                    f" taxdump archive listing may be incomplete"
                )

        archive_col_names = ["archive_id", "taxdump_id"]
        archive_meta_df = (
            pd.DataFrame.from_records(archive_meta_recs, columns=archive_col_names)
            .sort_values("taxdump_id")
            .set_index("archive_id", drop=False)
        )

        if not (hasattr(self, "_archive_meta") and archive_meta_df.equals(self._archive_meta)):
            os.makedirs(self._taxdump_archive_meta_file.parent, exist_ok=True)
            archive_meta_df.to_parquet(self._taxdump_archive_meta_file, compression="zstd")

        return archive_meta_df

    def _write_dataframe_to_parquet(self, _dataframe: pd.DataFrame, parquet_file_path: Path) -> None:
        with duckdb.connect(config={"threads": self._threads}) as conn:
            conn.execute("CREATE TABLE t AS SELECT * FROM _dataframe")
            conn.execute(f"COPY t TO '{parquet_file_path}' (FORMAT PARQUET, COMPRESSION 'zstd')")

    def available_archive_ids(self) -> list[str]:
        """Return list of available NCBI Taxonomy archives."""
        archive_meta_df = self._get_archive_metadata()
        return archive_meta_df["archive_id"].tolist()

    def cache_archives(self, first_archive_id: str, final_archive_id: str) -> dict[str, str]:
        """Cache NCBI Taxonomy taxdump archives between (and including) those specified.

        Args:
            first_archive_id: Taxdump ID of first archived taxdump to cache.
            final_archive_id: Taxdump ID of final archived taxdump to cache.

        Returns:
            A dictionary mapping, for each of the cached archives,
            the relevant archive ID to its corresponding taxdump ID.
        """
        avail_archive_ids = self.available_archive_ids()
        archive_idx_range = inclusive_range_between_values(
            avail_archive_ids, first_archive_id, final_archive_id
        )

        archives_cached = {}
        for archive_idx in archive_idx_range:
            archive_id = avail_archive_ids[archive_idx]
            taxdump_id = self.cache_taxdump(archive_id=archive_id)
            archives_cached[archive_id] = taxdump_id

        return archives_cached

    def cache_taxdump(self, archive_id: Optional[str] = None) -> str:
        """Cache the current NCBI Taxonomy taxdump, or the given archive, if specified."""

        if archive_id:
            taxdump_base_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/"
            taxdump_file_name = f"taxdmp_{archive_id}.zip"
        else:
            taxdump_base_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
            taxdump_file_name = "taxdmp.zip"

        taxdump_url = taxdump_base_url + taxdump_file_name

        with ExitStack() as stack:
            with requests.get(taxdump_url, stream=True, timeout=60) as response:
                response.raise_for_status()

                last_modified = response.headers["Last-Modified"]
                taxdump_id = _taxdump_id_from_last_modified(last_modified)
                taxdump_dir_path = self._get_taxdump_dir_path(taxdump_id)

                taxdump_dir_path.parent.mkdir(mode=0o775, parents=True, exist_ok=True)
                tmp_dir = stack.enter_context(
                    TemporaryDirectory(dir=taxdump_dir_path.parent, prefix=".taxdump_")
                )
                tmp_dir_path = Path(tmp_dir)
                tmp_taxdump_path = tmp_dir_path / taxdump_file_name
                with open(tmp_taxdump_path, "wb") as out_file_obj:
                    shutil.copyfileobj(response.raw, out_file_obj)

                exp_file_size = int(response.headers["Content-Length"])
                obs_file_size = os.path.getsize(tmp_taxdump_path)
                if obs_file_size != exp_file_size:
                    raise RuntimeError(
                        f"observed taxdump size ({obs_file_size} bytes) does not"
                        f" match expected content length ({exp_file_size} bytes)"
                    )

            file_pairs = []
            taxon_ids_by_table = {}
            for table_name in _NCBI_TAXA_TABLE_NAMES:
                pqt_file_name = f"{table_name}.parquet"
                tmp_pqt_file_path = tmp_dir_path / pqt_file_name
                out_pqt_file_path = taxdump_dir_path / pqt_file_name
                df = self._load_data_from_taxdump(
                    tmp_taxdump_path, table_name, taxdump_id, archive_id=archive_id
                )
                self._write_dataframe_to_parquet(df, tmp_pqt_file_path)
                file_pairs.append((tmp_pqt_file_path, out_pqt_file_path))
                taxon_ids_by_table[table_name] = set(df["taxon_id"].unique())

            name_taxon_ids = taxon_ids_by_table["ncbi_taxa_name"]
            node_taxon_ids = taxon_ids_by_table["ncbi_taxa_node"]
            if name_taxon_ids != node_taxon_ids:
                raise ValueError(
                    f"Taxonomy ID mismatch: 'ncbi_taxa_name' ({len(name_taxon_ids)})"
                    f" vs 'ncbi_taxa_node' ({len(node_taxon_ids)})"
                )

            try:
                taxdump_dir_path.mkdir(mode=0o775, exist_ok=False)
            except FileExistsError:
                warn(f"skipping cache of taxdump '{taxdump_id}' as it is already cached")
            else:
                for tmp_file_path, dst_file_path in file_pairs:
                    shutil.move(tmp_file_path, dst_file_path)

        return taxdump_id

    def cached_taxdump_ids(self) -> list[str]:
        """List cached NCBI Taxonomy taxdumps.

        Returns:
            A list containing the taxdump IDs of cached taxdumps.
        """
        cached_taxdump_ids = []
        for node_file_path in self._taxdump_root_path.glob("**/ncbi_taxa_node.parquet"):
            taxdump_id = self._taxdump_id_from_path_part(node_file_path.parent.name)
            cached_taxdump_ids.append(taxdump_id)
        return sorted(cached_taxdump_ids)

    def taxon_db(self, taxdump_id: str, threads: int = 1) -> TaxonDB:
        """Access the NCBI Taxonomy database specified by the given taxdump ID.

        Args:
            taxdump_id: Taxdump ID of taxonomy database to access.
            threads: Number of threads to use.

        Returns:
            TaxonDB providing access to the specified taxonomy database.
        """
        taxdump_dir_path = self._get_taxdump_dir_path(taxdump_id)
        return TaxonDB(taxdump_dir_path, threads=threads)
