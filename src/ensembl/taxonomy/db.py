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
"""Module for working with an NCBI Taxonomy database."""

__all__ = [
    "TaxonDB",
    "TaxonCacheError",
    "UnknownTaxaError",
]

from collections.abc import Iterable
from contextlib import contextmanager, ExitStack
from itertools import islice
from numbers import Integral
import operator
import os
from pathlib import Path
import re
from typing import Any, Iterator, Optional

import duckdb
import numpy as np
import pandas as pd
from numpy.typing import NDArray

from ensembl.taxonomy.utils import group_by_array


_NCBI_TAXA_TABLE_NAMES = ["ncbi_taxa_name", "ncbi_taxa_node"]

_NCBI_TAXON_ID_RE = re.compile(r"(?:ncbitaxon:)?(?P<taxon_id>[0-9]+)")


class TaxonCacheError(Exception):
    """Taxonomy cache exception class."""


class UnknownTaxaError(Exception):
    """Exception class for taxonomy IDs that are not in the given taxonomy database."""


def _execute_multi_value_query(
    conn: duckdb.DuckDBPyConnection,
    base_query: str,
    value_lists: list[tuple[str, Any]],
    conditions: Optional[list[str]] = None,
    final_clauses: Optional[list[str]] = None,
    params: Optional[dict] = None,
) -> dict[str, NDArray | pd.Categorical]:

    conditions = list(conditions) if conditions is not None else []
    final_clauses = list(final_clauses) if final_clauses is not None else []
    params = dict(params) if params is not None else {}

    with ExitStack() as stack:
        join_clauses = []
        for value_list_num, (col_name, value_list) in enumerate(value_lists, start=1):
            if len(value_list) > 400:
                tmp_table_name = f"t{value_list_num}"
                stack.enter_context(_temp_value_table(conn, tmp_table_name, value_list))
                join_clauses.append(f"JOIN {tmp_table_name} ON {tmp_table_name}.column0 = {col_name}")
            else:
                param_name = f"t{value_list_num}"
                conditions.append(f"{col_name} IN ${param_name}")
                params[param_name] = value_list

        conditional_clauses = []
        if conditions:
            conditional_clauses.append(f"WHERE ({conditions[0]})")
            for condition in conditions[1:]:
                conditional_clauses.append(f"AND ({condition})")

        query = " ".join([base_query] + join_clauses + conditional_clauses + final_clauses)
        res_arrays = conn.execute(query, parameters=params).fetchnumpy()

    return res_arrays


def _execute_single_taxon_query(
    conn: duckdb.DuckDBPyConnection,
    query: str,
    fetch_method: str,
    taxon_id: int,
    params: Optional[dict] = None,
) -> Any:

    if not params:
        params = {}

    try:
        parameters = {"taxon_id": taxon_id} | params
        result = conn.execute(query, parameters=parameters)
    except duckdb.ConversionException as exc:  # pylint: disable=c-extension-no-member
        raise ValueError(f"invalid taxon_id: {taxon_id!r}") from exc

    fetch = operator.methodcaller(fetch_method)
    return fetch(result)


@contextmanager
def _temp_value_table(
    conn: duckdb.DuckDBPyConnection, table_name: str, table_values: Iterable[Any]
) -> Iterator[str]:
    _value_arr = np.array(table_values)
    conn.execute(f"CREATE TEMP TABLE {table_name} AS SELECT * FROM _value_arr")
    try:
        yield table_name
    finally:
        conn.execute(f"DROP TABLE IF EXISTS {table_name}")


class TaxonDB:
    """An NCBI Taxonomy database class."""

    def __init__(self, path: os.PathLike | str, threads=1) -> None:

        self._db_dir_path = Path(path)
        self._conn = duckdb.connect(config={"threads": threads})
        self._conn.execute("SET enable_progress_bar = false")
        self._conn.execute(f"SET file_search_path = '{self._db_dir_path}'")

        for table_name in _NCBI_TAXA_TABLE_NAMES:
            try:
                self._conn.execute(f"CREATE TABLE {table_name} AS SELECT * FROM '{table_name}.parquet'")
            except duckdb.IOException as exc:  # pylint: disable=c-extension-no-member
                missing_pqt_err_re = re.compile(r'IO Error: No files found that match the pattern ".+"')
                msg = exc.args[0]
                if missing_pqt_err_re.fullmatch(msg):
                    raise TaxonCacheError(f"taxonomy database not found at '{path}'") from exc

    def __contains__(self, taxon_id: int | str) -> bool:
        # Tests whether this NCBI Taxonomy database contains the given Taxonomy ID.
        int_taxon_id = self._preprocess_taxon_id(taxon_id)
        query = "SELECT $taxon_id IN (SELECT taxon_id FROM ncbi_taxa_node)"
        row = _execute_single_taxon_query(self._conn, query, "fetchone", int_taxon_id)
        return row[0]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self._conn.close()

    def __getattr__(self, name: str) -> Any:

        attr_name_to_query = {
            "_archive_id": "SELECT name FROM ncbi_taxa_name WHERE name_class == 'archive_id' LIMIT 1",
            "_import_date": "SELECT name FROM ncbi_taxa_name WHERE name_class == 'import date' LIMIT 1",
            "_name_classes": "SELECT DISTINCT name_class FROM ncbi_taxa_name ORDER BY name_class",
            "_rank_names": "SELECT DISTINCT rank FROM ncbi_taxa_node ORDER BY rank",
            "_root_id": "SELECT root_id FROM ncbi_taxa_node LIMIT 1",
            "_taxdump_id": "SELECT name FROM ncbi_taxa_name WHERE name_class == 'taxdump_id' LIMIT 1",
        }

        try:
            query = attr_name_to_query[name]
        except KeyError as exc:
            raise AttributeError(f"{type(self)} object has no attribute '{name}'") from exc

        if name in ("_archive_id", "_import_date", "_root_id", "_taxdump_id"):
            res = self._conn.execute(query).fetchone()
            value = res[0] if res is not None else None
        elif name == "_name_classes":
            value = self._conn.execute(query).fetchnumpy()["name_class"].tolist()
        elif name == "_rank_names":
            value = self._conn.execute(query).fetchnumpy()["rank"].tolist()
        else:
            assert name in attr_name_to_query, f"attribute name '{name}' must have a taxonomy database query"

        setattr(self, name, value)
        return value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self._taxdump_id}')"

    def _preprocess_taxon_id(self, taxon_id: int | str) -> int:
        return self._preprocess_taxon_ids([taxon_id])[0]

    def _preprocess_taxon_ids(self, taxon_ids: Iterable[int | str]) -> list[int]:
        int_taxon_ids = []
        for taxon_id in taxon_ids:
            if isinstance(taxon_id, Integral):
                if isinstance(taxon_id, np.integer):
                    int_taxon_id = taxon_id.item()
                else:
                    int_taxon_id = taxon_id
            elif isinstance(taxon_id, str):
                if match := _NCBI_TAXON_ID_RE.fullmatch(taxon_id):
                    int_taxon_id = int(match["taxon_id"])
                else:
                    raise ValueError(f"invalid taxon_id: '{taxon_id}'")
            else:
                raise TypeError(f"taxon_id {taxon_id!r} is of unsupported type {type(taxon_id)}")
            int_taxon_ids.append(int_taxon_id)
        return int_taxon_ids

    def _divergent_lineage_ids(self, leaf_ids: list[int]) -> list[int]:

        root_query = """
            SELECT
                n1.left_index,
                n1.right_index
            FROM
                ncbi_taxa_node n1
            WHERE
                n1.left_index <= $min_left_index
            AND
                n1.right_index >= $max_right_index
            ORDER BY
                n1.left_index DESC
            LIMIT 1
        """
        min_left_index, max_right_index = self._get_leftright_index_extrema(leaf_ids)
        root_params = {"min_left_index": min_left_index, "max_right_index": max_right_index}
        root_left_right_indices = self._conn.execute(root_query, parameters=root_params).fetchone()
        assert (
            root_left_right_indices
        ), f"taxonomy must have left-right indices containing interval ({min_left_index}, {max_right_index})"
        root_left_index, root_right_index = root_left_right_indices

        min_lr_diff_query = """
            SELECT
                MIN(right_index - left_index) AS min_lr_diff
            FROM
                ncbi_taxa_node n1
        """
        leaf_taxon_id_sets = [("n1.taxon_id", leaf_ids)]
        min_lr_diff_arrays = _execute_multi_value_query(self._conn, min_lr_diff_query, leaf_taxon_id_sets)
        min_lr_diff = min_lr_diff_arrays["min_lr_diff"].tolist()[0]

        lineage_query = """
            SELECT
                DISTINCT n2.taxon_id
            FROM
                ncbi_taxa_node n1
            JOIN
                ncbi_taxa_node n2
            ON
                n2.left_index < n1.left_index
            AND
                n2.right_index > n1.right_index
        """

        lineage_conditions = [
            "n2.left_index >= $root_left_index",
            "n2.right_index <= $root_right_index",
            "n2.right_index - n2.left_index > $min_lr_diff",
        ]
        lineage_params = {
            "root_left_index": root_left_index,
            "root_right_index": root_right_index,
            "min_lr_diff": min_lr_diff,
        }
        lineage_arrays = _execute_multi_value_query(
            self._conn,
            lineage_query,
            leaf_taxon_id_sets,
            conditions=lineage_conditions,
            params=lineage_params,
        )

        # We return the union set of the leaf and lineage taxon_ids, in case
        # there are lineal relationships between any of the input taxon_ids.
        return list(set(leaf_ids) | set(lineage_arrays["taxon_id"].tolist()))

    def _get_leftright_index_extrema(self, taxon_ids: list[int]) -> tuple[int, int]:

        query = """
            SELECT
                left_index,
                right_index
            FROM
                ncbi_taxa_node n1
        """

        taxon_id_sets = [("n1.taxon_id", taxon_ids)]
        res_arrays = _execute_multi_value_query(self._conn, query, taxon_id_sets)
        left_index_arr = res_arrays["left_index"]
        if len(left_index_arr) != len(taxon_ids):
            unknown_taxon_ids = [x for x in taxon_ids if x not in self]
            unknown_taxon_id_str = ", ".join(map(str, unknown_taxon_ids))
            raise UnknownTaxaError(
                f"{len(unknown_taxon_ids)} taxon_ids not found in"
                f" taxdump '{self._taxdump_id}': {unknown_taxon_id_str}"
            )

        return left_index_arr.min().item(), res_arrays["right_index"].max().item()

    def _query_taxon_attribute(
        self,
        query_method_name: str,
        taxon_attr_name: str,
        query_taxon_id: int | str,
        convert_method_name: str,
        required: bool = False,
    ):

        int_taxon_id = self._preprocess_taxon_id(query_taxon_id)

        if query_method_name == "query_taxa_attributes":
            dataframes = list(self.query_taxa_attributes([taxon_attr_name], [int_taxon_id]))
        else:
            assert (
                query_method_name == "query_taxa_attribute"
            ), "taxon query method must be 'query_taxa_attribute' or 'query_taxa_attributes'"
            dataframes = list(self.query_taxa_attribute(taxon_attr_name, [int_taxon_id]))

        try:
            df = dataframes[0]
        except IndexError as exc:
            if int_taxon_id not in self:
                raise UnknownTaxaError(
                    f"ncbitaxon:{int_taxon_id} not found in taxonomy database '{self.taxdump_id}'"
                ) from exc
            if required:
                raise RuntimeError(
                    f"{taxon_attr_name} of ncbitaxon:{int_taxon_id}"
                    f" could not be obtained from taxonomy database '{self.taxdump_id}'"
                ) from exc
            taxon_attr_value = None
        else:
            convert = operator.methodcaller(convert_method_name)
            res_arr = df[taxon_attr_name].to_numpy()
            taxon_attr_value = convert(res_arr[0])

        return taxon_attr_value

    def _validate_name_classes(self, name_classes: Iterable[str], label: str = "name class") -> None:
        unk_name_classes = set(name_classes) - set(self._name_classes)
        if unk_name_classes:
            unk_name_class_str = "', '".join(unk_name_classes)
            raise ValueError(f"unknown {label} value(s): '{unk_name_class_str}'")

    def ancestor_ids(self, taxon_id):
        """Get taxonomy IDs of the ancestors of the given taxon.

        Args:
            taxon_id: Get the ancestors of the taxon with this taxonomy ID.

        Returns:
            Taxonomy IDs of the ancestors of the given taxon.
        """
        return self.lineage_ids(taxon_id)[:-1]

    @property
    def archive_id(self) -> str | None:
        """Return date of this NCBI Taxonomy archive."""
        return self._archive_id

    def available_name_classes(self) -> list[str]:
        """Return name classes in this taxonomy database."""
        return self._name_classes

    def available_ranks(self) -> list[str]:
        """Return names of ranks in this taxonomy database."""
        return self._rank_names

    def child_ids(self, taxon_id: int | str) -> list:
        """Return child taxonomy IDs of the given taxon.

        Args:
            taxon_id: Get the children of the taxon with this taxonomy ID.

        Returns:
            Unordered list of taxonomy IDs of the children of the given taxon.
        """
        return self._query_taxon_attribute(
            "query_taxa_attribute", "child_ids", taxon_id, "tolist", required=True
        )

    def common_lineage_ids(self, taxon_ids) -> list[int]:
        """Return taxonomy IDs of common lineage of the given taxa.

        Args:
            taxon_ids: Get the common lineage of the taxa with these taxonomy IDs.

        Returns:
            Taxonomy IDs of the common lineage of the given taxa.
        """
        int_taxon_ids = self._preprocess_taxon_ids(taxon_ids)

        query = """
            SELECT
                taxon_id
            FROM
                ncbi_taxa_node
            WHERE
                left_index <= $min_left_index
            AND
                right_index >= $max_right_index
            ORDER BY
                left_index
        """

        min_left_index, max_right_index = self._get_leftright_index_extrema(int_taxon_ids)
        params = {"min_left_index": min_left_index, "max_right_index": max_right_index}

        res_arrays = self._conn.execute(query, parameters=params).fetchnumpy()
        return res_arrays["taxon_id"].astype(int).tolist()

    def descendant_ids(self, taxon_id: int | str) -> list[int]:
        """Get taxonomy IDs of descendants of the given taxon.

        Args:
            taxon_id: Get the descendants of the taxon with this taxonomy ID.

        Returns:
            Taxonomy IDs of the descendants of the given taxon.
        """
        return self._query_taxon_attribute(
            "query_taxa_attribute", "descendant_ids", taxon_id, "tolist", required=True
        )

    def export_subtree(
        self, *, subtree_root_id: Optional[int | str] = None, leaf_ids: Optional[list[int | str]] = None
    ) -> str:
        """Export taxonomy subtree, as specified by subtree root or leaves.

        Args:
            subtree_root_id: Export a subtree containing this taxon and all its descendants.
            leaf_ids: Export a subtree with these taxa as leaves.

        Returns:
            A Newick string representing the specified subtree.
        """

        if subtree_root_id and leaf_ids:
            raise RuntimeError("TaxonDB.export_subtree() takes a subtree_root_id or leaf_ids, but not both")

        if subtree_root_id:
            int_subtree_root_id = self._preprocess_taxon_id(subtree_root_id)
            tree_taxon_ids = [int_subtree_root_id] + self.descendant_ids(int_subtree_root_id)
        elif leaf_ids:
            int_leaf_ids = self._preprocess_taxon_ids(leaf_ids)
            tree_taxon_ids = self._divergent_lineage_ids(int_leaf_ids)
        else:
            raise RuntimeError("TaxonDB.export_subtree() requires either a subtree_root_id or leaf_ids")

        query = "SELECT n1.taxon_id, n1.parent_id, n1.left_index, n1.right_index FROM ncbi_taxa_node n1"
        tree_taxon_id_sets = [("n1.taxon_id", tree_taxon_ids)]
        res_arrays = _execute_multi_value_query(self._conn, query, tree_taxon_id_sets)

        taxon_id_arr = res_arrays["taxon_id"]
        parent_id_arr = res_arrays["parent_id"]
        left_index_arr = res_arrays["left_index"]
        right_index_arr = res_arrays["right_index"]

        lr_index_arr = np.concatenate((left_index_arr, right_index_arr))
        sort_idxs = np.argsort(lr_index_arr)

        lr_taxon_arr = np.tile(taxon_id_arr, 2)[sort_idxs]
        lr_parent_arr = np.tile(parent_id_arr, 2)[sort_idxs]

        tree_leaf_ids = lr_taxon_arr[np.diff(lr_taxon_arr, prepend=0) == 0]
        lr_leaf_mask = np.isin(lr_taxon_arr, tree_leaf_ids)

        right_index_mask = np.full_like(taxon_id_arr, True, dtype=np.bool_)
        lr_right_mask = np.concatenate((~right_index_mask, right_index_mask))[sort_idxs]

        max_name_length = len(str(lr_taxon_arr.max()))
        max_nwk_part_length = max_name_length + 1  # add 1 for closing paren: ')'
        nwk_arr = np.full_like(lr_taxon_arr, "", dtype=f"<U{max_nwk_part_length}")

        diff_taxon_mask = np.diff(lr_taxon_arr, prepend=lr_taxon_arr[0]) != 0
        same_parent_mask = np.diff(lr_parent_arr, prepend=-1) == 0
        comma_mask = same_parent_mask & diff_taxon_mask
        nwk_arr[comma_mask] = ","

        paren1_mask = ~lr_right_mask & ~lr_leaf_mask
        nwk_arr[paren1_mask] = np.char.add(nwk_arr[paren1_mask], "(")

        paren2_mask = lr_right_mask & ~lr_leaf_mask
        nwk_arr[paren2_mask] = ")"
        nwk_arr[lr_right_mask] = np.char.add(
            nwk_arr[lr_right_mask], lr_taxon_arr[lr_right_mask].astype(nwk_arr.dtype)
        )

        return "".join(nwk_arr) + ";"

    def genbank_common_name(self, taxon_id: int | str) -> str | None:
        """Return GenBank common name of the given taxon.

        Args:
            taxon_id: Get the GenBank common name of the taxon with this taxonomy ID.

        Returns:
            If available, GenBank common name of the given taxon; otherwise returns None.
        """
        return self._query_taxon_attribute(
            "query_taxa_attribute", "genbank common name", taxon_id, "item", required=True
        )

    @property
    def import_date(self) -> str:
        """Return import date/time of this NCBI Taxonomy database."""
        return self._import_date

    def last_common_taxon_id(self, taxon_ids: Iterable[int | str]) -> int:
        """Return Taxonomy ID of the last common taxon of the given taxa."""
        return self.common_lineage_ids(taxon_ids)[-1]

    def lineage_ids(self, taxon_id: int | str) -> list[int]:
        """Return Taxonomy IDs of the lineage of the given taxon.

        Args:
            taxon_id: Get the lineage of the taxon with this Taxonomy ID.

        Returns:
            Taxonomy IDs of the lineage of the given taxon.
        """
        return self._query_taxon_attribute(
            "query_taxa_attribute", "lineage_ids", taxon_id, "tolist", required=True
        )

    def num_ancestors(self, taxon_id: int | str) -> int:
        """Return number of ancestors of the given taxon."""
        return self._query_taxon_attribute(
            "query_taxa_attribute", "num_ancestors", taxon_id, "item", required=True
        )

    def num_descendants(self, taxon_id: int | str) -> int:
        """Return number of descendants of the given taxon."""
        return self._query_taxon_attribute(
            "query_taxa_attribute", "num_descendants", taxon_id, "item", required=True
        )

    def parent_id(self, taxon_id: int | str) -> int | None:
        """Return Taxonomy ID of the parent of the given taxon.

        Args:
            taxon_id: Get the parent of the taxon with this Taxonomy ID.

        Returns:
            Taxonomy ID of the parent of the given taxon.
        """
        parent_id = self._query_taxon_attribute(
            "query_taxa_attributes", "parent_id", taxon_id, "item", required=True
        )
        return parent_id or None

    def query_taxa_attribute(
        self, attr_name: str, taxon_ids: Optional[Iterable[int | str]] = None
    ) -> Iterator[pd.DataFrame]:
        """Query attribute for the given Taxonomy IDs.

        Args:
            attr_name: Name of the attribute to query. In addition to those supported by
                ``query_taxa_attributes``, this method supports attributes ``child_ids``,
                ``descendant_ids``, ``lineage_ids``, ``num_ancestors`` and ``num_descendants``.
            taxon_ids: Query attribute for these Taxonomy IDs. By default,
                the specified attribute is returned for all taxa.

        Yields:
            Dataframe chunk of query results.
        """

        attr_name_to_query = {
            "child_ids": """
                SELECT
                    n1.taxon_id,
                    list(n2.taxon_id) AS child_ids
                FROM
                    ncbi_taxa_node n1
                JOIN
                    ncbi_taxa_node n2
                ON
                    n2.parent_id = n1.taxon_id
            """,
            "descendant_ids": """
                SELECT
                    n1.taxon_id,
                    list(n2.taxon_id) AS descendant_ids
                FROM
                    ncbi_taxa_node n1
                JOIN
                    ncbi_taxa_node n2
                ON
                    n2.left_index > n1.left_index
                AND
                    n2.right_index < n1.right_index
            """,
            "lineage_ids": """
                SELECT
                    n1.taxon_id,
                    list(n2.taxon_id ORDER BY n2.left_index) AS lineage_ids
                FROM
                    ncbi_taxa_node n1
                JOIN
                    ncbi_taxa_node n2
                ON
                    n2.left_index <= n1.left_index
                AND
                    n2.right_index >= n1.right_index
            """,
            "num_ancestors": """
                SELECT
                    n1.taxon_id,
                    count(*) - 1 AS num_ancestors
                FROM
                    ncbi_taxa_node n1
                JOIN
                    ncbi_taxa_node n2
                ON
                    n2.left_index <= n1.left_index
                AND
                    n2.right_index >= n1.right_index
            """,
            "num_descendants": """
                SELECT
                    n1.taxon_id,
                    CAST((n1.right_index - n1.left_index - 1) / 2 AS INTEGER) AS num_descendants
                FROM
                    ncbi_taxa_node n1
            """,
        }

        try:
            base_query = attr_name_to_query[attr_name]
        except KeyError:
            yield from self.query_taxa_attributes([attr_name], taxon_ids=taxon_ids)
        except TypeError as exc:
            raise TypeError(
                f"taxon attribute name {attr_name!r} is of unsupported type {type(attr_name)}"
            ) from exc
        else:
            if attr_name == "num_descendants":
                final_clauses = []
            else:
                final_clauses = ["GROUP BY n1.taxon_id"]

            if taxon_ids:
                int_taxon_ids = self._preprocess_taxon_ids(taxon_ids)
            else:
                int_taxon_ids = self.taxon_ids()

            taxon_id_batch_size = 512
            taxon_id_iterator = iter(int_taxon_ids)
            # From Python 3.12 onwards, we can use 'itertools.batched' here.
            while batch_taxon_ids := tuple(islice(taxon_id_iterator, taxon_id_batch_size)):
                params = {"t1": batch_taxon_ids}
                conditional_clause = "WHERE (n1.taxon_id IN $t1)"
                query = " ".join([base_query, conditional_clause] + final_clauses)
                df_chunk = self._conn.execute(query, parameters=params).df()
                if df_chunk.empty:
                    continue
                yield df_chunk

    def query_taxa_attributes(self, attr_names, taxon_ids=None):
        """Query attributes for the given Taxonomy IDs.

        Args:
            attr_names: Iterable of attributes to query. Supported attributes include ``parent_id``,
                ``rank``, ``genbank_hidden_flag`` and any available taxon name class.
            taxon_ids: Query attributes for these Taxonomy IDs. By default,
                the specified attributes are returned for all taxa.

        Yields:
            Dataframe chunk of query results.
        """
        attr_names = list(attr_names)

        if not attr_names:
            raise ValueError("no taxon attribute name specified; please specify at least one")

        known_node_attr_names = [
            "taxon_id",
            "parent_id",
            "rank",
            "genbank_hidden_flag",
        ]

        # We want taxon_id to be in the first column.
        if "taxon_id" in attr_names:
            attr_names.remove("taxon_id")
        attr_names.insert(0, "taxon_id")

        name_classes = []
        node_attr_names = []
        for attr_name in attr_names:
            if attr_name in known_node_attr_names:
                node_attr_names.append(attr_name)
            else:
                name_classes.append(attr_name)
        self._validate_name_classes(name_classes, label="attribute name")

        if len(node_attr_names) > 1:
            table_clauses = ["FROM ncbi_taxa_node"]
            if name_classes:
                table_clauses.append("JOIN ncbi_taxa_name USING (taxon_id)")
        else:
            table_clauses = ["FROM ncbi_taxa_name"]

        with ExitStack() as stack:

            params = {}
            conditions = []
            if taxon_ids:
                int_taxon_ids = self._preprocess_taxon_ids(taxon_ids)
                if len(taxon_ids) > 512:
                    tmp_table_name = "t1"
                    stack.enter_context(_temp_value_table(self._conn, tmp_table_name, int_taxon_ids))
                    table_clauses.append(f"JOIN {tmp_table_name} ON {tmp_table_name}.column0 = taxon_id")
                else:
                    param_name = "t1"
                    conditions.append(f"taxon_id IN ${param_name}")
                    params[param_name] = int_taxon_ids

            if name_classes:
                query_col_names = node_attr_names + ["name", "name_class"]
                conditions.append("name_class IN $name_classes")
                params["name_classes"] = name_classes
            else:
                query_col_names = node_attr_names
            select_clause = f"SELECT {', '.join(query_col_names)}"

            conditional_clauses = []
            if conditions:
                conditional_clauses.append(f"WHERE ({conditions[0]})")
                for condition in conditions[1:]:
                    conditional_clauses.append(f"AND ({condition})")

            query = " ".join([select_clause] + table_clauses + conditional_clauses)

            if name_classes:
                query = f"""
                    PIVOT ({query})
                    ON name_class IN ('{"', '".join(name_classes)}')
                    USING list(name)
                    GROUP BY {', '.join(node_attr_names)}
                """

            cursor = self._conn.execute(query, parameters=params)
            while True:
                df_chunk = cursor.fetch_df_chunk()
                if df_chunk.empty:
                    break
                yield df_chunk

    def query_taxa_by_name(
        self,
        pattern: str,
        name_classes: Optional[Iterable[str]] = None,
        case_sensitive: bool = False,
        full_match: bool = False,
    ) -> list[int]:
        """Query taxon names with the given pattern.

        Args:
            pattern: Query string containing a taxon name or RE2 regular expression.
            name_classes: Name classes to search with the given query string.
            case_sensitive: Search in a case-sensitive manner.
            full_match: Require a full match of the taxon name.

        Returns:
            Taxonomy IDs matching the search criteria.
        """
        regexp_func_name = "regexp_full_match" if full_match else "regexp_matches"
        options = "c" if case_sensitive else "i"

        query_parts = [
            f"""
            SELECT taxon_id
            FROM ncbi_taxa_name
            WHERE {regexp_func_name}(name, '{pattern}', '{options}')
        """
        ]

        params = {}
        if name_classes:
            name_classes = list(name_classes)
            self._validate_name_classes(name_classes)
            query_parts.append("AND name_class IN $name_classes")
            params["name_classes"] = name_classes

        query = " ".join(query_parts)

        res_arrays = self._conn.execute(query, parameters=params).fetchnumpy()
        return res_arrays["taxon_id"].astype(int).tolist()

    def rank(self, taxon_id: int | str) -> str:
        """Return rank of the given taxon."""
        return self._query_taxon_attribute("query_taxa_attributes", "rank", taxon_id, "item", required=True)

    @property
    def root_id(self) -> int:
        """Return root ID of this taxonomy."""
        return self._root_id

    def scientific_name(self, taxon_id: int | str) -> str:
        """Return scientific name of the given taxon.

        Args:
            taxon_id: Get the scientific name of the taxon with this Taxonomy ID.

        Returns:
            Scientific name of the given taxon.
        """
        return self._query_taxon_attribute(
            "query_taxa_attributes", "scientific name", taxon_id, "item", required=True
        )

    @property
    def taxdump_id(self) -> str:
        """Return NCBI taxdump ID."""
        return self._taxdump_id

    def taxon_ids(self) -> list:
        """Return all taxon_ids in this Taxonomy."""
        query = "SELECT taxon_id FROM ncbi_taxa_node"
        return self._conn.execute(query).fetchnumpy()["taxon_id"].tolist()

    def taxon_names(
        self, taxon_id: int | str, name_classes: Optional[Iterable[str]] = None
    ) -> dict[str, list[str]]:
        """Return names of the given taxon.

        Args:
            taxon_id: Get the names of the taxon with this Taxonomy ID.
            name_classes: Name classes for which names should be obtained.

        Returns:
            Mapping containing the names of the given taxon, grouped by name class.
        """
        int_taxon_id = self._preprocess_taxon_id(taxon_id)

        query_parts = ["SELECT name_class, name FROM ncbi_taxa_name WHERE taxon_id = $taxon_id"]
        sci_name_is_sentinel = False
        params = {}
        if name_classes:
            name_classes = list(name_classes)
            self._validate_name_classes(name_classes)

            if "scientific name" not in name_classes:
                # Because all taxa are expected to have a scientific name, we use it as a sentinel value.
                # When we later check the query result, if we do not find a scientific name, we assume the
                # taxonomy ID is not in the current taxonomy database, and we flag it as an unknown taxon.
                name_classes.append("scientific name")
                sci_name_is_sentinel = True

            query_parts.append("AND name_class IN $name_classes")
            params["name_classes"] = name_classes
        query = " ".join(query_parts)

        res_arrays = _execute_single_taxon_query(self._conn, query, "fetchnumpy", int_taxon_id, params=params)
        taxon_names_by_class = group_by_array(res_arrays["name"], res_arrays["name_class"])

        if "scientific name" not in taxon_names_by_class:
            raise UnknownTaxaError(
                f"ncbitaxon:{int_taxon_id} not found" f" in taxonomy database '{self.taxdump_id}'"
            )

        if sci_name_is_sentinel:
            taxon_names_by_class.pop("scientific name")

        return taxon_names_by_class
