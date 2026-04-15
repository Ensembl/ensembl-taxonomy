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
"""Taxonomy utilities module."""

__all__ = [
    "group_by_array",
    "inclusive_range_between_values",
    "NcbiTaxdumpDialect",
]

from collections.abc import Iterable, Sequence
import csv
from typing import Any

import numpy as np
from numpy.typing import NDArray


class NcbiTaxdumpDialect(csv.Dialect):
    """CSV dialect for NCBI taxdump files."""

    delimiter = "\t"
    doublequote = False
    escapechar = None
    lineterminator = "\n"
    quotechar = '"'
    quoting = csv.QUOTE_NONE
    skipinitialspace = False
    strict = True


def group_by_array(data_arr: NDArray, group_arr: NDArray, sort_values: bool = False) -> dict:
    """Group the data array by the values in the given group array.

    Adapted under a CC BY-SA 4.0 license ( https://creativecommons.org/licenses/by-sa/4.0/ )
    from solution by Rubikkon (2021-07-30) Group values from one numpy array based on values
    from another. Stack Overflow. https://stackoverflow.com/a/68596754

    Args:
        data_arr: A 1D NumPy array containing values to be grouped.
        group_arr: A NumPy array of the same shape as ``data_arr``,
            in which each element contains the grouping key of the
            corresponding element of ``data_arr``.
        sort_values: True if group values should be sorted,
            False otherwise.

    Returns:
        A dictionary mapping group key to values.
    """
    group_sort_idxs = np.argsort(group_arr)
    group_arr = group_arr[group_sort_idxs]
    data_arr = data_arr[group_sort_idxs]
    group_ids, group_start_idxs = np.unique(group_arr, return_index=True)
    data_by_group = np.split(data_arr, group_start_idxs[1:])

    data_by_group_id = {}
    for group_id, group_values in zip(group_ids.tolist(), data_by_group):
        group_values = group_values.tolist()
        if sort_values:
            group_values.sort()
        data_by_group_id[group_id] = group_values

    return data_by_group_id


def inclusive_range_between_values(
    sequence: Sequence[Any], first_value: Any, final_value: Any
) -> Iterable[Any]:
    """Return the inclusive range between the specified values in the given sequence.

    Args:
        sequence: A sequence of values.
        first_value: Value with which the inclusive range will start.
        final_value: Value with which the inclusive range will end.

    Returns:
        A range which can be used to iterate over the given sequence,
        with the first and final values being those specified.
    """
    first_idx = sequence.index(first_value)
    final_idx = sequence.index(final_value)
    if first_idx <= final_idx:
        inclusive_range = range(first_idx, final_idx + 1, 1)
    else:
        inclusive_range = range(first_idx, final_idx - 1, -1)
    return inclusive_range
