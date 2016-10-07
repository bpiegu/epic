import pytest

from collections import namedtuple

import pandas as pd

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"

MockNamespace = namedtuple(
    "MockNamespace",
    ["number_cores", "genome", "keep_duplicates", "window_size",
     "fragment_size", "paired_end", "gaps_allowed",
     "false_discovery_rate_cutoff", "effective_genome_length", "store_matrix",
     "bigwig", "bed", "sum_bigwig",
     "treatment", "control"])

egs = 2290813547.4  # this is the effective genome size used by the original sicer for hg19


@pytest.fixture(scope="session")
def args_200_fast():
    return MockNamespace(25, "hg19", False, 200, 150, False, 3, 1, egs, False, False, False,
                         False, ["examples/test.bed"],
                         ["examples/control.bed"])


@pytest.fixture(scope="session")
def args_200():
    return MockNamespace(1, "hg19", False, 200, 150, False, 3, 0.05, egs,
                         False, False, False, False,
                         ["examples/test.bed"], ["examples/control.bed"])


@pytest.fixture(scope="session")
def args_50():
    return MockNamespace(1, "hg19", False, 50, 150, False, 3, 0.05, egs, False, False, False,
                         False, ["examples/test.bed"],
                         ["examples/control.bed"])


@pytest.fixture(scope="session")
def expected_result_example_input():
    return pd.read_table("examples/expected_results.csv", sep=" ", skiprows=1)


@pytest.fixture(scope="session")
def epic_overlap_files():
    return ["/local/home/endrebak/code/epic/examples/epic-overlaps/0h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/3h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/6h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/9h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/12h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/15h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/18h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/21h_H3K4me3.regions",
"/local/home/endrebak/code/epic/examples/epic-overlaps/24h_H3K4me3.regions"]



@pytest.fixture(scope="session")
def epic_overlap_intermediate_nucleotide_matrixes():
    return ["/local/home/endrebak/code/epic/examples/epic-overlaps/0h_H3K4me3.matrix",
            "/local/home/endrebak/code/epic/examples/epic-overlaps/3h_H3K4me3.matrix",
            "/local/home/endrebak/code/epic/examples/epic-overlaps/6h_H3K4me3.matrix"]


@pytest.fixture(scope="session")
def epic_overlap_intermediate_region_matrixes():
    return ["/local/home/endrebak/code/epic/examples/epic-overlaps/0h_H3K4me3_region.matrix",
            "/local/home/endrebak/code/epic/examples/epic-overlaps/3h_H3K4me3_region.matrix",
            "/local/home/endrebak/code/epic/examples/epic-overlaps/6h_H3K4me3_region.matrix"]
