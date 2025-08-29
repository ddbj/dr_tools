from pathlib import Path
from typing import Dict

import pytest

from dr_tools.MSS import MSS

HERE = Path(__file__).parent.resolve()
REPO_ROOT = HERE.parent
EXAMPLES_DIR = REPO_ROOT.joinpath("examples")
EG_COMPLETE = {
    "ann": EXAMPLES_DIR.joinpath("complete_genome.ann"),
    "seq": EXAMPLES_DIR.joinpath("complete_genome.fa"),
    "v1_json": EXAMPLES_DIR.joinpath("complete_genome.json"),
    "v2_json": EXAMPLES_DIR.joinpath("complete_genome.v2.json"),
}
EG_VRL = {
    "ann": EXAMPLES_DIR.joinpath("vrl_result.ann"),
    "seq": EXAMPLES_DIR.joinpath("vrl_result.fa"),
    "v1_json": EXAMPLES_DIR.joinpath("vrl_result.json"),
    "v2_json": EXAMPLES_DIR.joinpath("vrl_result.v2.json"),
}


@pytest.fixture
def eg_complete() -> Dict[str, Path]:
    return EG_COMPLETE


@pytest.fixture
def eg_vrl() -> Dict[str, Path]:
    return EG_VRL


@pytest.fixture
def mss_complete() -> MSS:
    return MSS.parse(EG_COMPLETE["ann"], EG_COMPLETE["seq"])


@pytest.fixture
def mss_vrl() -> MSS:
    return MSS.parse(EG_VRL["ann"], EG_VRL["seq"])
