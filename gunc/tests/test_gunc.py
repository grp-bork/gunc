import pytest
import logging
from ..gunc import *
from .._version import get_versions
from pkg_resources import resource_filename

gunc_gene_counts = resource_filename(__name__, "test_data/tiny_test.gene_counts.json")
test_data_dir = resource_filename(__name__, "test_data/")


def test_get_genecount_from_gunc_output():
    assert get_genecount_from_gunc_output(gunc_gene_counts, "tiny_test.faa") == 35


def test_start_checks():
    logger = logging.getLogger()
    with pytest.raises(SystemExit):
        start_checks()


def test_parse_args():
    with pytest.raises(SystemExit):
        parse_args(["-h"])
    with pytest.raises(SystemExit):
        parse_args(["-f"])
    parser = parse_args(["run", "-f", "test_path"])
    assert parser.input_file == "test_path"
    assert parser.sensitive is False


def test_get_files_in_dir_with_suffix():
    assert len(get_files_in_dir_with_suffix(test_data_dir, ".json")) == 1
