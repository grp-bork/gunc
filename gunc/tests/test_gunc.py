import pytest
from ..gunc import *
from .._version import get_versions
from pkg_resources import resource_filename

gunc_gene_counts = resource_filename(__name__,
                                     'test_data/tiny_test.gene_counts.json')


def test_get_genecount_from_gunc_output():
    assert get_genecount_from_gunc_output(gunc_gene_counts,
                                          'tiny_test.faa') == 35


def test_start_checks():
    with pytest.raises(SystemExit):
        start_checks()
