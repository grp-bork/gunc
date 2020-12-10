from ..external_tools import *
from pkg_resources import resource_filename

fasta_file = resource_filename(__name__,
                               'test_data/tiny_test.faa.gz')


def test_get_record_count_in_fasta():
    assert get_record_count_in_fasta(fasta_file) == 35
