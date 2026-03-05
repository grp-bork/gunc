from ..external_tools import *
from importlib.resources import files as _pkg_files

def resource_filename(package, resource):
    return str(_pkg_files(package).joinpath(resource))

fasta_file = resource_filename(__name__, "test_data/tiny_test.faa.gz")


def test_get_record_count_in_fasta():
    assert get_record_count_in_fasta(fasta_file) == 35
