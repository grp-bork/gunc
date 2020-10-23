from ..gunc_database import *
from pkg_resources import resource_filename


diamond_output = resource_filename(__name__,
                                   'test_data/tiny_test.diamond.out')


def test_md5sum_file():
    assert md5sum_file(diamond_output) == '3fb9214bfa31f8dda65011efcba56a3a'
