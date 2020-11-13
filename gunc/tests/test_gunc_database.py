from ..gunc_database import *
from pkg_resources import resource_filename


diamond_output = resource_filename(__name__,
                                   'test_data/tiny_test.diamond.out')
base_url = 'https://swifter.embl.de/~fullam/gunc/'
file_name = 'gunc_db_2.0.4.dmnd.gz'
gz_file_url = f'{base_url}{file_name}'
md5_file_url = f'{base_url}{file_name}.md5'


def test_md5sum_file():
    assert md5sum_file(diamond_output) == '3fb9214bfa31f8dda65011efcba56a3a'


def test_download_file(tmp_path):
    tmp_outfile = tmp_path / 'md5file'
    download_file(md5_file_url, tmp_outfile)
    assert tmp_outfile.read_text() == 'bc93a855e0760aad5c4e5f2d0e26da46  gunc_db_2.0.4.dmnd.gz\n'


def test_get_md5_from_url():
    assert get_md5_from_url(gz_file_url) == 'bc93a855e0760aad5c4e5f2d0e26da46'


def test_check_md5(mocker):
    mocker.patch('gunc.gunc_database.get_md5_from_url', return_value='3fb9214bfa31f8dda65011efcba56a3a')
    check_md5('x', diamond_output)
