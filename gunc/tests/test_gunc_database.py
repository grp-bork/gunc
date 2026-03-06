import pytest
from ..gunc_database import md5sum_file, download_file, get_md5_from_url
from importlib.resources import files as _pkg_files

def resource_filename(package, resource):
    return str(_pkg_files(package).joinpath(resource))


diamond_output = resource_filename("gunc.tests", "test_data/tiny_test.diamond.out")
base_url = "https://swifter.embl.de/~fullam/gunc/"
file_name = "gunc_db_progenomes2.1.dmnd.gz"
gz_file_url = f"{base_url}{file_name}"
md5_file_url = f"{base_url}{file_name}.md5"


def test_md5sum_file():
    assert md5sum_file(diamond_output) == "3fb9214bfa31f8dda65011efcba56a3a"


@pytest.mark.integration
def test_download_file(tmp_path):
    tmp_outfile = tmp_path / "md5file"
    download_file(md5_file_url, tmp_outfile)
    assert (
        tmp_outfile.read_text()
        == "bc93a855e0760aad5c4e5f2d0e26da46  gunc_db_progenomes2.1.dmnd.gz\n"
    )


@pytest.mark.integration
def test_get_md5_from_url():
    assert get_md5_from_url(gz_file_url) == "bc93a855e0760aad5c4e5f2d0e26da46"


