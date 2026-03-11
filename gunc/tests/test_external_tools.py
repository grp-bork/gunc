from ..external_tools import (
    get_record_count_in_fasta,
    check_diamond_version_correct,
    REQUIRED_DIAMOND_VERSION,
)
from importlib.resources import files as _pkg_files


def resource_filename(package, resource):
    return str(_pkg_files(package).joinpath(resource))


fasta_file = resource_filename("gunc.tests", "test_data/tiny_test.faa.gz")


def test_get_record_count_in_fasta():
    assert get_record_count_in_fasta(fasta_file) == 35


def test_check_diamond_version_correct_matching(monkeypatch):
    from .. import external_tools
    monkeypatch.setattr(external_tools, "check_diamond_version", lambda: REQUIRED_DIAMOND_VERSION)
    assert check_diamond_version_correct() is True


def test_check_diamond_version_correct_wrong_version(monkeypatch):
    from .. import external_tools
    monkeypatch.setattr(external_tools, "check_diamond_version", lambda: "0.0.0")
    assert check_diamond_version_correct() is False


def test_check_diamond_version_correct_diamond_missing(monkeypatch):
    from .. import external_tools

    def raise_error():
        raise FileNotFoundError("diamond not found")
    monkeypatch.setattr(external_tools, "check_diamond_version", raise_error)
    assert check_diamond_version_correct() is False
