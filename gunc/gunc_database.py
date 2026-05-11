#!/usr/bin/env python3
import os
import sys
import gzip
import time
import shutil
import hashlib
import logging
import requests

logger = logging.getLogger(__name__)


EMBL_BASE_URL = "https://black.embl.de/~fullam/gunc/"

DB_DOWNLOADS = {
    "progenomes_2.1": {
        "urls": [
            "https://zenodo.org/records/19631736/files/gunc_db_progenomes2.1.dmnd.gz",
            f"{EMBL_BASE_URL}gunc_db_progenomes2.1.dmnd.gz",
        ],
        "file_name": "gunc_db_progenomes2.1.dmnd.gz",
        "gz_md5": "bc93a855e0760aad5c4e5f2d0e26da46",
        "dmnd_md5": "447c9330056b02f29f30fe81fe4af4eb",
    },
    "progenomes_3": {
        "urls": [
            "https://zenodo.org/records/19631883/files/gunc_db_progenomes3.dmnd.gz",
            f"{EMBL_BASE_URL}gunc_db_progenomes3.dmnd.gz",
        ],
        "file_name": "gunc_db_progenomes3.dmnd.gz",
        "gz_md5": "a26b3be496017f291eece27740d0f37f",
        "dmnd_md5": "b1a9347d219a632e9ecf4652dee9bdd1",
    },
    "gtdb_95": {
        "urls": [
            "https://zenodo.org/records/19631804/files/gunc_db_gtdb95.dmnd.gz",
            f"{EMBL_BASE_URL}gunc_db_gtdb95.dmnd.gz",
        ],
        "file_name": "gunc_db_gtdb95.dmnd.gz",
        "gz_md5": "14ed95a2fb1360925e28b7b55df14574",
        "dmnd_md5": "3b502dc047efaafb9831a5eaec7617fd",
    },
    "gtdb_214": {
        "urls": [
            "https://zenodo.org/records/19632326/files/gunc_db_gtdb214.dmnd.gz",
            f"{EMBL_BASE_URL}gunc_db_gtdb214.dmnd.gz",
        ],
        "file_name": "gunc_db_gtdb214.dmnd.gz",
        "gz_md5": "3933007d83a7a85e4672295ee1f9d91f",
        "dmnd_md5": "ac4d6304cab3d62a396703eeb039d3b7",
    },
}

TEST_DATA_BASE_URLS = [
    "https://zenodo.org/records/19631420/files/",
    f"{EMBL_BASE_URL}test_data/",
]
TEST_DATA_FILES = {
    "ci_test.dmnd": "1b040a71b351f4767d1994ca5a6f1a54",
    "chimeric.faa": "4137809237c834348017170e4ee6eb1f",
    "clean.faa": "a28f35b61a5d2adad21ee22e68765405",
    "genome2taxonomy.tsv": "7d074a2063907f1d1c9613043ed0c68a",
}


def md5sum_file(file):
    """Computes MD5 sum of file.

    Arguments:
        file (str): Path of file to md5sum

    Returns:
        str: md5sum
    """
    block_size = 8192
    m = hashlib.md5()
    with open(file, "rb") as f:
        while True:
            data = f.read(block_size)
            if not data:
                return m.hexdigest()
            m.update(data)


def download_file(file_url, out_file, attempts_per_url=3, retry_delay=5):
    """Download a file to disk

    Streams a file from URL to disk. Accepts a single URL or a list of URLs;
    each URL is tried up to attempts_per_url times before moving on to the
    next, and the first successful download wins.

    Arguments:
        file_url (str | list[str]): URL or list of URLs to try
        out_file (str): Target file path
        attempts_per_url (int): Number of attempts per URL before falling back
        retry_delay (int): Seconds to wait between retries of the same URL
    """
    urls = [file_url] if isinstance(file_url, str) else list(file_url)
    last_error = None
    for i, url in enumerate(urls):
        for attempt in range(1, attempts_per_url + 1):
            try:
                with requests.get(url, stream=True, timeout=60) as r:
                    r.raise_for_status()
                    with open(out_file, "wb") as f:
                        shutil.copyfileobj(r.raw, f)
                return
            except requests.exceptions.RequestException as e:
                last_error = e
                if attempt < attempts_per_url:
                    logger.warning(
                        f"Download from {url} failed ({e}). "
                        f"Retrying ({attempt}/{attempts_per_url - 1})..."
                    )
                    time.sleep(retry_delay)
                elif i < len(urls) - 1:
                    logger.warning(f"Download from {url} failed ({e}). Trying fallback...")
    logger.error(f"Download failed: {last_error}")
    sys.exit(1)


def decompress_gzip_file(gz_file, out_file):
    """Decompress a gzip file.

    Uncompressed given gzip file to out_file.

    Arguments:
        gz_file (str): Path of gzip file
        out_file (str): Path of target uncompressed out file
    """
    try:
        with gzip.open(gz_file, "rb") as f_in:
            with open(out_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        if os.path.exists(out_file):
            os.unlink(out_file)
        logger.error(f"Decompression of {gz_file} failed: {e}")
        sys.exit(1)


def check_md5(file_path, expected_md5):
    """Check md5 and remove file if incorrect.

    Arguments:
        file_path (str): Path of file
        expected_md5 (str): Expected md5sum
    """
    downloaded_md5 = md5sum_file(file_path)
    if downloaded_md5 != expected_md5:
        os.unlink(file_path)
        logger.error(f"MD5 check failed for {file_path}. File removed. Try downloading again.")
        sys.exit(1)


def get_test_data(base_dir):
    """Download GUNC test data.

    Downloads a minimal diamond database and two test genomes (chimeric and
    clean) that can be used to verify a GUNC installation end-to-end.

    Arguments:
        base_dir (str): Path of output directory
    """
    if not os.path.isdir(base_dir):
        logger.error(f"Output directory {base_dir} doesn't exist.")
        sys.exit(1)

    for file_name, expected_md5 in TEST_DATA_FILES.items():
        urls = [f"{base}{file_name}" for base in TEST_DATA_BASE_URLS]
        out_path = os.path.join(base_dir, file_name)
        logger.info(f"Downloading {file_name}...")
        download_file(urls, out_path)
        check_md5(out_path, expected_md5)

    logger.info("Test data downloaded successfully.")
    logger.info(f"Files saved to: {base_dir}")
    print("\nTo verify your installation, run:")
    print(
        f"  mkdir gunc_test_out \n"
        f"  gunc run --gene_calls \\\n"
        f"    --input_dir {base_dir} \\\n"
        f"    --file_suffix .faa \\\n"
        f"    --db_file {base_dir}/ci_test.dmnd \\\n"
        f"    --custom_genome2taxonomy {base_dir}/genome2taxonomy.tsv \\\n"
        f"    --out_dir ./gunc_test_out"
        f"  \n"
    )
    print("Expected: chimeric -> pass.GUNC=False, clean -> pass.GUNC=True\n\n")


def get_db(base_dir, db="progenomes_2.1"):
    """Download GUNC DB.

    Arguments:
        base_dir (str): Path of output directory
        db (str): Which db to download. Allowed: progenomes_2.1, progenomes_3,
                  gtdb_95, gtdb_214, test_data
    """
    if db == "test_data":
        get_test_data(base_dir)
        return

    if db not in DB_DOWNLOADS:
        logger.error(f"DB {db} unknown. Allowed: progenomes_2.1, progenomes_3, gtdb_95, gtdb_214, test_data")
        sys.exit(1)
    entry = DB_DOWNLOADS[db]
    gz_file_urls = entry["urls"]
    gz_file_path = os.path.join(base_dir, entry["file_name"])
    out_file = gz_file_path[:-3]  # strip .gz

    if not os.path.isdir(base_dir):
        logger.error(f"Output directory {base_dir} doesn't exist.")
        sys.exit(1)

    logger.info("DB downloading...")

    download_file(gz_file_urls, gz_file_path)

    logger.info("DB download successful.")
    logger.info("Computing DB md5sum...")

    check_md5(gz_file_path, entry["gz_md5"])

    logger.info("md5sum check successful.")
    logger.info("Uncompressing file...")

    decompress_gzip_file(gz_file_path, out_file)

    logger.info("Decompression complete.")
    logger.info("Computing DB md5sum...")

    check_md5(out_file, entry["dmnd_md5"])
    os.unlink(gz_file_path)

    logger.info("md5sum check successful.")
    logger.info("DB download successful.")
    logger.info(f"DB saved to: {out_file}. Use with: gunc run -r {out_file}")
