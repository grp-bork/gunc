#!/usr/bin/env python3
import os
import sys
import gzip
import shutil
import hashlib
import logging
import requests

logger = logging.getLogger(__name__)


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


def download_file(file_url, out_file):
    """Download a file to disk

    Streams a file from URL to disk.

    Arguments:
        file_url (str): URL of file to download
        out_file (str): Target file path
    """
    try:
        with requests.get(file_url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(out_file, "wb") as f:
                shutil.copyfileobj(r.raw, f)
    except requests.exceptions.ConnectionError:
        logger.error("Could not connect to server. Check your internet connection.")
        sys.exit(1)
    except requests.exceptions.Timeout:
        logger.error(f"Connection timed out while downloading {file_url}.")
        sys.exit(1)
    except requests.exceptions.HTTPError as e:
        logger.error(f"HTTP error downloading {file_url}: {e}")
        sys.exit(1)
    except requests.exceptions.RequestException as e:
        logger.error(f"Download failed: {e}")
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


def get_md5_from_url(file_url):
    """Read md5sum file from url

    Assumes the md5 file is the file url with .md5 appended

    Arguments:
        file_url (str): URL of file

    Returns:
        str: md5sum
    """
    try:
        r = requests.get(f"{file_url}.md5", timeout=30)
        r.raise_for_status()
        parts = r.text.split(" ")
        if not parts or not parts[0]:
            logger.error("Unexpected MD5 response from server.")
            sys.exit(1)
        return parts[0]
    except requests.exceptions.RequestException as e:
        logger.error(f"Could not retrieve MD5 checksum: {e}")
        sys.exit(1)


def check_md5(file_url, file_path):
    """Check md5 and remove file if incorrect.

    Arguments:
        file_url (str): URL of file
        file_path (str): Path of file
    """
    expected_md5 = get_md5_from_url(file_url)
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
    base_url = "https://black.embl.de/~fullam/gunc/test_data/"
    files = ["ci_test.dmnd", "chimeric.faa", "clean.faa", "genome2taxonomy.tsv"]

    if not os.path.isdir(base_dir):
        logger.error(f"Output directory {base_dir} doesn't exist.")
        sys.exit(1)

    for file_name in files:
        url = f"{base_url}{file_name}"
        out_path = os.path.join(base_dir, file_name)
        logger.info(f"Downloading {file_name}...")
        download_file(url, out_path)
        check_md5(url, out_path)

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

    base_url = "https://black.embl.de/~fullam/gunc/"
    if db == "progenomes_2.1":
        file_name = "gunc_db_progenomes2.1.dmnd.gz"
    elif db == "progenomes_3":
        file_name = "gunc_db_progenomes3.dmnd.gz"
    elif db == "gtdb_95":
        file_name = "gunc_db_gtdb95.dmnd.gz"
    elif db == "gtdb_214":
        file_name = "gunc_db_gtdb214.dmnd.gz"
    else:
        logger.error(f"DB {db} unknown. Allowed: progenomes_2.1, progenomes_3, gtdb_95, gtdb_214, test_data")
        sys.exit(1)
    gz_file_url = f"{base_url}{file_name}"
    gz_file_path = os.path.join(base_dir, file_name)
    file_url = gz_file_url.replace(".gz", "")
    out_file = gz_file_path.replace(".gz", "")

    if not os.path.isdir(base_dir):
        logger.error(f"Output directory {base_dir} doesn't exist.")
        sys.exit(1)

    logger.info("DB downloading...")

    download_file(gz_file_url, gz_file_path)

    logger.info("DB download successful.")
    logger.info("Computing DB md5sum...")

    check_md5(gz_file_url, gz_file_path)

    logger.info("md5sum check successful.")
    logger.info("Uncompressing file...")

    decompress_gzip_file(gz_file_path, out_file)

    logger.info("Decompression complete.")
    logger.info("Computing DB md5sum...")

    check_md5(file_url, out_file)
    os.unlink(gz_file_path)

    logger.info("md5sum check successful.")
    logger.info("DB download successful.")
    logger.info(f"DB saved to: {out_file}. Use with: gunc run -r {out_file}")
