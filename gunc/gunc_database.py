#!/usr/bin/env python3
import os
import sys
import gzip
import shutil
import hashlib
import logging
import requests

logger = logging.getLogger("gunc_database.py")


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
        sys.exit("[ERROR] Could not connect to server. Check your internet connection.")
    except requests.exceptions.Timeout:
        sys.exit(f"[ERROR] Connection timed out while downloading {file_url}.")
    except requests.exceptions.HTTPError as e:
        sys.exit(f"[ERROR] HTTP error downloading {file_url}: {e}")
    except requests.exceptions.RequestException as e:
        sys.exit(f"[ERROR] Download failed: {e}")


def decompress_gzip_file(gz_file, out_file):
    """Decompress a gzip file.

    Uncompressed given gzip file to out_file.

    Arguments:
        gz_file (str): Path of gzip file
        out_file (str): Path of target uncompressed out file
    """
    with gzip.open(gz_file, "rb") as f_in:
        with open(out_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


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
            sys.exit("[ERROR] Unexpected MD5 response from server.")
        return parts[0]
    except requests.exceptions.RequestException as e:
        sys.exit(f"[ERROR] Could not retrieve MD5 checksum: {e}")


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


def get_db(base_dir, db="progenomes_2.1"):
    """Download GUNC DB.

    Arguments:
        base_dir (str): Path of output directory
        db (str): Which db to download. Allowed: progenomes_2.1, progenomes_3, gtdb_95, gtdb_214
    """
    base_url = "https://swifter.embl.de/~fullam/gunc/"
    if db == "progenomes_2.1":
        file_name = "gunc_db_progenomes2.1.dmnd.gz"
    elif db == "progenomes_3":
        file_name = "gunc_db_progenomes3.dmnd.gz"
    elif db == "gtdb_95":
        file_name = "gunc_db_gtdb95.dmnd.gz"
    elif db == "gtdb_214":
        file_name = "gunc_db_gtdb214.dmnd.gz"
    else:
        sys.exit(f"[ERROR] DB {db} unknown. Allowed: progenomes_2.1, progenomes_3, gtdb_95, gtdb_214")
    gz_file_url = f"{base_url}{file_name}"
    gz_file_path = os.path.join(base_dir, file_name)
    file_url = gz_file_url.replace(".gz", "")
    out_file = gz_file_path.replace(".gz", "")

    if not os.path.isdir(base_dir):
        sys.exit(f"[ERROR] Output directory {base_dir} doesn't exist.")

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
