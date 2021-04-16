#!/usr/bin/env python3
import os
import sys
import gzip
import shutil
import hashlib
import requests


def md5sum_file(file):
    """Computes MD5 sum of file.

    Arguments:
        file (str): Path of file to md5sum

    Returns:
        str: md5sum
    """
    block_size = 8192
    m = hashlib.md5()
    with open(file, 'rb') as f:
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
    with requests.get(file_url, stream=True) as r:
        with open(out_file, 'wb') as f:
            shutil.copyfileobj(r.raw, f)


def decompress_gzip_file(gz_file, out_file):
    """Decompress a gzip file.

    Uncompressed given gzip file to out_file.

    Arguments:
        gz_file (str): Path of gzip file
        out_file (str): Path of target uncompressed out file
    """
    with gzip.open(gz_file, 'rb') as f_in:
        with open(out_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def get_md5_from_url(file_url):
    """Read md5sum file from url

    Assumes the md5 file is the file url with .md5 appended

    Arguments:
        file_url (str): URL of file

    Returns:
        str: md5sum
    """
    return requests.get(f'{file_url}.md5').text.split(' ')[0]


def check_md5(file_url, file_path):
    """Check md5 and remove file if incorect.

    Arguments:
        file_url (str): URL of file
        file_path (str): Path of file
    """
    expected_md5 = get_md5_from_url(file_url)
    downloaded_md5 = md5sum_file(file_path)
    if downloaded_md5 != expected_md5:
        os.unlink(file_path)
        sys.exit(f'[ERROR] MD5 check failed, removing {file_path}.')


def get_db(base_dir, db='progenomes'):
    """Download GUNC DB.

    Arguments:
        base_dir (str): Path of output directory
    """
    base_url = 'https://swifter.embl.de/~fullam/gunc/'
    if db == 'progenomes':
        file_name = 'gunc_db_progenomes2.1.dmnd.gz'
    elif db == 'gtdb':
        file_name = 'gunc_db_gtdb95.dmnd.gz'
    else:
        sys.exit(f'[ERROR] DB {db} unknown. Allowed: progenomes, gtdb')
    gz_file_url = f'{base_url}{file_name}'
    gz_file_path = os.path.join(base_dir, file_name)
    file_url = gz_file_url.replace('.gz', '')
    out_file = gz_file_path.replace('.gz', '')

    if not os.path.isdir(base_dir):
        sys.exit(f'[ERROR] Output Directory {base_dir} doesnt exist.')

    print('[INFO] DB downloading...')

    download_file(gz_file_url, gz_file_path)

    print('[INFO] DB download successful.')
    print('[INFO] Computing DB md5sum...')

    check_md5(gz_file_url, gz_file_path)

    print('[INFO] md5sum check successful.')
    print('[INFO] Uncompressing file...')

    decompress_gzip_file(gz_file_path, out_file)

    print('[INFO] Decompression complete.')
    print('[INFO] Computing DB md5sum...')

    check_md5(file_url, out_file)
    os.unlink(gz_file_path)

    print('[INFO] md5sum check successful.')
    print('[INFO] DB download successful.')
    print(f'[INFO] DB path: {out_file}')
