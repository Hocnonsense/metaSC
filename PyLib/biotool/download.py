# -*- coding: utf-8 -*-
"""
 * @Date: 2021-02-03 11:09:20
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-02-03 11:42:53
 * @FilePath: /HScripts/Python/mylib/biotool/download.py
 * @Description:
        download genome from net
"""

import sys
import os
from urllib.request import urlretrieve, urlopen
import time
import ftplib
import json
import xml.etree.ElementTree as ET
from mylib.tool.path import makedirs


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return (byte / 1048576)


class ReportHook():
    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """

        if blocknum == 0:
            self.start_time = time.time()

            if total_size > 0:
                print("Downloading file of size: {:.2f} MB\n".format(byte_to_megabyte(total_size)))
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = total_downloaded * 100.0 / total_size
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += ("{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec "
                           .format(percent_downloaded, byte_to_megabyte(download_rate),
                                   estimated_minutes, estimated_seconds))

            status += "        \r"
            print(status)


def download(url, download_file, overwrite=False, verbose=False):
    """
    Download a file from a url
    """

    if (not os.path.isfile(download_file)) or overwrite:
        try:
            if verbose:
                print('Downloading "{}" to "{}"\n'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
            print('\n')
        except EnvironmentError as e:
            sys.stderr.write('unable to download "{}"'.format(url))
            if verbose:
                print(e)
            exit()
        except Exception as e:
            sys.stderr.write('unable to download "{}"'.format(url))
            print("Fault!")
            print(e)

    elif verbose:
        print('File "{}" present\n'.format(download_file))


def retrieve_refseq_url(gcx_id):
    refseq_base_ftp_url = 'ftp.ncbi.nlm.nih.gov'
    refseq_genomes_url = 'genomes/all'

    gcx, number = gcx_id.split('.')[0].split('_')
    gcx_url = '/'.join([gcx] + [number[i:i + 3] for i in range(0, len(number), 3)])

    ftp = ftplib.FTP(refseq_base_ftp_url)
    _ = ftp.login()
    _ = ftp.cwd(refseq_genomes_url + '/' + gcx_url)
    folder = ftp.nlst()[0]
    _ = ftp.cwd(folder)
    files = ftp.nlst()
    _ = ftp.quit()

    for ff in files:
        if (folder + '_genomic.fna.gz') == ff:
            return 'https://' + '/'.join([refseq_base_ftp_url, refseq_genomes_url, gcx_url, folder, ff])


def retrive_gwh_url(gwh_id):
    gwh_base_ftp_url = 'download.big.ac.cn/gwh'

    gwh_base_url = 'https://bigd.big.ac.cn/gwh/api/public/assembly'
    gwh_genome_url = json.loads(urlopen(gwh_base_url + "/" + gwh_id).read())['ftpPathDna']
    return "ftp://" + "/".join([gwh_base_ftp_url, gwh_genome_url])


def img_login(name: str, passwd: str, cookies='cookies'):
    os.system(
        "curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode "
        f"'login={name}' --data-urlencode 'password={passwd}' -c {cookies}"
    )


def retrive_img_url(img_id, cookies='cookies'):
    """
     * @description:
     * @param {str} img_id: IMGXXXXXX
     * @return {str} url of sequence of IMG (if sequence exists), however, should download with 'curl'
    """
    if not os.path.exists(cookies):
        img_login('', '', cookies)
    img_base_url = 'https://genome.jgi.doe.gov'

    img_search_url = 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism='
    os.system(
        f"curl '{img_search_url}{img_id}' -b {cookies} > {img_id}.xml"
    )
    tree = ET.parse(img_id + '.xml')
    root = tree.getroot()
    IMG_Data = [i for i in root if i.attrib['name'] == 'IMG Data'][0]
    for i in IMG_Data:
        if i.attrib['filename'].endswith('.fna'):
            return img_base_url + i.attrib['url']


def download_fna(sequence_id: str, output, verbose=True, cookies='cookies'):
    makedirs(output)
    print(time.time())

    if sequence_id.startswith('GCA') or sequence_id.startswith('GCF'):
        download(retrieve_refseq_url(sequence_id), os.path.join(output, sequence_id + ".fna.gz"))
    elif sequence_id.startswith('GWH'):
        download(retrive_gwh_url(sequence_id), os.path.join(output, sequence_id + ".fna.gz"), verbose=True)
    elif sequence_id.startswith('IMG'):
        img_url = retrive_img_url(sequence_id, cookies)
        os.system(f"curl '{img_url}' -b {cookies} > {os.path.join(output, sequence_id + '.fna')}")
