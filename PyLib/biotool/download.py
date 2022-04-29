# -*- coding: utf-8 -*-
"""
 * @Date: 2021-02-03 11:09:20
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 11:16:17
 * @FilePath: /metaSC/PyLib/biotool/download.py
 * @Description:
        download genome from net
"""

import ftplib
import json
import os
import sys
import time
from typing import Union
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.request import urlopen, urlretrieve

from PyLib.PyLibTool.file_info import basicConfig, verbose_import
from PyLib.tool.shell import runsh_safe

logger = verbose_import(__name__, __doc__)


def byte_to_megabyte(byte):
    """
    Convert byte value to megabyte
    """

    return byte / 1048576


class ReportHook:
    def __init__(self):
        self.start_time = time.time()

    def report(self, blocknum, block_size, total_size):
        """
        Print download progress message
        """
        basicConfig()

        if blocknum == 0:
            self.start_time = time.time()

            if total_size > 0:
                logger.info(
                    "Downloading file of size: {:.2f} MB".format(
                        byte_to_megabyte(total_size)
                    )
                )
        else:
            total_downloaded = blocknum * block_size
            status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

            if total_size > 0:
                percent_downloaded = min(total_downloaded * 100.0 / total_size, 100.0)
                # use carriage return plus sys.stderr to overwrite stderr
                download_rate = total_downloaded / (time.time() - self.start_time)
                estimated_time = (total_size - total_downloaded) / download_rate
                estimated_minutes = int(estimated_time / 60.0)
                estimated_seconds = estimated_time - estimated_minutes * 60.0
                status += "{:3.2f} %  {:5.2f} MB/sec {:2.0f} min {:2.0f} sec ".format(
                    percent_downloaded,
                    byte_to_megabyte(download_rate),
                    estimated_minutes,
                    estimated_seconds,
                )

            status += "        \r"
            logger.debug(status)


def download(url, download_file=None, overwrite=False):
    """
    Download a file from a url
    """
    basicConfig()
    if download_file is None:
        download_file = url.split("/")[-1]

    if (not os.path.isfile(download_file)) or overwrite:
        try:
            logger.info('Downloading "{}" to "{}"'.format(url, download_file))

            urlretrieve(url, download_file, reporthook=ReportHook().report)
            return download_file
        except EnvironmentError as e:
            sys.stderr.write('unable to download "{}"'.format(url))
            logger.error(e)
            exit()
        except Exception as e:
            sys.stderr.write('unable to download "{}"'.format(url))
            logger.error("Fault!")
            logger.error(e)
            return ""

    else:
        logger.warning('File "{}" present'.format(download_file))
        return download_file


def retrieve_refseq_url(gcx_id):
    refseq_base_ftp_url = "ftp.ncbi.nlm.nih.gov"
    refseq_genomes_url = "genomes/all"

    gcx, number = gcx_id.split(".")[0].split("_")
    gcx_url = "/".join([gcx] + [number[i : i + 3] for i in range(0, len(number), 3)])

    ftp = ftplib.FTP(refseq_base_ftp_url)
    _ = ftp.login()
    _ = ftp.cwd(refseq_genomes_url + "/" + gcx_url)
    folder = ftp.nlst()[0]
    _ = ftp.cwd(folder)
    files = ftp.nlst()
    _ = ftp.quit()

    for ff in files:
        if (folder + "_genomic.fna.gz") == ff:
            return "https://" + "/".join(
                [refseq_base_ftp_url, refseq_genomes_url, gcx_url, folder, ff]
            )


def retrive_gwh_url(gwh_id):
    gwh_base_ftp_url = "download.big.ac.cn/gwh"

    gwh_base_url = "https://bigd.big.ac.cn/gwh/api/public/assembly"
    gwh_genome_url = json.loads(urlopen(gwh_base_url + "/" + gwh_id).read())[
        "ftpPathDna"
    ]
    return "ftp://" + "/".join([gwh_base_ftp_url, gwh_genome_url])


def img_login(name: str, passwd: str, cookies="cookies"):
    os.system(
        "curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode "
        f"'login={name}' --data-urlencode 'password={passwd}' -c {cookies}"
    )


def retrive_img_url(img_id, cookies="cookies"):
    """
    * @description:
    * @param {str} img_id: IMGXXXXXX
    * @return {str} url of sequence of IMG (if sequence exists), however, should download with 'curl'
    """
    if not os.path.exists(cookies):
        img_login("", "", cookies)
    img_base_url = "https://genome.jgi.doe.gov"

    img_search_url = (
        "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism="
    )
    os.system(f"curl '{img_search_url}{img_id}' -b {cookies} > {img_id}.xml")
    tree = ET.parse(img_id + ".xml")
    root = tree.getroot()
    IMG_Data = [i for i in root if i.attrib["name"] == "IMG Data"][0]
    for i in IMG_Data:
        if i.attrib["filename"].endswith(".fna"):
            return img_base_url + i.attrib["url"]


def download_fna(
    sequence_id: str,
    output: Union[str, Path] = "./",
    cookies="cookies",
    overwrite=False,
):
    basicConfig()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    fna_file = output / f"{sequence_id}.fna"
    logger.info(f"download to {fna_file}")

    if overwrite or not os.path.isfile(fna_file):
        if sequence_id.startswith("GCA") or sequence_id.startswith("GCF"):
            download(retrieve_refseq_url(sequence_id), f"{fna_file}.gz")
            runsh_safe(f"gunzip {fna_file}.gz")
        elif sequence_id.startswith("GWH"):
            download(retrive_gwh_url(sequence_id), f"{fna_file}.gz")
            runsh_safe(f"gunzip {fna_file}.gz")
        elif sequence_id.startswith("IMG"):
            img_url = retrive_img_url(sequence_id, cookies)
            os.system(f"curl '{img_url}' -b {cookies} > {fna_file}")
    else:
        logger.warning(f"{fna_file} already exists, skip.")

    return fna_file


def retrieve_refseq_url_ls(gcx_id: str, file_suffix: str = ""):
    refseq_base_ftp_url = "ftp.ncbi.nlm.nih.gov"
    refseq_genomes_url = "genomes/all"

    gcx, number = gcx_id.split(".")[0].split("_")
    gcx_url = "/".join([gcx] + [number[i : i + 3] for i in range(0, len(number), 3)])

    ftp = ftplib.FTP(refseq_base_ftp_url)
    _ = ftp.login()
    _ = ftp.cwd(refseq_genomes_url + "/" + gcx_url)
    folder = ftp.nlst()[0]
    _ = ftp.cwd(folder)
    files = ftp.nlst()
    _ = ftp.quit()

    gcx_download_fmt = "https://" + "/".join(
        [refseq_base_ftp_url, refseq_genomes_url, gcx_url, folder, "{}"]
    )
    download_files = [ff for ff in files if file_suffix in ff]
    return [gcx_download_fmt.format(ff) for ff in download_files]
