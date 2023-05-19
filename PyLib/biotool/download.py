# -*- coding: utf-8 -*-
"""
 * @Date: 2021-02-03 11:09:20
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-05-19 11:19:29
 * @FilePath: /metaSC/PyLib/biotool/download.py
 * @Description:
        download genome from net
"""

import ftplib
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Union
from urllib.request import urlopen, urlretrieve

from Bio import Entrez

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
            logger.debug('Downloading "{}" to "{}"'.format(url, download_file))

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


def check_entrez_email(email=None):
    """try to use git config user.email"""
    if email:
        Entrez.email = email
        return
    if Entrez.email is None:
        out, err = runsh_safe("git config user.email")
        Entrez.email = out


class RetriveUrl:
    @staticmethod
    def verify_format(sequence_id: str) -> bool:
        return False

    @classmethod
    def _retrieve_url(cls, sequence_id: str) -> str:
        raise NotImplementedError

    @classmethod
    def retrieve_url(cls, sequence_id: str):
        if cls.verify_format(sequence_id):
            return cls._retrieve_url(sequence_id)

    @classmethod
    def download(cls, sequence_id: str, filename: Path, overwrite=True):
        return download(cls._retrieve_url(sequence_id), filename, overwrite)


class RefSeqURL(RetriveUrl):
    @staticmethod
    def verify_format(sequence_id):
        return sequence_id.startswith("GCA") or sequence_id.startswith("GCF")

    refseq_base_ftp_url = "ftp.ncbi.nlm.nih.gov"
    refseq_genomes_url = "genomes/all"

    @classmethod
    def _retrieve_files(cls, sequence_id: str):
        gcx_id = sequence_id

        gcx, number = gcx_id.split(".")[0].split("_")
        gcx_url = "/".join(
            [gcx] + [number[i : i + 3] for i in range(0, len(number), 3)]
        )

        ftp = ftplib.FTP(cls.refseq_base_ftp_url)
        _ = ftp.login()
        _ = ftp.cwd(cls.refseq_genomes_url + "/" + gcx_url)
        folder = ftp.nlst()[0]
        _ = ftp.cwd(folder)
        files = ftp.nlst()
        _ = ftp.quit()

        gcx_download_fmt = "https://" + "/".join(
            [cls.refseq_base_ftp_url, cls.refseq_genomes_url, gcx_url, folder, "{}"]
        )
        return gcx_download_fmt, folder, files

    @classmethod
    def ls_url(cls, sequence_id: str, file_suffix: str = ""):
        gcx_download_fmt, folder, files = cls._retrieve_files(sequence_id)

        download_files = [ff for ff in files if file_suffix in ff]
        return [gcx_download_fmt.format(ff) for ff in download_files]

    @classmethod
    def _retrieve_url(cls, sequence_id: str) -> str:
        gcx_download_fmt, folder, files = cls._retrieve_files(sequence_id)

        for ff in files:
            if (folder + "_genomic.fna.gz") == ff:
                return gcx_download_fmt.format(ff)
        raise KeyError("no item found")

    @classmethod
    def download(cls, sequence_id: str, filename: Path, overwrite=True):
        download(cls._retrieve_url(sequence_id), f"{filename}.gz", overwrite=overwrite)
        runsh_safe(f"gunzip {filename}.gz")
        return filename


class GWHSeqURL(RetriveUrl):
    @staticmethod
    def verify_format(sequence_id):
        return sequence_id.startswith("GWH")

    @staticmethod
    def _retrieve_url(sequence_id: str) -> str:
        gwh_id = sequence_id
        gwh_base_ftp_url = "download.big.ac.cn/gwh"

        gwh_base_url = "https://bigd.big.ac.cn/gwh/api/public/assembly"
        gwh_genome_url = json.loads(urlopen(gwh_base_url + "/" + gwh_id).read())[
            "ftpPathDna"
        ]
        return "ftp://" + "/".join([gwh_base_ftp_url, gwh_genome_url])

    @classmethod
    def download(cls, sequence_id: str, filename: Path, overwrite=True):
        download(cls._retrieve_url(sequence_id), f"{filename}.gz", overwrite=overwrite)
        runsh_safe(f"gunzip {filename}.gz")
        return filename


class IMGSeqURL(RetriveUrl):
    __cookies: Path = None  # type: ignore

    @classmethod
    def img_login(
        cls, name: str = "", passwd: str = "", cookies="cookies", overwrite=False
    ):
        if overwrite or not os.path.exists(cookies):
            os.system(
                "curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode "
                f"'login={name}' --data-urlencode 'password={passwd}' -c {cookies}"
            )
        cls.__cookies = Path(cookies)

    @classmethod
    @property
    def cookies(cls):
        if not isinstance(cls.__cookies, Path):
            raise FileNotFoundError(
                f"you should generate a cookie via {type(cls)}.img_login"
            )
        return str(cls.__cookies)

    @staticmethod
    def verify_format(sequence_id):
        return sequence_id.startswith("GWH")

    @classmethod
    def _retrieve_url(cls, sequence_id: str) -> str:
        """
        * @description:
        * @param {str} img_id: IMGXXXXXX
        * @return {str} url of sequence of IMG (if sequence exists), however, should download with 'curl'
        """
        img_id = sequence_id
        img_base_url = "https://genome.jgi.doe.gov"

        img_search_url = "https://genome.jgi.doe.gov/portal/ext-api/downloads/get-directory?organism="
        os.system(f"curl '{img_search_url}{img_id}' -b {cls.__cookies} > {img_id}.xml")
        tree = ET.parse(img_id + ".xml")
        root = tree.getroot()
        IMG_Data = [i for i in root if i.attrib["name"] == "IMG Data"][0]
        for i in IMG_Data:
            if i.attrib["filename"].endswith(".fna"):
                return img_base_url + i.attrib["url"]
        raise KeyError("no item found")

    @classmethod
    def download(cls, sequence_id: str, filename: Path, overwrite=True):
        if (not os.path.isfile(filename)) or overwrite:
            img_url = cls._retrieve_url(sequence_id)
            os.system(f"curl '{img_url}' -b {cls.__cookies} > {filename}")
        return filename


class GenBankSeqURL(RetriveUrl):
    """
    reference:
        - https://support.nlm.nih.gov/knowledgebase/article/KA-03436/en-us
    """

    pattrens = {
        k: re.compile(v)
        for k, v in {
            "1": r"^([A-Z]\d{5})(.\d+)?",
            "2": r"^([A-Z]{2}\d{6})(.\d+)?",
        }.items()
    }

    @classmethod
    def verify_format(cls, sequence_id):
        matches = {k: re.match(v, sequence_id) for k, v in cls.pattrens.items()}
        return any(i is not None for i in matches.values())

    @staticmethod
    def _retrieve_url(sequence_id: str) -> str:
        genebank_id = sequence_id
        genebank_base_url = (
            "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id="
        )
        return genebank_base_url + genebank_id

    @classmethod
    def download(cls, sequence_id: str, filename: Path, overwrite=True):
        download(cls._retrieve_url(sequence_id), f"{filename}", overwrite=overwrite)
        return filename


class WGSSeqURL(RefSeqURL):
    """
    reference:
        - https://www.ncbi.nlm.nih.gov/genbank/wgs/
    """

    pattrens = {
        k: re.compile(v)
        for k, v in {
            "1": r"^[A-Z]{4}\d{8}(.\d*)?",
            "1901": r"^[A-Z]{6}\d{9}(.\d*)?",
        }.items()
    }
    RETMAX = 100

    @classmethod
    def verify_format(cls, sequence_id):
        matches = {k: re.match(v, sequence_id) for k, v in cls.pattrens.items()}
        return any(i is not None for i in matches.values())

    @classmethod
    def _retrieve_url(cls, sequence_id: str) -> str:
        check_entrez_email()

        wgs_id = sequence_id
        with Entrez.esearch(db="genome", term=wgs_id, retmax=cls.RETMAX) as handle:
            search = Entrez.read(handle)

        if len(search["IdList"]) > 1:
            raise KeyError("more than one record returned")
        with Entrez.esummary(db="genome", id=search["IdList"][0]) as handle:
            summary = Entrez.read(handle)
        assert isinstance(summary, Entrez.Parser.ListElement)
        if len(summary) == 1:
            refseq_id = summary[0]["Assembly_Accession"]
        else:
            raise KeyError("more than one record returned")
        return super()._retrieve_url(refseq_id)


def _download_fna(sequence_id: str, output: Path, cookies):
    """
    @param sequence_id: name of reference genome, possibly starts with [GCA, GCF, GWH, IMG]
    @param output: Path to store output
                   - if it ends with ".fna", then output.name will be used as filename
                   - otherwise, file will be stored at output / f"{sequence_id}.fna"
    @param cookies: for downloading IMG genomes
    @return: path to downloaded genomes
    """
    if output.name.endswith(".fna"):
        fna_file = output
        output = output.parent
    else:
        fna_file = output / f"{sequence_id}.fna"

    output.mkdir(parents=True, exist_ok=True)
    logger.info(f"download to {fna_file}")
    for rurl in (
        RefSeqURL,
        GWHSeqURL,
        IMGSeqURL,
        WGSSeqURL,
        GenBankSeqURL,
    ):
        if rurl.verify_format(sequence_id):
            logger.info(f"detect format {rurl.__name__}")
            rurl.download(sequence_id, fna_file)
            break
    else:
        raise NotImplementedError
    return fna_file


def download_fna(
    sequence_id: str,
    output: Union[str, Path] = "./",
    cookies="cookies",
    overwrite=False,
    retry=0,
):
    basicConfig()
    fna_file = Path(output) / f"{sequence_id}.fna"

    if overwrite or not fna_file.is_file():
        for i in range(retry, -1, -1):
            try:
                fna_file = _download_fna(sequence_id, fna_file, cookies)
            except ConnectionRefusedError as e:
                logger.warning(f"{sequence_id} failed at {e.strerror}, retry ({i})")
                if fna_file.is_file():
                    runsh_safe(f"/bin/rm {fna_file}")
                if i == 0:
                    logger.warning(f"{sequence_id} failed after {retry} tries")
                    raise e
                time.sleep(5)
                continue
            else:
                logger.info(f"{fna_file} downloaded successfully.")
                break
        else:
            raise Exception(f"failed to download {fna_file}.")
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
