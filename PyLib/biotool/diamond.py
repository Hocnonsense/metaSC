# -*- coding: utf-8 -*-
"""
 * @Date: 2021-09-23 17:29:53
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-10 12:21:35
 * @FilePath: /metaSC/PyLib/biotool/diamond.py
 * @Description: Run diamond with python
"""

import os
from pathlib import Path
from PyLib.PyLibTool.file_info import basicConfig, verbose_import
import tempfile
import subprocess

logger = verbose_import(__name__, __doc__)


def remove_suffix(name: str, suffix=""):
    return name[: -len(suffix)] if name.endswith(suffix) else name


def diamond_check_db(fastapath: Path, outdir: Path, threads=1):
    """Make database with DIAMOND."""
    diamond_db = outdir / remove_suffix(fastapath.name, ".faa")
    cmd = (
        f"diamond makedb "
        f"    --in {fastapath} "
        f"    -d {diamond_db} "
        f"    -p {threads} --quiet"
    )
    logger.info(cmd)
    if not os.access(diamond_db, os.F_OK):
        subprocess.check_call(cmd, shell=True)
    return diamond_db


def diamond_blastp(
    query: Path, db: Path, level="sensitive", outdir: Path = Path("."), threads=1
):
    """DIAMOND blastp alignment."""
    level = "" if level == "fast" else f"--{level}"

    # make database with DIAMOND
    diamond_db = diamond_check_db(db, outdir, threads)

    blast_out = outdir / (
        remove_suffix(query.name, ".faa")
        + "_"
        + remove_suffix(db.name, ".faa")
        + ".blast"
    )
    cmd = (
        f"diamond blastp "
        f"    -d {diamond_db} "
        f"    -q {query} "
        f"    -o {blast_out} "
        f"    -k 10 -e 0.001 --id 30 --query-cover 70 --subject-cover 70 "
        f"    -b 10 --dbsize 1000000000 -p {threads} "
        f"    --quiet {level}"
    )
    logger.info(cmd)
    if not os.access(blast_out, os.F_OK):
        subprocess.check_call(cmd, shell=True)
    return blast_out


def diamond_1to1(
    genome1: Path, genome2: Path, out: Path, sensitive_level="sensitive", threads=1
):
    with tempfile.TemporaryDirectory() as tempdir_:
        tempdir = Path(tempdir_)
        basicConfig()
        logger.debug(f"diamand: ({tempdir}) : {genome1}, {genome2}")
        for i in range(1):
            try:
                out1 = diamond_blastp(
                    genome1, genome2, sensitive_level, tempdir, threads
                )
                out2 = diamond_blastp(
                    genome2, genome1, sensitive_level, tempdir, threads
                )
                logger.debug(
                    f"filter diamand finished: {tempdir} : {genome1}, {genome2}"
                )
                break
            except subprocess.CalledProcessError as e:
                logger.warning(f"{tempdir} : {genome1}, {genome2} failed for {i}")
                print(e.returncode, e.cmd, e.output, e.stderr, sep="\n\n")
                e1 = e
        else:
            raise Exception(e1)
        with (out.open("w") as fout, out1.open() as f1, out2.open() as f2):
            for fi in [f1, f2]:
                # fi.rollover()
                # fi.seek(0)
                while True:
                    inblock = fi.read(1024 * 1024)
                    if not inblock:
                        break
                    fout.write(inblock)
    logger.info(f"genemap wrote to {out}")
    return out
