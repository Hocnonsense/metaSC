# -*- coding: utf-8 -*-
"""
 * @Date: 2021-09-23 17:29:53
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-10-11 12:33:30
 * @FilePath: /metaSC/PyLib/biotool/diamond.py
 * @Description: Run diamond with python
"""

import os
import subprocess


def diamond_check_db(fastadir, fastaid, suffix, outdir):
    """Make database with DIAMOND."""
    diamond_db = os.path.join(outdir, fastaid)
    if not os.access(diamond_db, os.F_OK):
        subprocess.check_call(
            (
                "diamond makedb "
                "--in " + os.path.join(fastadir, fastaid + suffix) + " "
                "-d " + diamond_db + " "
                "-p 1 --quiet"
            ),
            shell=True,
        )
    return diamond_db


def diamond_blastp(fastadir, query, db, suffix, level, outdir):
    """DIAMOND blastp alignment."""
    level = "" if level == "fast" else f"--{level}"

    # make database with DIAMOND
    diamond_db = diamond_check_db(fastadir, db, suffix, outdir)

    blast_out = os.path.join(outdir, f"{query}_{db}.blast")
    if not os.access(blast_out, os.F_OK):
        subprocess.check_call(
            (
                "diamond blastp "
                "-d " + diamond_db + " "  # database
                "-q " + os.path.join(fastadir, query + suffix) + " "  # query like "<dir>/faa/<query>.faa"
                "-o " + blast_out + " "
                "-k 10 -e 0.001 --id 30 --query-cover 70 --subject-cover 70 "
                "-b 10 --dbsize 1000000000 -p 1 "
                "--quiet " + level
            ),
            shell=True,
        )
    return blast_out
