# -*- coding: utf-8 -*-
"""
 * @Date: 2022-08-19 22:31:47
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-08-28 14:48:23
 * @FilePath: /metaSC/PyLib/seqPipe/gene_predict_cut.py
 * @Description:
"""

import os
from pathlib import Path
from typing import TextIO

import click

from PyLib.PyLibTool.file_info import verbose_import
from PyLib.seqPipe import x03_gene_cut as gene_cut


logger = verbose_import(__name__, __doc__)


def prodigal(
    genome: Path,
    out_dirname: Path = None,
    out_basename="",
    filter=False,
    meta=False,
    min_faa_len=33,
):
    if out_dirname is None:
        out_dirname = genome.parent
    out_prefix = out_dirname / (out_basename or genome.name)
    if filter:
        out_prefix_filter, out_prefix = out_prefix, out_prefix / "tmp"
    if out_prefix.exists():
        os.system(f"rm -r {out_prefix}")
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    os.system(
        f"prodigal -i {genome}"
        f"         {'-p meta' if meta else ''}"
        f"         -d {out_prefix}.fna"
        f"         -a {out_prefix}.faa"
        f"         -f gff "
        f"         -o {out_prefix}.gff"
        f"         -q",
    )
    if filter:
        out_prefix = prodigal_filter(
            out_prefix, out_prefix_filter, min_faa_len=min_faa_len
        )
        logger.info(f"finish writing to {out_prefix}")
        os.system(f"rm -r {out_prefix}")
    return out_prefix


def prodigal_filter(out_prefix: Path, filter_prefix: Path = None, min_faa_len=33):
    in_file_prefix = gene_cut.PrefixFiles(out_prefix)
    if filter_prefix is None:
        filter_prefix = Path(f"{in_file_prefix}_cut")
    with (
        open(f"{filter_prefix}.faa", "w") as faao,
        open(f"{filter_prefix}.fna", "w") as fnao,
        open(f"{filter_prefix}.gff", "w") as gffo,
    ):
        gene_cut.main(in_file_prefix, min_faa_len, (faao, fnao, gffo))
        faao.flush()
        fnao.flush()
        gffo.flush()

    return filter_prefix


def prodigal_faa2gff_iter(pfain: TextIO):
    for line in pfain:
        if line.startswith(">"):
            head, start, end, strand, other_ = line.strip().split(" # ")
            contigId = head[1:].rsplit("_", 1)[0]
            other = "ID=" + head[1:] + ";" + other_.split(";", 1)[1]
            yield (
                contigId,
                "Prodigal_v2.6.3-modify",
                "CDS",
                start,
                end,
                ".",
                strand[:-1] or "+",
                0,
                other,
            )


@click.command()
@click.option(
    "--gene-prefix",
    type=Path,
    help="path to gene prefix",
)
def prodigal_faa2gff_cmd(
    gene_prefix: Path,
):
    faa, _, gff = [Path(f"{gene_prefix}.{suffix}") for suffix in ("faa", "fna", "gff")]
    assert faa.is_file()
    with (
        open(faa, "w") as faai,
        open(gff, "w") as gffo,
    ):
        for item in prodigal_faa2gff_iter(faai):
            print(*item, sep="\t", file=gffo)
        gffo.flush()
    #
    return 0


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option(
    "--genome",
    type=Path,
    help="path to genome",
)
@click.option(
    "--gene-prefix",
    type=str,
    default="",
    help="path to gene output prefix",
)
@click.option(
    "--filter",
    is_flag=True,
    help="filter genome",
)
@click.option(
    "--meta",
    is_flag=True,
    help="meta mode",
)
@click.option(
    "--min-faa-len",
    default=33,
    help="min gene length",
)
def run(
    loglevel: str,
    genome: Path,
    gene_prefix: str,
    filter=False,
    meta=False,
    min_faa_len=33,
):
    logger.setLevel(level=loglevel.upper())  # info

    if not gene_prefix:
        out_dirname = None
        out_basename = ""
    else:
        gene_prefix_ = Path(gene_prefix)
        out_dirname, out_basename = gene_prefix_.parent, gene_prefix_.name

    logger.warning(">>> job start")
    prodigal(genome, out_dirname, out_basename, filter, meta, min_faa_len)
    logger.warning(">>> job finish")
    logger.info("success!")


if __name__ == "__main__":
    run()
