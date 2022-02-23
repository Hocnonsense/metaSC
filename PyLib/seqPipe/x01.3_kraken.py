# -*- coding: utf-8 -*-
"""
 * @Date: 2021-07-01 20:30:00
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-23 17:01:26
 * @FilePath: /metaSC/PyLib/seqPipe/x01.3_kraken.py
 * @Description:
"""

import os
import click
import logging
import glob
import pandas as pd
from typing import List, Literal
from tqdm import tqdm
from joblib import Parallel, delayed

from PyLib.biotool.kraken import kraken_rarefaction, kraken_level_filter


def report_kraken_level(
    level: Literal["k", "p", "c", "o", "f", "g", "s"], kraken_reports, threads=1
):
    def get_kraken_report_level(kraken_report: str):
        with open(kraken_report) as fi:
            return kraken_report, {k: v for k, v in kraken_level_filter(fi, level)}
            # if k.startswith("k__Bacteria") or k.startswith("k__Archaea")

    kraken_reports_level = dict(
        Parallel(threads, verbose=11)(
            delayed(get_kraken_report_level)(kraken_report)
            for kraken_report in kraken_reports
        )
    )

    return pd.DataFrame(kraken_reports_level, dtype=int).fillna(0)


@click.group()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option("-o", "--output-prefix", default="./", help="output prifix for file")
@click.argument("kraken_report_pattern", type=str, nargs=1)
@click.pass_context
def easy_click(ctx, loglevel: str, output_prefix: str, kraken_report_pattern: str):
    """deal with kraken report to summarize the metagenome sample"""
    logging.basicConfig(level=loglevel.upper())  # info
    kraken_reports = glob.glob(kraken_report_pattern)
    # output_prefix
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    ctx.obj["output_prefix"] = output_prefix
    # kraken_reports
    if not kraken_reports:
        raise FileNotFoundError(
            f"pattern {kraken_report_pattern} do not match any file"
        )
    else:
        ctx.obj["kraken_reports"] = kraken_reports


@easy_click.command()
@click.pass_context
@click.argument(
    "taxon", type=click.Choice(["k", "p", "c", "o", "f", "g", "s"]), default="g"
)
@click.option(
    "-t",
    "--threads",
    default=1,
    help="Threads to speed up. ",
)
def bar(ctx, taxon, threads):
    """Filter kraken to given TAXON; output a table of these genomes"""
    kraken_reports = ctx.obj["kraken_reports"]
    output = ctx.obj["output_prefix"] + f"kraken_{taxon}.csv"
    report_kraken_level(taxon, kraken_reports, threads).to_csv(output)


@easy_click.command()
@click.pass_context
@click.option(
    "-s",
    "--step",
    default=2e6,
    help="Step size for sample sizes in rarefaction curves.",
)
@click.option(
    "-r",
    "--repeat",
    default=10,
    help="Repeat time at each step. ",
)
@click.option(
    "-t",
    "--threads",
    default=1,
    help="Threads to speed up. ",
)
def rare(ctx, step, repeat, threads):
    """generate rarefaction report from KRAKEN_REPORTs"""
    kraken_reports: List[str] = ctx.obj["kraken_reports"]
    output = ctx.obj["output_prefix"] + f"kraken_rare.csv"
    kraken_rare = pd.concat(
        [
            Parallel(threads, verbose=11)(
                delayed(kraken_rarefaction)(report, int(step), repeat)
                for report in kraken_reports
            )
        ],
        axis=0,
    )
    kraken_rare.to_csv(output)


if __name__ == "__main__":
    easy_click(obj={})
