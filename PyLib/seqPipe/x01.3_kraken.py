# -*- coding: utf-8 -*-
"""
 * @Date: 2021-07-01 20:30:00
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-15 17:06:18
 * @FilePath: /metaSC/PyLib/seqPipe/x01.3_kraken.py
 * @Description:
"""

import os
import re
import click
import logging
import glob
from ast import literal_eval as eval
import sys
import pandas as pd
from typing import Callable, Dict, List, Union, TextIO, Literal

from PyLib.biotool.kraken import kraken_rarefaction, kraken_level_filter


def report_kraken_level(
    level: Literal["k", "p", "c", "o", "f", "g", "s"], kraken_reports
):
    kraken_reports_level = {}
    for kraken_report in kraken_reports:
        with open(kraken_report) as fi:
            kraken_reports_level[kraken_report] = {
                k: v
                for k, v in kraken_level_filter(fi, level)
                if k.startswith("k__Bacteria") or k.startswith("k__Archaea")
            }

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
def bar(ctx, taxon):
    """Filter kraken to given TAXON; output a table of these genomes"""
    kraken_reports = ctx.obj["kraken_reports"]
    output = ctx.obj["output_prefix"] + f"kraken_{taxon}.csv"
    report_kraken_level(taxon, kraken_reports).to_csv(output)


@easy_click.command()
@click.pass_context
@click.option(
    "-s",
    "--step",
    type=int,
    default=2e6,
    help="Step size for sample sizes in rarefaction curves.",
)
def rare(ctx, step):
    """generate rarefaction report from KRAKEN_REPORTs"""
    kraken_reports: List[str] = ctx.obj["kraken_reports"]
    output = ctx.obj["output_prefix"] + f"kraken_rare.csv"
    kraken_rare = pd.concat(
        [kraken_rarefaction(report, step) for report in kraken_reports], axis=0
    )
    kraken_rare.to_csv(output)


if __name__ == "__main__":
    easy_click(obj={})
