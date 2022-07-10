# -*- coding: utf-8 -*-
"""
 * @Date: 2022-02-11 15:44:32
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-07-10 12:49:25
 * @FilePath: /metaSC/PyLib/biotool/kraken.py
 * @Description:
"""

from pathlib import Path
from typing import Dict, List, Final, Generator, Literal, TextIO, Tuple

import pandas as pd

from PyLib.PyLibTool.file_info import verbose_import


logger = verbose_import(__name__, __doc__)


main_lvls: Final = ["R", "K", "D", "P", "C", "O", "F", "G", "S"]
levels: Final[List[Literal["k", "p", "c", "o", "f", "g", "s"]]] = [
    "k",
    "p",
    "c",
    "o",
    "f",
    "g",
    "s",
]


def process_kraken_report_line(curr_str: str):
    split_str = curr_str.strip().split("\t")
    try:
        int(split_str[1])
    except ValueError:
        # If header line, skip
        return
    percents = float(split_str[0])
    all_reads = int(split_str[1])
    exact_reads = int(split_str[2])
    level_type = split_str[3]
    name = split_str[-1].strip()
    # Determine level based on number of spaces
    level_num = int((len(split_str[-1]) - len(name)) / 2)
    return name, level_num, level_type, all_reads, exact_reads, percents


def kraken_iter(
    in_file: TextIO,
    include_x=False,
) -> Generator[Tuple[str, int], None, None]:
    def check_report(level):
        return include_x or level != "x"

    curr_path = []
    prev_lvl_num = -1
    for line in in_file:
        report_vals = process_kraken_report_line(line)
        # If header line, skip
        if report_vals is None or len(report_vals) < 5:
            continue
        # Get relevant information from the line
        name, level_num, level_type, all_reads, exact_reads, percents = report_vals
        if level_type == "U":
            continue
        # Create level name
        if level_type not in main_lvls:
            level_type = "x"
        elif level_type == "D":  # elif level_type == "K"
            level_type = "K"
        level_str = level_type.lower() + "__" + name
        # Determine full string to add
        if prev_lvl_num == -1:
            # First level
            prev_lvl_num = level_num
            curr_path.append(level_str)
        else:
            # Move back if needed
            while level_num != (prev_lvl_num + 1):
                prev_lvl_num -= 1
                curr_path.pop()
            # Update
            curr_path.append(level_str)
            prev_lvl_num = level_num
            # Print if at non-traditional level and that is requested
            if check_report(level_type):
                # Print all ancestors of current level followed by ;
                taxon = ";".join(
                    (
                        string
                        for string in curr_path
                        if check_report(string[0]) and string[0] != "r"
                    )
                )
                # Print final level and then number of reads
                yield taxon, all_reads


def format_taxon(taxon: str, formatleveli=7):
    formatlevels: Final = "k__", "p__", "c__", "o__", "f__", "g__", "s__"
    taxon_dict: Final = {i[:3]: i for i in taxon.split(";")}
    return ";".join(
        (taxon_dict.get(i, i) for i, _ in zip(formatlevels, range(formatleveli)))
    )


def kraken_level_filter(
    in_file: TextIO,
    level: Literal["k", "p", "c", "o", "f", "g", "s"],
    prokaryotes_only=True,
):
    leveli: Final = levels.index(level)
    last_taxon = "r__root"  # impossible in any taxon
    for taxon, all_reads in kraken_iter(in_file, include_x=False):
        if prokaryotes_only and (
            not (taxon.startswith("k__Archaea") or taxon.startswith("k__Bacteria"))
        ):
            continue
        taxoni = taxon.split(";")[-1][0]
        if taxoni in levels[:leveli]:
            continue
        if taxoni == level:
            yield taxon, all_reads
            last_taxon = taxon
        else:
            if taxon.startswith(last_taxon):
                continue
            yield taxon, all_reads
            last_taxon = taxon


def kraken_summary(outdir: Path, *files):
    sample_all: Dict[str, Tuple[int, int]] = {}
    all_reads = pd.DataFrame()
    for file in files:
        logger.info(f"reading {file}")
        with open(file) as fi:
            uroot = fi.readline().strip().split("\t")
            root = fi.readline().strip().split("\t")
            total_reads = int(uroot[1]) + int(root[1])
            sample_all[file] = total_reads, int(uroot[1])
        for level in levels:
            with open(file) as fi:
                taxon_reads_prok = pd.DataFrame(
                    kraken_level_filter(fi, level, prokaryotes_only=True),  # type: ignore
                    columns=["taxon", "reads"],
                )
            taxon_reads_prok["taxon"] = taxon_reads_prok["taxon"].apply(
                lambda x: "d" + x[1:]
            )
            taxon_reads_prok["level"] = level
            taxon_reads_prok["sample"] = file

            all_reads = pd.concat([all_reads, taxon_reads_prok])

    outdir.mkdir(exist_ok=True)
    pd.DataFrame(sample_all, index=["total_reads", "unroot_reads"]).T.to_csv(
        outdir / "kraken-root.csv"
    )
    for level, taxon_reads in all_reads.groupby("level"):
        taxon_reads.pivot_table(
            values="reads", index="taxon", columns="sample", aggfunc=sum, fill_value=0
        ).to_csv(outdir / f"kraken-{level}.csv")
