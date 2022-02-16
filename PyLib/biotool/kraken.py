# -*- coding: utf-8 -*-
"""
 * @Date: 2022-02-11 15:44:32
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-15 16:48:23
 * @FilePath: /metaSC/PyLib/biotool/kraken.py
 * @Description:
"""
from typing import Final, Generator, List, Literal, TextIO, Tuple, Union

import numpy as np
import pandas as pd
from tqdm import trange

main_lvls: Final = ["R", "K", "D", "P", "C", "O", "F", "G", "S"]


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
        if len(report_vals) < 5:
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
    in_file: TextIO, level: Literal["k", "p", "c", "o", "f", "g", "s"] = None
):
    """{
        'in_header': ['percents', 'all_reads', 'exact_reads', 'level_type', 'name'],
        'out_header': [format_taxon, reads_at_level]
    }"""

    if level is None:
        for taxon, all_reads in kraken_iter(in_file, include_x=False):
            yield format_taxon(taxon), all_reads
    else:
        levels: Final = ["k", "p", "c", "o", "f", "g", "s"]
        leveli: Final = levels.index(level) + 1  # type: ignore
        supresslevels: Final = [f"{i}__" for i in levels[leveli:]]

        last_reports: List[Tuple[str, int]] = [("", 0)]
        for taxon, all_reads in kraken_iter(in_file, include_x=False):
            if sum(supresslevel in taxon for supresslevel in supresslevels):
                continue
            while last_reports and not taxon.startswith(last_reports[-1][0]):
                last_taxon, last_reads = last_reports.pop()
                yield format_taxon(last_taxon, leveli), last_reads
            last_taxon, last_reads = last_reports[-1]
            last_reports[-1] = last_taxon, last_reads - all_reads
            last_reports.append((taxon, all_reads))
        total_reads = last_reports.pop(0)[1]
        while last_reports:
            last_taxon, last_reads = last_reports.pop()
            yield format_taxon(last_taxon, leveli), last_reads
        print("total reads:", -total_reads)


def rarefaction(otus, step, repeat=10, seed=0):
    prng = np.random.RandomState(seed)  # reproducible results

    otus = np.array(otus)
    noccur = np.sum(otus)  # number of occurrences for each sample
    nvar = otus.shape[0]
    step = int(step)

    return pd.DataFrame(
        {
            step_: [
                sum(
                    np.bincount(
                        prng.choice(nvar, step_, p=otus / float(noccur)), minlength=nvar
                    )
                    > 0
                )
                for _ in trange(repeat, desc=f"repeating at {step_}")
            ]
            for step_ in trange(step, otus.sum(), step)
        }
    )


def kraken_rarefaction(filepath: str, step=2e6):
    kraken_report: pd.DataFrame = pd.read_csv(
        filepath,
        sep="\t",
        names=["percents", "all_reads", "exact_reads", "level_type", "uid", "name"],
        # usecols=["exact_reads"],
    ).dropna()
    otus = kraken_report[
        kraken_report["level_type"].apply(lambda x: x[0] not in ("U", "R"))
    ]["exact_reads"]
    rare = rarefaction(otus, step)
    rare_long = rare.melt(var_name="sample_size", value_name="otus").append(
        {"sample_size": otus.sum(), "otus": otus.size}, ignore_index=True
    )
    rare_long["sample"] = filepath
    return rare_long
