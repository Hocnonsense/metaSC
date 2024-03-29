# -*- coding: utf-8 -*-
"""
 * @Date: 2020-10-02 22:32:18
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-10-04 10:45:46
 * @FilePath: /metaSC/PyLib/biotool/checkm/manual.py
 * @Description:
        Try to cluster each scaffold by the Mark Gene on it
        First, get a list of each Gene on each Scafflod
        Second, group these Sacffolds by Gene
        Third, change data format to a list sorted by distance

        As we can see, when calculating Contamination, It get count(gene) per gene_set per mark_gene_set
            So, I'll point these gene to scaffolds by different color, and we can see it more clearly
        However, there are only a linely color space
            So, I'll first group gene into gene sets, then draw trees based on these distance, Then change the tree to a list.
            Alternatively, I'll split those scaffolds into several different scaffolds for each gene set.

 * @TODU: function
"""

import argparse

import argparse
import os
from collections import OrderedDict
from sys import stderr

from .reload import (
    reload_binIdToBinMarkerSets,
    reload_resultsManagers,
)


def str_to_bool(s: str) -> bool:
    """
    @description: 识别 "是否" 类型的参数
    @param s:
        true, false, 1, or 0.
    @return {type}
    """
    s = s.lower()
    if s not in {"true", "false", "1", "0"}:
        raise ValueError("Expected a bool, not %s" % s)
    return s in {"true", "1"}


def define_parameters(parser: argparse.ArgumentParser) -> None:
    """@description: 参数设置
    Define command line parameters. This base method defines a --verbose
    flag. Overrides should call super.

    Examples include positional arguments
        `parser.add_argument('variant', nargs='?',
        help='Simulation variant.')`
    options
        `parser.add_argument('--seed', default='000000',
        help='Simulation seed.')`.
    and flags
        `parser.add_argument('--verbose', action='store_true',
        help='Enable verbose logging.')`.
    """
    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose logging."
    )


def define_option(
    parser: argparse.ArgumentParser,
    name: str,
    default,
    help: str,
    short_name: str = "",
) -> None:
    """@description: 可选参数设置
    Add an option with the given name and datatype to the parser.
    """
    datatype = type(default)
    if short_name:
        parser.add_argument(
            "-" + short_name,
            "--" + name,
            type=type(default),
            default=default,
            help="({}; {}) {}".format(datatype.__name__, default, help),
        )
    else:
        parser.add_argument(
            "--" + name,
            type=type(default),
            default=default,
            help="({}; {}) {}".format(datatype.__name__, default, help),
        )


def define_parameter_bool(
    parser: argparse.ArgumentParser, name: str, default, help: str
) -> None:
    """@description: "bool" 类型的参数设置. 有默认值
    Add a boolean option parameter to the parser. The CLI input can be
    `--name`, `--no_name`, `--name true`, `--name false`, `--name 1`,
    `--name 0`, `--name=true`, etc. The default can be True or False, and
    changing it won't affect any of those explicit input forms. This method
    adds the default value to the help text.
    """
    default = bool(default)
    examples = "true or 1" if default else "false or 0"
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--" + name,
        nargs="?",
        default=default,
        const="true",  # needed for nargs='?'
        type=str_to_bool,
        help="({}; {}) {}".format("bool", examples, help),
    )
    group.add_argument(
        "--no_" + name, dest=name, action="store_false", help="Like {}=0".format(name)
    )


def send_parser(description):
    """Send an argparse.ArgumentParser to user.
    @param description: description
    """
    return argparse.ArgumentParser(description=description)


def markerGeneToContig(markerSet, resultsManager, reportContig=True):
    """Calculate marker gene sets distribute in each contig
    @param markerSet         : _checkm.MarkerSet
    @param resultsManager    : _checkm.ResultsManager
    @param reportContig      : True to report contig, False to keep silence
    @return (contigs, ms_wight)
                             : contigs: dict -> {contig_name: set(index of marker set)}
                             : ms_weight: dict -> [multiCopy/length for each markerSet]
    """
    contigs = {}
    ms_wight = {}

    # ms: {'PF04760.10', }
    for i, ms in enumerate(markerSet.markerSet):
        # marker: 'PF04760.10'
        multiCopy = 0
        for marker in ms:
            for hit in resultsManager.markerHits.get(marker, []):
                contig_name = hit.target_name[: hit.target_name.rfind("_")]
                if contig_name not in contigs:
                    contigs[contig_name] = set()
                contigs[contig_name].add(i)
                multiCopy += 1
        ms_wight[i] = multiCopy - 1 if multiCopy else 0
        ms_wight[i] /= len(ms)

    if reportContig:
        report = {}
        for contig_name in contigs:
            count = len(contigs[contig_name])
            report[count] = report.get(count, 0) + 1
        print("UID", markerSet.UID, ":", report)

    return contigs, ms_wight


def generateAdditionalCsv(out_file, markerSet, contigs, ms_wight):
    """Generate additional csv in R.
    @useage: After get the file, in R:
         marksets <- paste0(mks_dir, bin, ".marker_gene.csv")  # ! output by checkMarkToR.py
         mm <- mmload(
             additional = read.csv(marksets)
         se = mmplot(
             color_by = 'cluster_id',
    """
    DUPLICATE = int(
        len(ms_wight) * 5 / 4
    )  # more than one marker genes on this contig. should be carefull

    sorted_markerSets = sorted(list(ms_wight), key=lambda i: ms_wight[i], reverse=True)
    out_markerSets = (
        sorted_markerSets[:] if len(sorted_markerSets) >= 10 else sorted_markerSets
    )
    print([(len(markerSet.markerSet[ms]), ms_wight[ms]) for ms in out_markerSets])

    contig_valuse = {}
    for contig_name, markerSets in contigs.items():
        # cluster_id represented for different markerSets
        cluster_id = DUPLICATE
        # non-duplicated contigs
        if len(markerSets) == 1:
            cluster_id = sorted_markerSets.index(markerSets.pop())

        # remarkable most duplicated markerSets.
        sorted_ms = [1 if i in markerSets else "" for i in out_markerSets]
        contig_valuse[contig_name] = [cluster_id, *sorted_ms]

    print("save to file: {}".format(out_file))
    with open(out_file, "w") as fout:
        # title
        print(
            "contig_id",
            "cluster_id",
            *list(range(len(out_markerSets))),
            sep=",",
            file=fout
        )
        for contig_name in contig_valuse:
            print(contig_name, *contig_valuse[contig_name], sep=",", file=fout)


def getBinsToModify(
    ckmap: OrderedDict, Compl_THRESHOLD: int = 70, Contai_THRESHOLD: int = 10
):
    """OK, now, find which are good and which should be modified."""
    goodbins = []
    mdf_bins = []
    uid_nums = []
    # Bin Id, Marker lineage (UID), Completeness, Contamination
    for binID, (UID, Compl, Contai) in ckmap.items():
        if Compl > Compl_THRESHOLD:
            if Contai < Contai_THRESHOLD:
                goodbins.append(binID)
            else:
                mdf_bins.append(binID)
                uid_nums.append(UID)
    else:
        print(
            """
            Now, you get:
                {} good bins
                {} bins to modify
                {} bins useless
            more information:
                {}
            """.format(
                len(goodbins),
                len(mdf_bins),
                len(ckmap) - len(goodbins) - len(mdf_bins),
                str(ckmap),
            ),
            file=stderr,
        )
    return goodbins, mdf_bins, uid_nums
    # print("Output list to txt", file=stderr)
    # with open("checkm-good.tmp", "w") as gout, open("checkm-modify.tmp", "w") as mout, open("checkm-to-r.tmp", "w") as uout:
    #    gout.write("\n".join(goodbins)+"\n")
    #    mout.write("\n".join(mdf_bins)+"\n")
    #    uout.write("\n".join([" ".join(match)
    #                          for match in zip(mdf_bins, uid_nums)])+"\n")


def main():
    parser = send_parser(__doc__)
    parser.add_argument(
        "checkM_result_dir",
        type=str,
        help="the dir of your checkM output. "
        "For example, if you run `checkm lineage_wf -x fa -t 40 ${dbins} ./02_checkm`, "
        "then your `checkM_result_dir` is `./02_checkm`",
    )
    parser.add_argument(
        "binId",
        nargs="?",
        type=str,
        default="",
        help="the binID you are going to modified. "
        "generally, Completeness > 70 and Contamination > 10. (view in checkM output)",
    )
    parser.add_argument(
        "UID",
        nargs="?",
        type=str,
        default="",
        help="the UID of the best match of given bin. " "given in row `Marker lineage`",
    )
    parser.add_argument(
        "out_dir",
        nargs="?",
        type=str,
        default="",
        help="output dir of the file. "
        "The output_file name is always: binID+`.marker_gene.csv`",
    )
    argv = parser.parse_args()
    checkM_result_dir = argv.checkM_result_dir
    binId = argv.binId
    UID = argv.UID
    out_dir = argv.out_dir

    binIdToBinMarkerSets = reload_binIdToBinMarkerSets(checkM_result_dir)
    resultsManagers = reload_resultsManagers(checkM_result_dir, binIdToBinMarkerSets)

    if UID:
        resultsManager = resultsManagers[binId]
        print("read given UID: {} in given bins: {}".format(UID, binId))
        for markerSet in binIdToBinMarkerSets[binId].markerSets:
            if markerSet.UID == UID:
                break
        contigs, ms_wight = markerGeneToContig(markerSet, resultsManager)

        if out_dir:
            out_file = os.path.join(out_dir, ".".join([binId, r"marker_gene.csv"]))
            generateAdditionalCsv(out_file, markerSet, contigs, ms_wight)
    elif binId:
        resultsManager = resultsManagers[binId]
        print("read all UID in given bins: {}".format(binId))
        for markerSet in binIdToBinMarkerSets[binId].markerSets:
            contigs, _ = markerGeneToContig(markerSet, resultsManager)
    else:
        print("read all UID in all bins")
        for binId, resultsManager in resultsManagers.items():
            print("binId:", binId)
            for markerSet in binIdToBinMarkerSets[binId].markerSets:
                contigs, _ = markerGeneToContig(markerSet, resultsManager)


if __name__ == "__main__":
    main()
# example:
# D:\Code\python\CheckM-master>tmp.py D:\Files\TODU\MetaWork\results bins.029_sub.contigs UID4444 D:\Files\TODU\MetaWork
