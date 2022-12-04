# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-15 13:56:44
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2022-12-04 09:21:02
 * @FilePath: /metaSC/PyLib/seqPipe/collect_gene_ko.py
 * @Description:

rule gene_clu:
    input:
        all_faa     = "{any}.faa",
    output:
        all_100     = "{any}-clu_100.tsv",
        all_clu_faa = "{any}-clu_rep.faa",
        all_clu     = "{any}-clu.tsv",
        clu_flat    = "{any}-clu_seq.faa.clu",
    threads: THREADS
    shell:
        '''
        source gene_clust.sh
        mmseq_ {input.all_faa} {threads}
        '''

rule genome_2_ko:
    input:
        all_100      = "{any}-clu_100.tsv",
        all_clu      = "{any}-clu.tsv",
        annot_kofam  = "{any}-clu_rep-kofam.tsv",
        annot_eggnog = "{any}-clu_rep-eggnog.tsv",
        annot_tcdb   = "{any}-clu_rep-tcdb.tsv",
        annot_mantis = "{any}-clu_rep-mantis.tsv",
    output:
        gene_annots  = "{any}-all_ko.csv",
    params:
        annot_prefix = "{any}-clu_rep-.*.tsv",
    shell:
        '''
            set +u
            source ~/.miniconda_init
            conda activate python39

        python -m workflow.remote.gene_annot \
            --all-100      {input.all_100} \
            --all-clu      {input.all_clu} \
            --annot-prefix {params.annot_prefix} \
            --gene-annots  {output.gene_annots} \
        '''
"""

import os
import re
from pathlib import Path
from typing import Union

import click
import pandas as pd
from PyLib.PyLibTool.file_info import verbose_import
from PyLib.reader.iters import read_table

logger = verbose_import(__name__, __doc__)


def drop_gene(
    gene: str, subsets  # type: tuple[Union[set, list, dict], bool]
) -> bool:
    subset, subset_is_contig = subsets
    if subset:
        if subset_is_contig:
            return gene.rsplit("_", 1)[0] not in subset
        else:
            return gene not in subset
    return False


class gene2KO:
    class gene_ko_iter:
        def __init__(self, filename: Path):
            self.filename = filename

        def __call__(self):
            raise NotImplementedError

    class ghost(gene_ko_iter):
        def __call__(self):
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[0:2]
                        yield gene, ko

    class kofam(gene_ko_iter):
        def __call__(self):
            with open(self.filename) as text:
                for values in read_table(text):
                    if len(values) > 1:
                        gene, ko = values[1:3]
                        yield gene, ko

    class eggnog(gene_ko_iter):
        def __call__(self):
            i_KEGG_ko = 11
            with open(self.filename) as text:
                for values in read_table(text):
                    kos = values[i_KEGG_ko]
                    if len(kos) > 1:
                        gene = values[0]
                        for ko in kos.split(","):
                            yield gene, ko[3:]

    class mantis(gene_ko_iter):
        KO_PATTERN = re.compile("(^|;)(K\\d{5})(;|$)")

        def __call__(self):
            with open(self.filename) as text:
                for values in read_table(text):
                    Ref_Hits = values[2]
                    kos = re.findall(self.KO_PATTERN, Ref_Hits)
                    if len(kos) >= 1:
                        gene = values[0]
                        for ko in kos:
                            yield gene, ko[1]

    ## collect gene KO
    annoters = [ghost, mantis, kofam, eggnog]

    def get_gene_KOs(self, annoters=None):
        """Only keep the first match:
        >>> gene_KOs.setdefault(gene, ko)"""
        if annoters is None:
            annoters = self.annoters
        gene_KOs = {}
        for annoter, file in zip(annoters, self.ann_files):
            if not os.path.isfile(file):
                continue

            gene_KOs_ = {}
            for gene, ko in annoter(file)():
                gene_KOs_.setdefault(gene, []).append(ko)
            for gene in gene_KOs_:
                if gene not in gene_KOs:
                    gene_KOs[gene] = ":".join(gene_KOs_[gene])

        return gene_KOs

    def get_gene_ko(self):
        return pd.Series(self.get_gene_KOs(), name="ko")

    def __init__(self, pattern):

        if pattern is None:
            ann_files = []
        else:
            if not isinstance(pattern, list):
                pattern = [pattern]
            ann_files_ = self._infer_ann_files(pattern)
            ann_files = [ann_files_.get(i, Path()) for i, _ in enumerate(self.annoters)]

            if not any(ann_files):
                raise FileNotFoundError(f"pattren(s) '{pattern}' donot match any file!")

        self.ann_files = ann_files

    def _infer_ann_files(self, patterns):
        ann_files = {i: Path() for i, _ in enumerate(self.annoters)}

        for pattern in reversed(patterns):
            pattern_re = re.compile(pattern.name)
            for file in pattern.parent.iterdir():
                if pattern_re.search(file.name):
                    for i, source in enumerate(self.annoters):
                        if source.__name__ in file.name.lower():
                            ann_files[i] = file
        return ann_files


def load_gene_abd(
    gene_count, method="rpb", melt=False, threshold_0=0.001
) -> pd.DataFrame:
    if isinstance(gene_count, pd.DataFrame):
        gene_count_raw = gene_count
    else:
        gene_count_raw = pd.read_csv(gene_count, index_col=0, sep="\t", comment="#")

    if method == "count":
        gene_abd = gene_count_raw.iloc[:, 5:]
    elif method == "rpb":
        gene_abd = gene_count_raw.pipe(lambda df: df.iloc[:, 5:].T / df.iloc[:, 4]).T
    elif method == "tpm":
        gene_abd = load_gene_abd(gene_count_raw, "rpb", threshold_0=threshold_0).pipe(
            lambda df: df * 1e6 / df.agg("sum")
        )
    else:
        raise ValueError("unknown method")

    gene_abd_cut: pd.DataFrame = gene_abd.pipe(lambda df: df * (df > threshold_0))

    if melt:
        return (
            gene_abd_cut.reset_index()
            .pipe(
                lambda df: df.melt(
                    id_vars=["Geneid"],
                    value_vars=df.columns,
                    var_name="bam",
                    value_name=method,
                )
            )
            .pipe(lambda df: df[df[method] > 0])
        )
    else:
        return gene_abd


def load_rep2all(all_100_: Path, all_clu_: Path):
    all_100 = pd.read_csv(all_100_, sep="\t", header=None, names=["rep100", "all"])
    all_clu = pd.read_csv(all_clu_, sep="\t", header=None, names=["rep", "rep100"])
    rep2all = all_100.merge(all_clu)[["rep", "all"]]
    return rep2all


def get_all_gene_ko(gene_ko: pd.Series, rep2all: pd.DataFrame):
    ko_exploded = (
        gene_ko.apply(lambda x: x.split(":") if x.startswith("K") else [])
        .explode()
        .dropna()
        .reset_index()
    )
    all_gene_ko = ko_exploded.merge(rep2all, left_on="index", right_on="rep")[
        ["all", "ko"]
    ].set_index("all")

    return all_gene_ko


def main(
    annot_prefix: Path,
    all_100_path: Path,
    all_clu_path: Path,
    all_gene_ko_path: Path,
):
    gene_ko = gene2KO(annot_prefix).get_gene_ko()
    rep2all = load_rep2all(all_100_path, all_clu_path)
    all_gene_ko = get_all_gene_ko(gene_ko, rep2all)
    all_gene_ko.to_csv(all_gene_ko_path)


@click.command()
@click.option("--loglevel", default="INFO", type=str, help="set level of logger")
@click.option("--all-100", type=Path, help="")
@click.option("--all-clu", type=Path, help="")
@click.option("--annot-prefix", type=Path, help="")
@click.option("--gene-annots", type=Path, help="")
def run(
    loglevel: str,
    all_100: Path,
    all_clu: Path,
    annot_prefix: Path,
    gene_annots: Path,
):
    logger.setLevel(level=loglevel.upper())  # info

    gene_annots.parent.mkdir(parents=True, exist_ok=True)

    logger.warning(">>> job start")
    state = main(annot_prefix, all_100, all_clu, gene_annots)
    logger.warning(">>> job finish")
    if state == 0:
        logger.info("success!")


if __name__ == "__main__":
    run()
