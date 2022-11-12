# -*- coding: utf-8 -*-
"""
 * @Date: 2022-11-12 20:31:16
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-12 21:19:41
 * @FilePath: /metaSC/PyLib/biotool/metabolic_overlap.py
 * @Description:
    extract from https://github.com/ericHester/metabolicOverlap/tree/master/scripts2

    for fasta in *.fna; do
        fna_to_rxn <$fasta | rxn_expandinfo >${fasta%.fna}.rxn
    done
    rxn_to_connections *.rxn | tee connections | lists_to_matrix -c3 >matrix

"""

import re
from collections import defaultdict
from pathlib import Path
from typing import Generator

ModelSEEDDatabase = Path("/home/hwrn/Data/Database2/ModelSEEDDatabase")


def fna_to_rxn():
    """
    give a genome, return gene and function annotation (pathway)
    """


rxn_item = tuple[str, tuple[list[str], str, list[str]]]


def rxn_expandinfo(reaction_file: Path) -> Generator[rxn_item, None, None]:
    """
    id, reaction
    "rxn48575", [["cpd37299[0]"], ">", ["cpd00794[0]", "cpd02725[0]"]]
    "rxn48575\tcpd37299[0]\t>\tcpd00794[0]+cpd02725[0]"
    """
    idx_id, idx_eq, idx_dir = (0, 5, 7)

    def cleanequation(raw, direction):
        def cleancompound(raw):
            if raw == " ":
                return ""
            return re.match(".*(cpd[0-9]+\[[^\]]*\])", raw).group(1)

        ins, outs = [
            [cleancompound(y) for y in x.split("+")] for x in re.split("<=>|=>|<=", raw)
        ]
        if direction == "<":  # Make sure it's always input > output
            direction = ">"
            ins, outs = outs, ins
        return (sorted(ins), direction, sorted(outs))

    reactiondb = {}
    filename = ModelSEEDDatabase / "Biochemistry" / "reactions.tsv"
    with open(filename) as file:
        reactiondb = {
            row[idx_id]: row[idx_id + 1 :]
            for row in (line.rstrip("\n").split("\t") for line in file)
        }

    with open(reaction_file) as file:
        reactions = (line.rstrip("\n") for line in file)
        return (
            (id, cleanequation(reactiondb[id][idx_eq], reactiondb[id][idx_dir]))
            for id in reactions
        )


def rxn_to_connections(
    genome_rxns: dict[str, Generator[rxn_item, None, None]]
) -> Generator[tuple[str, str, str], None, None]:
    organisms: dict[str, dict[str, set]] = {}
    for genome, genome_rxn_i in genome_rxns.items():
        organisms[genome] = {"ins": set(), "outs": set()}
        for rxn in genome_rxn_i:
            organisms[genome]["ins"].update(rxn[1][0])
            organisms[genome]["outs"].update(rxn[1][2])
            if rxn[1][1] == "=":
                organisms[genome]["ins"].update(rxn[1][2])
                organisms[genome]["outs"].update(rxn[1][0])

    for name, current in organisms.items():
        for cpd in current["outs"]:
            for othername, other in organisms.items():
                # if other == current: continue
                if cpd in other["ins"]:
                    yield (name, cpd, othername)


def lists_to_matrix(rxn_connections: Generator[tuple[str, str, str], None, None]):
    matrix: defaultdict = defaultdict(lambda: defaultdict(float))
    for name, cpd, othername in rxn_connections:
        matrix[name][othername] += 1

    # I made this. I have no clue how it works
    params = set([param for (_, params) in matrix.items() for param in params.keys()])
    return [(["name"] + list(params))] + [
        ([name] + [str(fndict[param]) for param in params])
        for (name, fndict) in matrix.items()
    ]
