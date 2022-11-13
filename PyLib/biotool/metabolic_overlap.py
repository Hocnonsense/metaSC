# -*- coding: utf-8 -*-
"""
 * @Date: 2022-11-12 20:31:16
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-11-13 10:35:53
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
from typing import Generator, NamedTuple, Literal

ModelSEEDDatabase = Path("/home/hwrn/Data/Database2/ModelSEEDDatabase")


def load_SEED_reactions(ModelSEEDDatabase=ModelSEEDDatabase):
    idx_id = 0
    filename = ModelSEEDDatabase / "Biochemistry" / "reactions.tsv"
    with open(filename) as file:
        reactiondb = {
            row[idx_id]: row[idx_id + 1 :]
            for row in (line.rstrip("\n").split("\t") for line in file)
        }
    return reactiondb


def fna_to_rxn():
    """
    give a genome, return gene and function annotation (pathway)
    """
    raise NotImplementedError


class Reaction(NamedTuple):
    id: str
    ins: list[str]
    outs: list[str]
    direction: Literal["=", ">", "<"]


def cleancompound(raw: str):
    """
    >>> cleancompound("cpd37299[0]")
    'cpd37299'
    """
    COMPOUND_RE = re.compile(".*(cpd[0-9]+\[[^\]]*\])")
    if (match := re.match(COMPOUND_RE, raw)) is not None:
        return match.group(1)
    return ""  # raw = " "


def cleanequation(raw: str, direction: Literal["=", ">", "<"]) -> Reaction:
    ins, outs = [
        [cleancompound(y) for y in x.split("+")] for x in re.split("<=>|=>|<=", raw)
    ]
    if direction == "<":  # Make sure it's always input > output
        direction = ">"
        ins, outs = outs, ins
    return Reaction("", sorted(ins), sorted(outs), direction)


def rxn_expandinfo(reaction_file: Path) -> Generator[Reaction, None, None]:
    """
    id, reaction
    "rxn48575", [["cpd37299[0]"], ">", ["cpd00794[0]", "cpd02725[0]"]]
    "rxn48575\tcpd37299[0]\t>\tcpd00794[0]+cpd02725[0]"
    """
    idx_eq, idx_dir = (5, 7)

    reactiondb = load_SEED_reactions()

    with open(reaction_file) as file:
        reactions = (line.rstrip("\n") for line in file)
        return (
            cleanequation(reactiondb[id][idx_eq], reactiondb[id][idx_dir])._replace(
                id=id
            )
            for id in reactions
        )


def rxn_to_connections(
    genome_rxns: dict[str, Generator[Reaction, None, None]]
) -> Generator[tuple[str, str, str], None, None]:
    organisms: dict[str, dict[str, set]] = {}
    for genome, genome_rxn_i in genome_rxns.items():
        organisms[genome] = {"ins": set(), "outs": set()}
        for rxn in genome_rxn_i:
            organisms[genome]["ins"].update(rxn.ins)
            organisms[genome]["outs"].update(rxn.outs)
            if rxn.direction == "=":
                organisms[genome]["ins"].update(rxn.outs)
                organisms[genome]["outs"].update(rxn.ins)

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
