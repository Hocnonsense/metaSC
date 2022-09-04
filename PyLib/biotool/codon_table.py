# -*- coding: utf-8 -*-
"""
 * @Date: 2017-06-02
 * @Editors: Jessica Bryant
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-04 12:26:41
 * @FilePath: /metaSC/PyLib/biotool/codon_table.py
 * @Description:
    This module holds dictionaries that contain codon and amino acid tables, and functions that
     caculate GC, N-ARSC, C-ARSC and Nc.

    Original citations for calculated metrics:

       Baudouin-Cornu P, Surdin-Kerjan Y, Marliere P, Thomas D. 2001. Molecular evolution of protein atomic
        composition. Science 293 297-300.

       Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87 23-29.

    Updated: 5/29/2017

    https://raw.githubusercontent.com/JessAwBryant/gene-characteristics/master/codon_table.py
"""

__author__ = "Jessica_Bryant"
__email__ = "jessawbryant@gmail.com"


import itertools as itertools
from operator import add as add
import random as random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import copy as cp


class AminoAcidInfo:
    """
    lists the synonymous family type (SF*), encoded amino acid (aa), ranking of nitrogen use of each codon relative
    to other the synonymous codons encoding the same aa (GC_rank), number of nitrogen atoms on each encoded aa side
    chain(N),  number of sulfur atoms on each encoded aa side chain (S), molecular weight of the encoded aa (MW).
    refs:
    Molecular Weight, N and S counts come from page 30 of 'Understanding Bioinformatics' by Zvelbil and Baum
    """

    def __init__(self, aa, SF, GC_rank, N, S, C, MW, codons):
        self.SF = SF
        self.aa = aa
        self.GC_rank = GC_rank
        self.N = N
        self.S = S
        self.C = C
        self.MW = MW
        self.codons = codons


codon_aas_11 = {
    "D": AminoAcidInfo("D", "SF2", 1, 0, 0, 2, 133.1032, {"GAT", "GAC"}),
    "E": AminoAcidInfo("E", "SF2", 1, 0, 0, 3, 147.1299, {"GAA", "GAG"}),
    "S": AminoAcidInfo(
        "S", "SF6", 1, 0, 0, 1, 105.0930, {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}
    ),
    "T": AminoAcidInfo("T", "SF4", 1, 0, 0, 2, 119.1197, {"ACT", "ACC", "ACA", "ACG"}),
    "Y": AminoAcidInfo("Y", "SF2", 1, 0, 0, 7, 181.1894, {"TAT", "TAC"}),
    "A": AminoAcidInfo("A", "SF4", 1, 0, 0, 1, 89.0935, {"GCT", "GCC", "GCA", "GCG"}),
    "V": AminoAcidInfo("V", "SF4", 0, 0, 0, 3, 117.1469, {"GTT", "GTC", "GTA", "GTG"}),
    "L": AminoAcidInfo(
        "L", "SF6", 1, 0, 0, 4, 131.1736, {"TTA", "TTG", "CTT", "CTC", "CTA", "CTG"}
    ),
    "I": AminoAcidInfo("I", "SF3", 1, 0, 0, 4, 131.1736, {"ATT", "ATC", "ATA"}),
    "P": AminoAcidInfo("P", "SF4", 1, 0, 0, 3, 115.1310, {"CCT", "CCC", "CCA", "CCG"}),
    "F": AminoAcidInfo("F", "SF2", 1, 0, 0, 7, 165.1900, {"TTT", "TTC"}),
    "G": AminoAcidInfo("G", "SF4", 1, 0, 0, 0, 75.0669, {"GGT", "GGC", "GGA", "GGG"}),
    "C": AminoAcidInfo("C", "SF2", 1, 0, 1, 1, 121.1590, {"TGT", "TGC"}),
    "M": AminoAcidInfo("M", "SF1", None, 0, 1, 3, 149.2124, {"ATG"}),
    "K": AminoAcidInfo("K", "SF2", 1, 1, 0, 4, 146.1882, {"AAA", "AAG"}),
    "W": AminoAcidInfo("W", "SF1", None, 1, 0, 9, 204.2262, {"TGG"}),
    "N": AminoAcidInfo("N", "SF2", 1, 1, 0, 2, 132.1184, {"AAT", "AAC"}),
    "Q": AminoAcidInfo("Q", "SF2", 1, 1, 0, 3, 146.1451, {"CAA", "CAG"}),
    "H": AminoAcidInfo("H", "SF2", 1, 2, 0, 4, 155.1552, {"CAT", "CAC"}),
    "R": AminoAcidInfo(
        "R", "SF6", 0.5, 3, 0, 4, 174.2017, {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}
    ),
}


codon_dict = {codon: aa for aa in codon_aas_11.values() for codon in aa.codons}


def calculate_SCU(gene, errorfile_handle):
    """
    This script takes a nucleotide sequence as input and returns measure of synonymous codon usage (SCU) Nc and
    average GC_rank (codons ranked by # nitrogen atoms relative to other synonymous codons) of the sequence.

    code_11 is a dictionary containing general information about each codon.
    """
    # This dictionary tallies codon and aa counts for the gene
    codon_usage_dictionary = {
        "codons": {codon: {"count": 0, "frequency": 0} for codon in codon_dict},
        "amino_acids": {
            aa: {"count": 0, "frequency": 0, "Ne": 0} for aa in codon_aas_11
        },
        "SF_type": {
            SF: [] for SF in ("absent", "SF1", "SF2", "SF3", "SF4", "SF6")
        },  # complete sub-dictionary indexed by 'SF types'
        "gene_codon_length": 0,
        "number_of_positions_with_GC_variability": 0,
        "sum_GC_rank": 0,
    }

    # gene parsing begins here.

    # remove start and stop codon from current DNA sequence
    current_gene_sequence = gene.seq[3 : (len(gene.seq) - 3)]

    # Parse through codons
    for i in xrange(0, len(current_gene_sequence), 3):
        # count each codon
        current_codon = current_gene_sequence[i : (i + 3)]

        # identify possible codons that will cause problems
        if current_codon in ["TAA", "TAG", "TGA"]:
            error1 = (
                gene.id.strip()
                + "\t"
                + "internal stop codon"
                + "\t"
                + current_codon
                + "\n"
            )
            errorfile_handle.writelines(error1)
            return ["Nan", "Nan"]

        if "N" in current_codon:
            error2 = gene.id.strip() + "\t" + "N present" + "\t" + current_codon + "\n"
            errorfile_handle.writelines(error2)
            return ["Nan", "Nan"]

        if current_codon not in codon_usage_dictionary["codons"]:
            print(current_codon, "hu", codon_usage_dictionary["codons"][current_codon])
            error3 = (
                gene.id.strip() + "\t" + "other issue" + "\t" + current_codon + "\n"
            )
            errorfile_handle.writelines(error3)
            return ["Nan", "Nan"]

        codon_usage_dictionary["codons"][current_codon]["count"] += 1
        codon_usage_dictionary["amino_acids"][codon_dict[current_codon].aa][
            "count"
        ] += 1
        codon_usage_dictionary["gene_codon_length"] += 1

        # count the number of codons with GC variability and their rank
        if codon_dict[current_codon].GC_rank != None:
            codon_usage_dictionary["number_of_positions_with_GC_variability"] += 1
            codon_usage_dictionary["sum_GC_rank"] += codon_dict[current_codon].GC_rank

    # calculate codon frequencies and update codon_usage_dictionary
    for k in codon_usage_dictionary["codons"]:
        amino_acid_total_usage = float(
            codon_usage_dictionary["amino_acids"][codon_dict[k].aa].count
        )

        if amino_acid_total_usage > 1:
            codon_usage_dictionary["codons"][k]["frequency"] = (
                float(codon_usage_dictionary["codons"][k]["count"])
                / amino_acid_total_usage
            )

        # if amino acid is rarely used (aka one or less times in protein)
        if amino_acid_total_usage < 2:
            codon_usage_dictionary["codons"][k]["frequency"] = "Rare"
            # codon_usage_dictionary['amino_acids'][codon_dict[k]['aa']]['Ne'] == 'absent'

    # Calculate Ne for each aa and record in codon_usage_dictionary['amino_acids']['aa']['Ne']
    for aa in codon_aas_11:
        if codon_usage_dictionary["amino_acids"][aa]["count"] < 2:
            codon_usage_dictionary["amino_acids"][aa]["Ne"] = "absent"
            codon_usage_dictionary["SF_type"]["absent"].append(aa)

        elif codon_usage_dictionary["amino_acids"][aa]["count"] > 1:
            squared_fequencies = []

            for each_codon in codon_aas_11[aa].codons:
                p = float(codon_usage_dictionary["codons"][each_codon]["frequency"])
                squared_fequencies.append(p * p)

            n = codon_usage_dictionary["amino_acids"][aa]["count"]
            F = ((n * (sum(squared_fequencies))) - 1) * 1.0 / (n - 1.0)

            # amino acid is rarely used if the numerator or denominator of equation for F == 0
            #
            if (((n * (sum(squared_fequencies))) - 1) * 1.0 / (n - 1.0)) == 0:
                codon_usage_dictionary["amino_acids"][aa]["Ne"] = "absent"
                codon_usage_dictionary["SF_type"]["absent"].append(aa)

            if round(F, 10) != 0:
                Ne = 1 / F
                codon_usage_dictionary["amino_acids"][aa]["Ne"] = Ne
                codon_usage_dictionary["SF_type"][codon_aas_11[aa].SF].append(aa)

        # print aa, n, F, Ne, squared_fequencies, codon_usage_dictionary['SF_type']

    # calculate GC rank sums
    GC_variability_summary = (
        float(codon_usage_dictionary["sum_GC_rank"])
        / float(codon_usage_dictionary["number_of_positions_with_GC_variability"])
        if codon_usage_dictionary["number_of_positions_with_GC_variability"] != 0
        else "Nan"
    )

    # parse through SF values to get final Nc
    Nc = 2.0
    # print gene.id

    for sf_type in ["SF2", "SF4", "SF6"]:
        aas_for_sf_type = codon_usage_dictionary["SF_type"][sf_type]
        sf_float = float(len(aas_for_sf_type))

        # if one of the groups other than SF3 isn't represented, then dont return an Nc value
        if sf_float == 0:
            return ["Nan", str(GC_variability_summary)]

        # This is actually F list
        F_list = [
            1 / codon_usage_dictionary["amino_acids"][x]["Ne"] for x in aas_for_sf_type
        ]
        av_F_for_sf_type = sum(F_list) / sf_float

        if sf_type == "SF2":
            Nc += 9 / av_F_for_sf_type
        if sf_type == "SF4":
            Nc += 5 / av_F_for_sf_type
        if sf_type == "SF6":
            Nc += 3 / av_F_for_sf_type

        codon_usage_dictionary["".join(sf_type + "_Ne")] = av_F_for_sf_type

    # SF3 is a special case because only I is SF3
    if len(codon_usage_dictionary["SF_type"]["SF3"]) == 1:
        Nc += codon_usage_dictionary["amino_acids"]["I"]["Ne"]

    if len(codon_usage_dictionary["SF_type"]["SF3"]) == 0:
        Nc += 1 / (
            (codon_usage_dictionary["SF2_Ne"] + codon_usage_dictionary["SF4_Ne"])
            / float(2)
        )

    if Nc > 61:
        Nc = 61
    return [
        str(round(Nc, 4)),
        str(GC_variability_summary),
    ]  # , codon_usage_dictionary, code_11


def codon_grouper(n, iterable, fillvalue=None):
    """grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"""
    args = [iter(iterable)] * n
    return itertools.izip_longest(fillvalue=fillvalue, *args)


def screen_out_ambiguous_codons(sequence):
    # removes codons with Ns
    return reduce(
        add,
        [
            reduce(add, x)
            for x in codon_grouper(3, sequence)
            if len(set(x) - set(["A", "T", "C", "G"])) < 1
        ],
    )


def random_permutation(iterable):
    """Random selection from itertools.permutations(iterable, r)"""
    pool = tuple(iterable)
    r = len(iterable)
    # check for a stop codon. if there is one, pick a new sequence.
    random_sequence = reduce(add, random.sample(pool, r))
    for codon in codon_grouper(3, random_sequence):
        codon = reduce(add, codon)
        if codon in ["TAA", "TAG", "TGA"]:
            return random_permutation(iterable)
    return random_sequence


def return_coding_sequence(fna_handle):
    """creates a generator that cycles through CDS features in gbff files"""
    x = 0
    g = SeqIO.parse(fna_handle, "genbank")
    while True:
        genome = g.next()
        for gene in genome.features:
            if gene.type != "CDS":
                continue
            if "pseudo" in gene.qualifiers:
                continue
            if "exception" in gene.qualifiers:
                continue
            is_entire_gene = gene.location.__repr__()
            if "BeforePosition" in is_entire_gene:
                continue
            if "AfterPosition" in is_entire_gene:
                continue
            sequence = gene.extract(genome.seq)
            if gene.qualifiers["codon_start"][0] != "1":
                print(gene.qualifiers["codon_start"][0])
                first_codon = int(gene.qualifiers["codon_start"][0]) - 1
                sequence = sequence[first_codon:]
                print(gene.qualifiers["locus_tag"][0])
            gene_seq = SeqRecord(
                sequence,
                id=gene.qualifiers["locus_tag"][0]
                + " "
                + gene.qualifiers["product"][0],
            )
            yield gene_seq


def scan_for_stop_codons(DNA_sequence):
    for codon in codon_grouper(3, DNA_sequence[0:-3]):
        codon = reduce(add, codon)
        if reduce(add, codon) in ["TAA", "TAG", "TGA"]:
            return "True"
    return False


def ARSC_and_MW_from_amino_acids(protein_sequence):
    """
    This functions takes an amino acid sequence coded in single letters and returns N/C ARSC and Molecular Weight
    N and S counts come from page 30 of 'Understanding Bioinformatics' by Zvelbil and Baum
    molecular weights from http://www.webqc.org/aminoacids.php
    """
    # remove whitespaces and '*' (termination) that prodigal adds to end of aa sequences
    protein_sequence = protein_sequence.strip().strip("*")

    # remove ambigous aas from the string
    protein_sequence_no_Xs = "".join(
        [x for x in protein_sequence if x != "X" if x != "-"]
    )

    # caculate ARSC and Molecular Weight
    N_ARSC = sum(map(lambda x: codon_aas_11[x].N, protein_sequence_no_Xs)) / float(
        len(protein_sequence_no_Xs)
    )
    C_ARSC = sum(map(lambda x: codon_aas_11[x].C, protein_sequence_no_Xs)) / float(
        len(protein_sequence_no_Xs)
    )
    av_molecular_weight = sum(
        map(lambda x: codon_aas_11[x].MW, protein_sequence_no_Xs)
    ) / float(len(protein_sequence_no_Xs))

    return [
        str(round(N_ARSC, 4)),
        str(round(av_molecular_weight, 4)),
        str(round(C_ARSC, 4)),
    ]


def ARSC_MW_from_nucleotides(sequence):
    """
    This functions takes a protein sequence coded in nucleotides and returns N-ARSC.
    Internal stop codons will cause this script problems
    """

    # remove whitespaces and stop codon at end of sequences
    sequence = sequence.strip()[0:-3]

    # caculate NARSC
    aa_sequence_length = len(sequence) / 3.0
    total_nitrogen_atoms = sum(
        map(lambda x: codon_dict[reduce(add, x)].N, codon_grouper(3, sequence))
    )
    av_ARSC = total_nitrogen_atoms / aa_sequence_length

    total_carbon_atoms = sum(
        map(lambda x: codon_dict[reduce(add, x)].C, codon_grouper(3, sequence))
    )
    av_C_ARSC = total_carbon_atoms / aa_sequence_length

    total_molecular_weight = sum(
        map(lambda x: codon_dict[reduce(add, x)].MW, codon_grouper(3, sequence))
    )
    av_molecular_weight = total_molecular_weight / aa_sequence_length

    return (
        str(round(av_ARSC, 4)),
        str(round(av_molecular_weight, 4)),
        str(round(av_C_ARSC, 4)),
        str(round(av_ARSC / av_C_ARSC, 4)),
    )


def weighted_by_abundance(
    gene_95_results_handle,
    abundance_table,
    metadata_handle,
    gene_cluster_map_handle=None,
    output_flag=None,
):
    """
    abundance weighted averages by samples
    """

    gene_results_table = open(gene_95_results_handle)
    abundance_table = pd.read_csv(abundance_table, sep="\t", index_col=0)
    metadata = pd.read_csv(metadata_handle, sep=",", index_col=0)

    # translate cluster handles
    if gene_cluster_map_handle is not None:
        mOTU_seq_to_clusters = pd.read_csv(
            gene_cluster_map_handle, index_col=0, sep="\t", names=["clusters"]
        )

    # itterate through lines (individual gene results)
    # assuming first entry is gene.id and second entry is the sample
    first_line = gene_results_table.readline().strip().split("\t")
    sample_names = metadata.index.values
    values_to_average = first_line[2 : len(first_line)]

    print(values_to_average, metadata.index)
    model = {x: {"sum": 0, "count": 0} for x in values_to_average}
    tally_dictionary = dict([(sample, cp.deepcopy(model)) for sample in sample_names])

    line = gene_results_table.readline()
    print(first_line)
    while line:
        data = line.strip().split("\t")
        line_values = dict(zip(first_line, data))

        if gene_cluster_map_handle is not None:
            cluster_sequence = mOTU_seq_to_clusters.loc[
                line_values["gene.id"], "clusters"
            ]
        else:
            cluster_sequence = line_values["gene.id"]

        # print( data[0], line_values['gene.id'], cluster_sequence)

        if cluster_sequence in abundance_table.index:
            for sample in sample_names:
                abundance = abundance_table.loc[cluster_sequence, sample]
                for key in values_to_average:
                    if line_values[key] != "Nan" and line_values[key] != "None":
                        tally_dictionary[sample][key]["sum"] += (
                            float(line_values[key]) * abundance
                        )
                        tally_dictionary[sample][key]["count"] += abundance
        else:
            print("cluster not in annotation file", cluster_sequence)
        line = gene_results_table.readline()

    # calculate averages
    final_dictionary = dict(
        [(sample, {x: "NA" for x in values_to_average}) for sample in sample_names]
    )

    for sample in tally_dictionary:
        for value in values_to_average:
            if tally_dictionary[sample][value]["count"] > 0:
                final_dictionary[sample][value] = (
                    tally_dictionary[sample][value]["sum"]
                    / tally_dictionary[sample][value]["count"]
                )
            final_dictionary[sample]["N_ARSC_count"] = tally_dictionary[sample][value][
                "sum"
            ]

    final_table = pd.DataFrame(final_dictionary)
    # flag the values with '_av_all_catalog' so I'll know where they came later:
    names = dict([(x, x + "_weigted_95per_clusters") for x in final_table.index])
    final_table = final_table.rename(names)

    # Concatinate results with metadata
    GC_meta_table = pd.concat([metadata, final_table.T], axis=1)

    output_handle = gene_95_results_handle + "_weighted_abundance"
    if output_flag is not None:
        output_handle = output_handle + output_flag
    GC_meta_table.to_csv(output_handle, sep="\t")

    return ()
