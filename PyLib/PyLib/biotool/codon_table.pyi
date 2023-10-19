from _typeshed import Incomplete as Incomplete
from collections.abc import Generator

class AminoAcidInfo:
    SF: Incomplete
    aa: Incomplete
    GC_rank: Incomplete
    N: Incomplete
    S: Incomplete
    C: Incomplete
    MW: Incomplete
    codons: Incomplete
    def __init__(self, aa, SF, GC_rank, N, S, C, MW, codons) -> None: ...

codon_aas_11: Incomplete
codon_dict: Incomplete

def calculate_SCU(gene, errorfile_handle) -> None: ...
def codon_grouper(n, iterable, fillvalue: Incomplete | None = ...): ...
def screen_out_ambiguous_codons(sequence) -> None: ...
def random_permutation(iterable) -> None: ...
def return_coding_sequence(fna_handle) -> Generator[Incomplete, None, None]: ...
def scan_for_stop_codons(DNA_sequence) -> None: ...
def ARSC_and_MW_from_amino_acids(protein_sequence) -> None: ...
def ARSC_MW_from_nucleotides(sequence) -> None: ...
def weighted_by_abundance(gene_95_results_handle, abundance_table, metadata_handle, gene_cluster_map_handle: Incomplete | None = ..., output_flag: Incomplete | None = ...): ...
