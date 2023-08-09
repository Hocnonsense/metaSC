# -*- coding: utf-8 -*-
"""
 * @Date: 2023-08-08 20:18:58
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-08-09 15:04:57
 * @FilePath: /metaSC/PyLib/biotool/kmer_assembly.py
 * @Description:
"""

from typing import Generator
from Bio import SeqIO, Seq


class KmerContig:
    __slots__ = ["kmer", "starts", "ends", "circle", "single"]

    def __init__(self, kmer) -> None:
        self.kmer = str(kmer)
        # contig starts with this kmer, and ends with another kmer
        self.starts: set[str] = set()
        # contig ends with this kmer, and starts with another kmer
        self.ends: set[str] = set()
        # circula contig both ends with this kmer
        self.circle: set[str] = set()
        # contigs with this kmer, while the other end do not overlap with any other seq
        self.single: set[str] = set()

    def __repr__(self) -> str:
        report = self.report_dict()
        report.pop("kmer")
        return str(report)

    def is_end(self):
        if not self.circle and (not self.starts or not self.ends):
            return [*self.starts | self.ends]

    def report_dict(self):
        return {
            "kmer": self.kmer,
            "starts": ",".join(sorted(self.starts)),
            "ends": ",".join(sorted(self.ends)),
            "circle": ",".join(sorted(self.circle)),
            "single": ",".join(sorted(self.single)),
        }


def get_kmer_str(seq: Seq.Seq):
    seq_, seq_rc = seq.upper(), seq.reverse_complement().upper()
    if seq_ < seq_rc:
        return str(seq_), ""
    else:
        return str(seq_rc), "rc"


class KmerAssembly:
    def __init__(self, file: str, k: int) -> None:
        self.file = file
        self.k = k
        self.kmers: dict[str, KmerContig] = {}
        # contigs with both ends unconnected
        self.isolate: set[str] = set()

    def parse(self) -> Generator[SeqIO.SeqRecord, None, None]:
        yield from SeqIO.parse(self.file, "fasta")

    def get_kmer(self):
        for seq in self.parse():
            seq_s, seq_e = seq[: self.k].seq, seq[-self.k :].seq
            if seq_s == seq_e:
                self[get_kmer_str(seq_s)[0]].circle.add(seq.id)
            else:
                seq_s_, is_rc = get_kmer_str(seq_s)
                (self[seq_s_].ends if is_rc else self[seq_s_].starts).add(seq.id)

                seq_e_, is_rc = get_kmer_str(seq_e)
                (self[seq_e_].starts if is_rc else self[seq_e_].ends).add(seq.id)

    def filter_ends(self) -> None:
        end_kmers: set[str] = set()
        end_contigs: set[str] = set()
        for kmer, kc in self.kmers.items():
            if contigs := kc.is_end():
                end_kmers.add(kmer)
                for contig in contigs:
                    # remove contigs occur twice (both ends are single)
                    if contig in end_contigs:
                        self.isolate.add(contig)
                        end_contigs.remove(contig)
                    else:
                        end_contigs.add(contig)
        for kmer in end_kmers:
            self.kmers.pop(kmer)
        for kc in self.kmers.values():
            for contigs in (kc.starts, kc.ends):
                if discards := contigs & end_contigs:
                    kc.single.update(discards)
                    end_contigs.difference_update(discards)

    def export(self, file=""):
        import pandas as pd

        if not file:
            file = str(self.file) + ".kmer.tsv"

        pd.DataFrame(
            [self.kmers[kmer].report_dict() for kmer in sorted(self.kmers)]
        ).to_csv(file, sep="\t", index=False)

    def __getitem__(self, kmer: str):
        if kmer not in self.kmers:
            self.kmers[kmer] = KmerContig(kmer)
        return self.kmers[kmer]
