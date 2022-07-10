# -*- coding: utf-8 -*-
"""
 * @Date: 2022-02-21 10:25:20
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-02-21 12:47:34
 * @FilePath: /metaSC/PyLib/biotool/fatree.py
 * @Description:
    Utils to build tree
"""
import os
from Bio import SeqIO, AlignIO
from PyLib.tool.shell import runsh_safe


def mafft_align_genes(fa_file: str):
    """
    ------------------------------------------------------------------------------
      MAFFT v7.487 (2021/Jul/25)
      https://mafft.cbrc.jp/alignment/software/
      MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
    ------------------------------------------------------------------------------
    High speed:
      % mafft in > out
      % mafft --retree 1 in > out (fast)

    High accuracy (for <~200 sequences x <~2,000 aa/nt):
      % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)
      % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out)
      % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out)

    If unsure which option to use:
      % mafft --auto in > out

    --op # :         Gap opening penalty, default: 1.53
    --ep # :         Offset (works like gap extension penalty), default: 0.0
    --maxiterate # : Maximum number of iterative refinement, default: 0
    --clustalout :   Output: clustal format, default: fasta
    --reorder :      Outorder: aligned, default: input order
    --quiet :        Do not report progress
    --thread # :     Number of threads (if unsure, --thread -1)
    --dash :         Add structural information (Rozewicki et al, submitted)
    """
    os.system(
        f"mafft "
        f"    --quiet "
        f"    --dash "
        f"    --auto "
        f"    {fa_file} "
        f"> {fa_file}.afa"
    )
    return f"{fa_file}.afa"


def muscle_clust_genes(fa_file: str):
    """
    MUSCLE v3.8.1551 by Robert C. Edgar

    http://www.drive5.com/muscle
    This software is donated to the public domain.
    Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.


    Basic usage

        muscle -in <inputfile> -out <outputfile>

    Common options (for a complete list please see the User Guide):

        -in <inputfile>    Input file in FASTA format (default stdin)
        -out <outputfile>  Output alignment in FASTA format (default stdout)
        -diags             Find diagonals (faster for similar sequences)
        -maxiters <n>      Maximum number of iterations (integer, default 16)
        -maxhours <h>      Maximum time to iterate in hours (default no limit)
        -html              Write output in HTML format (default FASTA)
        -msf               Write output in GCG MSF format (default FASTA)
        -clw               Write output in CLUSTALW format (default FASTA)
        -clwstrict         As -clw, with 'CLUSTAL W (1.81)' header
        -log[a] <logfile>  Log to file (append if -loga, overwrite if -log)
        -quiet             Do not write progress messages to stderr
        -version           Display version information and exit
    """
    runsh_safe(
        f"muscle " f"    --quiet " f"    -in {fa_file} " f"    -out {fa_file}.afa"
    )
    return f"{fa_file}.afa"
