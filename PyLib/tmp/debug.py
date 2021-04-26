#!/bin/python
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -o Out-err/%x-%j.out
#SBATCH -e Out-err/%x-%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
# -*- coding: utf-8 -*-
"""
    @Date: 2020-01-20 23:12:59
    @LastEditors: Hwrn
    @LastEditTime: 2020-09-21 21:26:31
    @FilePath: /mylib/debug.py
    @Description: to do checkm lineage_wf, select, checkmToR togather
"""


import os
from sys import path, argv, stderr, stdout
from _thread import start_new_thread


def run_checkmToR(ckm_dir, binID, UID, mdf_dir):
    for line in os.popen(
        "python /lustre/home/acct-ioozy/ioozy/Hwrn/Scripts/Python/checkMarkToR.py"
        + " " + ckm_dir
        + " " + binID
        + " " + UID
        + " " + mdf_dir
    ):
        print(line, file=stderr)
    for line in os.popen(
        "cp"
        + " " + os.path.join(bin_dir, binID+".fa")
        + " " + os.path.join(mdf_dir, binID+".fa")
    ):
        print(line, file=stderr)


__doc = """
    useage: sbatch checkMPipe.py bin_dir out_dir

    args:
        bin_dir:    the dir of the raw bin to be checked.
        out_dir:    the dir of checkm files output.
                    checkM resutls saved in out_dir/checkm
                    bins to modify and their marker_gene.csv are
                    saved in out_dir/toR
"""

if len(argv) == 1:
    print(__doc)
    exit(0)

# check: ckm_dir = "/lustre/home/acct-ioozy/ioozy/Hwrn/Pipeline/2020-09-MgAffact/F-06-MAG/02_checkm/"
bin_dir, ckm_dir = argv[1], argv[2]
# check: ckmout = "/lustre/home/acct-ioozy/ioozy/Hwrn/Pipeline/2020-09-MgAffact/Out-err/06.5-CheckM-1715533.out"
ckmout = os.path.join(ckm_dir, "checkm.out")
# check: mdf_dir = "/lustre/home/acct-ioozy/ioozy/Hwrn/Pipeline/2020-09-MgAffact/F-06-MAG/03_modify/0_raw/"
mdf_dir = os.path.join(ckm_dir, "modify")
os.makedirs(ckm_dir, exist_ok=True)
os.makedirs(mdf_dir, exist_ok=True)

shstdout = os.popen(
    "checkm lineage_wf -x fa -t 39 "
    + bin_dir + " " + ckm_dir
)

with open(ckmout, "w") as sout:
    for line in shstdout.readlines():
        print(line, file=stdout)
        sout.write(line)


# read ckmout
ckmtable = []
with open(ckmout) as fin:
    for line in fin:
        print(line)
        if line.strip().split()[0] == "Bin":
            break
    else:
        print("""
            we read the file like:
                ```
                [2020-09-19 16:44:36] INFO: CheckM v1.1.2
                [2020-09-19 16:44:36] INFO: checkm lineage ...
                ...                                        ...
                [2020-09-19 17:33:36] INFO: Parsing HMM hi ...
                ------------------------------------------ ...
             ->   Bin Id                            Marker ...
                ------------------------------------------ ...
                  metabat2_90_60.124          o__Cytophaga ...
                ```
            and we just reach the line of `->` table head

            So, what's your problem?
            """, file=stderr)
    fin.readline()  # read the secend line between table body and head
    for line in fin:
        if line[0] == "-":
            break
        values = line.strip().split()
        ckmtable.append([
            # Bin Id, Marker lineage (UID), Completeness, Contamination
            values[0], values[2][1:-1], float(values[-3]), float(values[-2])
        ])
    else:
        print("""
            the table will end with a line of "-"

            So, what's your problem?
            """, file=stderr)

print("OK, now, find which are good and which should be modified:", file=stderr)
goodbins = []
mdf_bins = []
# Bin Id, Marker lineage (UID), Completeness, Contamination, Strain heterogeneity
for binID, UID, Compl, Contai in ckmtable:
    if Compl > 70:
        if Contai < 10:
            goodbins.append(binID)
        else:
            mdf_bins.append(binID)
            start_new_thread(run_checkmToR,
                             (ckm_dir, binID, UID, mdf_dir))
else:
    print("""
        Now, you get :
            {} good bins
            {} bins to modify
            {} bins useless
        """.format(len(goodbins), len(mdf_bins),
                   len(ckmtable)-len(goodbins)-len(mdf_bins)),
          file=stderr)

print("Output list to txt")
with open("l-06-checkm-good.txt", "w") as gout \
        , open("l-06-checkm-good.txt", "w") as mout \
        :
    gout.write("\n".join(goodbins))
    mout.write("\n".join(mdf_bins))
