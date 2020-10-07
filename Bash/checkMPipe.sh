#!/bin/bash
#SBATCH -p cpu
#SBATCH --mail-type=end
#SBATCH --mail-user=Hwrn.aou@sjtu.edu.cn
#SBATCH -o Out-err/%x-%j.out
#SBATCH -e Out-err/%x-%j.err
#SBATCH -N 1
#SBATCH --ntasks-per-node=40
:<<!EOF!
 * @Date: 2020-09-20 15:33:13
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-06 21:56:41
 * @FilePath: /HScripts/Bash/checkMPipe.sh
 * @Description:
!EOF!
set -e && echo "$0 $*" >&2 && source ~/.bashrc


if (( $# == 0 ));
then
    cat >&2 <<!EOF!
    useage: sbatch checkMPipe.py bin_dir tor_dir

    args:
        bin_dir:    the dir of the raw bin to be checked.
        tor_dir:    the dir of checkm files output.
            tor_dir/checkm:
                    checkM resutls saved in
            tor_dir/modify:
                    bins to modify and their marker_gene.csv
            tor_dir/checkm.out:
                    checkm output
            tor_dir/checkm-good.list: good bins
            tor_dir/checkm-modify.list: bins can be modified

    *tmp files:
        ./checkm-good.tmp
        ./checkm-modify.tmp
        ${skmout}
        ./checkm-to-r.tmp

!EOF!
    exit 1
fi


echo "your args are:"
for var in $*; do echo "$var" >&2; done


conda activate python36 || echo "Does your environment has checkm?"

bin_dir=$1
tor_dir=$2
ckm_dir="${tor_dir}/checkm"
mdf_dir="${tor_dir}/modify"
ckm_out="${tor_dir}/checkm.out"

mkdir ${tor_dir} || ls -l ${tor_dir} >&2
mkdir ${ckm_dir} || ls -l ${ckm_dir} >&2
mkdir ${mdf_dir} || ls -l ${mdf_dir} >&2

if (( $# == 3 )) && [ -r ${ckm_dir} ]
then
    echo "using exist checkm results: ${ckm_out}" >&2
    ckm_out=$3
else
    echo "checkm stdout > ${ckm_out}" >&2
    checkm lineage_wf \
        -x fa -t 39 \
        ${bin_dir} ${ckm_dir} \
    > ${ckm_out}

    echo "checkm finished" >&2
fi

skmout=./skmout.tmp
cp ${ckm_out} ${skmout} \
|| exit 1

python <<!EOF!
from collections import OrderedDict
from sys import stderr

ckmap = OrderedDict()
with open("skmout.tmp") as fin:
    for line in fin:
        print(line)
        if line.strip().split()[0] == "Bin":
            break
    else:
        print("""
            we read the file like:
                \`\`\`
                [2020-09-19 16:44:36] INFO: CheckM v1.1.2
                [2020-09-19 16:44:36] INFO: checkm lineage ...
                ...                                        ...
                [2020-09-19 17:33:36] INFO: Parsing HMM hi ...
                ------------------------------------------ ...
             ->   Bin Id                            Marker ...
                ------------------------------------------ ...
                  metabat2_90_60.124          o__Cytophaga ...
                \`\`\`
            and we just reach the line of \`->\` table head

            So, what's your problem?
            """, file=stderr)
        exit(1)
    fin.readline()  # read the secend line between table body and head
    for line in fin:
        if line[0] == "-":
            break
        values = line.strip().split()
        # Bin Id: [Marker lineage (UID), Completeness, Contamination]
        ckmap[values[0]]=[
            values[2][1:-1], float(values[-3]), float(values[-2])
        ]
    else:
        print("""
            the table will end with a line of "-"

            So, what's your problem?
            """, file=stderr)

print("OK, now, find which are good and which should be modified:", file=stderr)
goodbins = []
mdf_bins = []
uid_nums = []
# Bin Id, Marker lineage (UID), Completeness, Contamination
for binID, (UID, Compl, Contai) in ckmap.items():
    if Compl > 70:
        if Contai < 10:
            goodbins.append(binID)
        else:
            mdf_bins.append(binID)
            uid_nums.append(UID)
else:
    print("""
        Now, you get:
            {} good bins
            {} bins to modify
            {} bins useless
        more information:
            {}
        """.format(len(goodbins),
                   len(mdf_bins),
                   len(ckmap)-len(goodbins)-len(mdf_bins),
                   str(ckmap)),
          file=stderr)

print("Output list to txt", file=stderr)
with open("checkm-good.tmp", "w") as gout, open("checkm-modify.tmp", "w") as mout, open("checkm-to-r.tmp", "w") as uout:
    gout.write("\n".join(goodbins)+"\n")
    mout.write("\n".join(mdf_bins)+"\n")
    uout.write("\n".join([" ".join(match)
                          for match in zip(mdf_bins, uid_nums)])+"\n")
!EOF!

cat ./checkm-to-r.tmp | while read match
do
{
    echo "${match} start" >&2
    python \
        ~/Users/Hwrn/Scripts/Python/checkMarkToR.py \
        ${ckm_dir} \
        ${match} \
        ${mdf_dir} \
    || echo "${match} fault" >&2
    echo "${match} done" >&2
} &
done

wait
echo "checkmToR ended" >&2

echo "move tmp files to ${tor_dir}" >&2
mv ./checkm-good.tmp "${tor_dir}/checkm-good.list"
mv ./checkm-modify.tmp "${tor_dir}/checkm-modify.list"
rm -f ${skmout}
rm -f ./checkm-to-r.tmp

echo "exit successfully" >&2
