#!/bin/bash
:<<!EOF!
 * @Editors: WangJing
 * @Date: 2020-12-30 21:03:03
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-03-28 23:14:23
 * @FilePath: /metaSC/RiboTaxa/qiime_step0.sh
 * @Description:
!EOF!


function show_usage (){
    echo ""
    echo "Usage: $(basename $0) manifest [-o outpath] [-p forward_primer,reverse_primer]"
    echo ""
    echo 'Example$'" $(basename $0) result/manifest.tsv \\ "
    echo "    -o result -p GTGYCAGCMGCCGCGGTAA,GGACTACNVGGGTWTCTAAT"
    echo ""
    echo "Required arguments:"
    echo "    manifest | A tsv table listing the fastq files for each sample"
    echo "               In the file there is three cols: "
    echo "                   sample-id	forward-absolute-filepath	reverse-absolute-filepath"
    echo "               Absolute path is recommended"
    echo "     * sequenced samples"
    echo "        sequenced samples ARE GIVEN in FILE manifest"
    echo "        However, file IN fastq format SHOULD exist IN given path"
    echo ""
    echo "Optional arguments:"
    echo "    -o outpath | Path to store the output files"
    echo "    -t threads | Threads if needed"
    echo "    -p primer | Pair of primers joined with ','"
    echo "        forward_primer | Forward PCR primer, default: GTGYCAGCMGCCGCGGTAA"
    echo "        reverse_primer | Forward PCR primer, default: GGACTACNVGGGTWTCTAAT"
    echo ""
    echo "FILE output:"
    echo "    * {outpath}/01.1-demux.qza"
    echo "   {outpath}/01.2-demux-trimmed.qza"
    echo "       The qza file is used FOR qiime_step2"
    echo "   {outpath}/01.2-demux-trimmed.qzv"
    echo "       The qzv file is used TO check length of reads FOR qiime_step2"
    echo ""

    return 0
}
sequence=sequences.qza
ref=ref-taxonomy.qza
label=''
threads=4


echo "$0 $*" >&2
set -- `getopt o:s:r:l:t:h "$@"`


while [ -n "$1" ]
do
    case "$1" in
    -h)
        show_usage
        exit ;;
    -o)
        outpath=$2
        mkdir -p $outpath
        echo "output files store in"
        echo "    $(cd $(dirname $outpath); pwd)"/"$(basename $outpath)"
        shift ;;
    -s)
        sequence=$2
        shift ;;
    -r)
        ref=$2
        shift ;;
    -l)
        label="-$2"
        shift ;;
    -t)
        threads=$2
        shift ;;
    --) ;;
    *)
        primer_f=${1%%,*}
        primer_r=${1##*,}
        echo "using primer: f: $primer_f"
        echo "              r: $primer_r"
        break ;;
    esac
    shift
done

if [ -z $manifest ]; then
    echo "Error: manifest must be provided!"
    show_usage
    exit
fi


qiime feature-classifier extract-reads \
    --i-sequences $sequences \
    --p-f-primer $primer_f \
    --p-r-primer $primer_r \
    --o-reads $outpath/ref_seqs${label}.qza \
    --p-n-jobs $threads \
    --verbose

qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads $outpath/ref_seqs-${label}.qza \
        --i-reference-taxonomy $ref \
        --o-classifier $outpath/classifier${label}.qza \
        --verbose

echo "Done."
