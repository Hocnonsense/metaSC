#!/bin/bash
:<<!EOF!
 * @Editors: WangJing
 * @Date: 2020-12-30 21:03:03
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-03-28 21:19:28
 * @FilePath: /metaSC/RiboTaxa/qiime_step1.sh
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
primer_f="GTGYCAGCMGCCGCGGTAA"
primer_r="GGACTACNVGGGTWTCTAAT"
threads=4


echo "$0 $*" >&2
set -- `getopt o:p:t:h "$@"`


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
    -p)
        primer_f=${2%%,*}
        primer_r=${2##*,}
        echo "using primer: f: $primer_f"
        echo "              r: $primer_r"
        shift ;;
    -t)
        threads=$2
        shift ;;
    --) ;;
    *)
        manifest=$1
        break ;;
    esac
    shift
done

if [ -z $manifest ]; then
    echo "Error: manifest must be provided!"
    show_usage
    exit
fi


# Import raw fastq files
echo "Importing FASTQ files..."
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path $manifest \
    --input-format PairedEndFastqManifestPhred33V2 \
    --output-path $outpath/01.1-demux.qza

# Cut PCR primer from the sequences
echo "Cutting PCR primers..."
qiime cutadapt trim-paired \
    --i-demultiplexed-sequences $outpath/01.1-demux.qza \
    --p-front-f $primer_f \
    --p-front-r $primer_r \
    --p-cores $threads \
    --o-trimmed-sequences $outpath/01.2-demux-trimmed.qza

# Summarize the trimmed sequences
echo "Summarizing sequence quality..."
qiime demux summarize \
    --i-data $outpath/01.2-demux-trimmed.qza \
    --o-visualization $outpath/01.2-demux-trimmed.qzv

echo "Done."
