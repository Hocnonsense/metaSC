#!/bin/bash
:<<!EOF!
 * @Editors: WangJing
 * @Date: 2020-12-30 21:03:03
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-08-16 15:14:04
 * @FilePath: /metaSC/RiboTaxa/qiime_step2.sh
 * @Description:
!EOF!


function show_usage (){
    echo ""
    echo "Usage: $(basename $0) seqs metadata length_forward,length_reverse [-p trim_forward,trim_reverse] [-o outpath] [-f reference]"
    echo ""
    echo 'Example$'" $(basename $0) result/01.2-demux-trimmed.qza \\ "
    echo "    result/00_metadata.tsv 200,180 \\ "
    echo "    -o result -f 00_data/classifier.qza -p 0,0 -t 4"
    echo ""
    echo "Required arguments:"
    echo "    seqs | Demultiplexed sequence file"
    echo "           File generated FROM qiime_step1"
    echo "    metadata | A tsv table listing the grouping method of each sample"
    echo "               For example:"
    echo "                   #SampleID	BarcodeSequence	LinkerPrimerSequence"
    echo "               #SampleID SHOULD be the same as IN manifest file"
    echo "    length_forward,length_reverse | The required length OF forward and reverse sequences"
    echo "                                    These two number IS evalued from qiime_step1"
    echo ""
    echo "Optional arguments:"
    echo "    -o outpath | Path to store the output files"
    echo "    -t threads | Threads if needed"
    echo "    -f reference | The reference for taxonomy classification, default: reference for 515F-806R"
    echo "                   This is generated by qiime_step0 where a SILVA reference is accessible"
    echo "    -p trim_forward,trim_reverse | The first nucleotides that should be trimmed from forward and reverse sequences"
    echo "                                   Default IS 0,0"
    echo ""
    echo "FILE output:"
    echo "   {outpath}/02.1-dada2_table.qzv"
    echo "       View and check total reads detected"
    echo "   {outpath}/{domain}/02.4-{domain}_table.qza"
    echo "       Evalue 'Feature Count' FOR beta diversity OF given domain IN qiime_step3"
         # 可以取 Interactive Sample Detail -> Feature Count 下所需的值进行 Beta 多样性对比
    echo "   {outpath}/02.2-rooted_tree.qza &"
    echo "   {outpath}/02.3-taxonomy.qza &"
    echo "   {outpath}/{domain}/02.4-{domain}_table.qza"
    echo "       Used in qiime_step3"
    echo "   {outpath}/02.3-dada2_table_taxa.qzv &"
    echo "   {outpath}/{domain}/02.4-{domain}_table_taxa.qzv"
    echo "       View alpha diversity"
    echo ""
}
trim_f=0
trim_r=0
reference="/software/DataBase/qiime2/silva-138-99-amplicon-classifier.qza"


echo "$0 $*" >&2
set -- `getopt o:p:t:f:h "$@"`


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
        trim_f=${2%%,*}
        trim_r=${2##*,}
        echo "using trim-left: f: $trim_f"
        echo "                 r: $trim_r"
        shift ;;
    -t)
        threads=$2
        shift ;;
     -f)
        reference=$2
        shift ;;
     --) ;;
     *)
        seqs=$1
        metadata=$2
        len_f=${3%%,*}
        len_r=${3##*,}
        break ;;
    esac
    shift
done

if [ -z $seqs ]; then
    echo "Error: seqs must be provided!\n"
    show_usage
    exit
fi
if [ -z $metadata ]; then
    echo "Error: metadata must be provided\n"
    show_usage
    exit
fi
if [ -z $len_f ]; then
    echo "Error: length for trunc must be provided\n"
    show_usage
    exit
fi
if [ -z $outpath ]; then
    outpath="./"
    echo "output files store in"
    echo "    `pwd`"
fi

# for stupid feature-classifier
export TMPDIR=$outpath/TMP
mkdir -p $TMPDIR

# Denoise with DADA2
# 数据去噪，生成OTU，--p-trunc-len-f是指把R1里的reads切到多长，这个值不要太高。因为短于这个值的序列会被删掉，让你本来测了很多数据，结果一下子只用了2%...--p-trunc-len-r同理。
# --p-trim-left-f是切除R1前端的，看你的数据情况，如果质量都挺高的，就留着，写0
echo "Denoising..."
qiime dada2 denoise-paired \
    --i-demultiplexed-seqs $seqs \
    --p-trunc-len-f $len_f \
    --p-trunc-len-r $len_r \
    --p-trim-left-f $trim_f \
    --p-trim-left-r $trim_r \
    --p-pooling-method pseudo \
    --p-n-threads $threads \
    --o-table $outpath/02.1-dada2_table.qza \
    --o-representative-sequences $outpath/02.1-dada2_rep_seq.qza \
    --o-denoising-stats $outpath/02.1-dada2_stats.qza

# Tabulate the rep seqs
qiime feature-table tabulate-seqs \
    --i-data $outpath/02.1-dada2_rep_seq.qza \
    --o-visualization $outpath/02.1-dada2_rep_seq.qzv

# Summarize the ASV table
#对dada2的结果可视化，看一下这一步的结果，看看你得到了多少OTU，用到了多少序列
echo "Summarizing ASV table..."
qiime feature-table summarize \
    --i-table $outpath/02.1-dada2_table.qza \
    --m-sample-metadata-file $metadata \
    --o-visualization $outpath/02.1-dada2_table.qzv

# Visualize the denoising stats
echo "Summarizing denoise stats"
qiime metadata tabulate \
    --m-input-file $outpath/02.1-dada2_stats.qza \
    --o-visualization $outpath/02.1-dada2_stats.qzv

# Make phylogenetic tree
echo "Making phylogeny..."
qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences $outpath/02.1-dada2_rep_seq.qza \
    --p-n-threads $threads \
    --o-alignment $outpath/02.2-aligned_dada2_rep_seq.qza \
    --o-masked-alignment $outpath/02.2-masked_aligned_dada2_rep_seq.qza \
    --o-tree $outpath/02.2-unrooted_tree.qza \
    --o-rooted-tree $outpath/02.2-rooted_tree.qza

# Taxonomic classification
# 将物种信息已被放入 02.3-taxonomy.qza, 接下来据此再次清除序列
echo "Classifying ASVs..."
qiime feature-classifier classify-sklearn \
    --i-reads $outpath/02.1-dada2_rep_seq.qza \
    --i-classifier $reference \
    --p-reads-per-batch 10000 \
    --p-n-jobs $threads \
    --o-classification $outpath/02.3-taxonomy.qza

# Report the taxonomy assignment
qiime metadata tabulate \
    --m-input-file $outpath/02.3-taxonomy.qza \
    --o-visualization $outpath/02.3-taxonomy.qzv

# Make taxonomy barplot
# 显示分组后的 barplot
qiime taxa barplot \
    --i-table $outpath/02.1-dada2_table.qza \
    --i-taxonomy $outpath/02.3-taxonomy.qza \
    --m-metadata-file $metadata \
    --o-visualization $outpath/02.3-dada2_table_taxa.qzv

for domain in Archaea Bacteria
do
    # Split the table to group
    echo "Extracting $domain..."
    mkdir -p $outpath/$domain

    qiime taxa filter-table \
        --i-table $outpath/02.1-dada2_table.qza \
        --i-taxonomy $outpath/02.3-taxonomy.qza \
        --p-include d_0__$domain \
        --o-filtered-table $outpath/$domain/02.4-${domain}_table.qza

    qiime feature-table summarize \
        --i-table $outpath/$domain/02.4-${domain}_table.qza \
        --m-sample-metadata-file $metadata \
        --o-visualization $outpath/$domain/02.4-${domain}_table.qzv

    qiime taxa barplot \
        --i-table $outpath/$domain/02.4-${domain}_table.qza \
        --i-taxonomy $outpath/02.3-taxonomy.qza \
        --m-metadata-file $metadata \
        --o-visualization $outpath/$domain/02.4-${domain}_table_taxa.qzv
done

echo "Done."
