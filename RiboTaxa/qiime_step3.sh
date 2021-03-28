#!/bin/bash
:<<!EOF!
 * @Editors: WangJing
 * @Date: 2020-12-30 21:03:03
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-03-28 22:59:22
 * @FilePath: /metaSC/RiboTaxa/qiime_step3.sh
 * @Description:
!EOF!


function show_usage (){
    echo ""
    echo "Usage: $(basename $0) table metadata tree subsampling reference [-o outpath] [group [...]]"
    echo ""
    echo 'Example$'" $(basename $0) result/01.2-demux-trimmed.qza \\ "
    echo "    results/Archaea/Archaea-table.qza results/00_metadata.tsv \\ "
    echo "    results/rooted-tree.qza 1000 \\ "
    echo "    -o result/Archaea -t 4"
    echo ""
    echo "Required parameters:"
    echo "    table | The ASV table file"
    echo "            File generated FROM qiime_step2"
    echo "    metadata | A tsv table listing the grouping method of each sample"
    echo "               The same in qiime_step2"
    echo "    tree | The phylogenetic tree of ASVs"
    echo "           File generated FROM qiime_step2"
    echo "    subsampling | Number of sequences per sample to be subsampled to"
    echo "           Evalued FROM qiime_step2"
    echo ""
    echo "Optional parameters:"
    echo "    -o outpath | The path to store the output files"
    echo "    -t threads | Threads if needed"
    echo "    group | Column names in the metadata that will test group significance on"
    echo ""
    echo "FILE output:"
    echo "   {outpath}/03.1-alpha_rarefaction.qzv"
    echo "       View and check alpha-diversity"
    echo ""
}


echo "$0 $*" >&2
set -- `getopt o:t:h "$@"`


while [ -n "$1" ]
do
    case "$1" in
     -h)
        show_usage
        exit ;;
     -o)
        outpath=$2
        mkdir -p $outpath
        echo "output files store in $outpath"
        shift ;;
    -t)
        threads=$2
        shift ;;
     --) ;;
     *)
        table=$1
        metadata=$2
        tree=$3
        subsample=$4
        break ;;
    esac
    shift
done

if [ -z $table ]; then
    echo "Error: table must be provided!"
    show_usage
    exit
fi
if [ -z $metadata ]; then
    echo "Error: metadata must be provided"
    show_usage
    exit
fi
if [ -z $tree ]; then
    echo "Error: tree must be provided"
    show_usage
    exit
fi
if [ -z $subsample ]; then
    echo "Error: subsampling must be provided"
    show_usage
    exit
fi

# Alpha rarefaction
#这里的--p-max-depth就写03_dada2_table.qzv里sequence count的最大值就行 (Frequency per feature->Maximum frequency, 如果报错则用报错给的值)
echo "Doing alpha-rarefaction"
qiime diversity alpha-rarefaction \
    --i-table $table \
    --i-phylogeny $tree \
    --p-max-depth $subsample \
    --m-metadata-file $metadata \
    --o-visualization $outpath/03.1-alpha_rarefaction.qzv || exit

# Diversity analysis
echo "Calculating core diversity metrices"
qiime diversity core-metrics-phylogenetic \
    --i-table $table \
    --i-phylogeny $tree \
    --p-sampling-depth $subsample \
    --m-metadata-file $metadata \
    --p-n-jobs-or-threads $threads \
    --output-dir $outpath/03-Diversity

# Statistics tests on alpha diversity
echo "Statistics tests on alpha diversities..."
function alpha_statistics (){
    qiime diversity alpha-group-significance \
        --i-alpha-diversity $outpath/03-Diversity/$index.qza \
        --m-metadata-file $metadata \
        --o-visualization $outpath/03-Diversity/$index-group-significance.qzv

    qiime diversity alpha-correlation \
        --i-alpha-diversity $outpath/03-Diversity/$index.qza \
        --m-metadata-file $metadata \
        --o-visualization $outpath/03-Diversity/$index-correlation.qzv
}

for index in faith_pd_vector observed_features_vector shannon_vector evenness_vector; do
  alpha_statistics
done

# Statistics tests on beta diversity
echo "Statistics tests on beta diversities..."
function beta_group_significance (){
    # qiime diversity beta-correlation \
        # --i-distance-matrix $outpath/03-Diversity/$index\_distance_matrix.qza\
        # --m-metadata-file $metadata \
        # --m-metadata-column $column \
        # --o-metadata-distance-matrix $outpath/03-Diversity/metadata_$column\_distance_matrix.qza\
        # --o-mantel-scatter-visualization $outpath/03-Diversity/$index-$coloumn-correlation.qzv

    qiime diversity beta-group-significance \
        --i-distance-matrix $outpath/03-Diversity/$index\_distance_matrix.qza\
        --m-metadata-file $metadata \
        --m-metadata-column $column \
        --p-pairwise \
        --o-visualization $outpath/03-Diversity/$index-$column-group-significance.qzv
}

while [ ! -z $6 ]; do
    column=$6
    shift
    for index in weighted_unifrac unweighted_unifrac jaccard bray_curtis; do
        beta_group_significance
    done
done

function beta_bioenv (){
    qiime diversity bioenv \
        --i-distance-matrix $outpath/03-Diversity/$index\_distance_matrix.qza\
        --m-metadata-file $metadata \
        --o-visualization $outpath/03-Diversity/$index-bioenv.qzv
}

for index in weighted_unifrac unweighted_unifrac jaccard bray_curtis; do
    beta_bioenv
done

echo "Done."
