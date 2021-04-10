<!--
 * @Date: 2021-04-08 20:44:28
 * @LastEditors: Hwrn
 * @LastEditTime: 2021-04-10 14:18:08
 * @FilePath: /metaSC/RiboTaxa/README.md
 * @Description:
-->
QIIME2 16S RNA diversity analyzing helper
===

- use `2>&1 |tee qiime.log`

---
## qiime_step0
- This script is used to generate classifier for
    - any primer pairs
    - and given seqs & ref-tax file
        - If seqs & ref-tax file should be changed, please re-run `qiime tools import` to generate new files

- output file:
    - ref_seqs.qza
    - classifier.qza
- To refer:
    - `label='-`
    - `ln -s ref_seqs${label}.qza your_dir/00_data/`
    - `ln -s classifier${label}.qza your_dir/00_data/`

## qiime_step1
- To refer:
    - `ln -s results.A/01.2-demux-trimmed.qzv Analyze/01.2-demux-trimmed.A.qzv`
    - `ln -s results.B/01.2-demux-trimmed.qzv Analyze/01.2-demux-trimmed.B.qzv`

## qiime_step2
- To refer:
    - `ln -s results.A/02.1-dada2_table.qzv Analyze/02.1-dada2_table.A.qzv`
    - `ln -s results.B/02.1-dada2_table.qzv Analyze/02.1-dada2_table.B.qzv`
    - `ln -s results.A/Archaea/02.4-Archaea_table_taxa.qzv Analyze/02.4-table_taxa.A.qzv`
    - `ln -s results.B/Bacteria/02.4-Bacteria_table_taxa.qzv Analyze/02.4-table_taxa.B.qzv`

## qiime_step3
- To refer:
    - `ln -s results.A/Archaea/03.1-alpha_rarefaction.qzv Analyze/03.1-alpha_rarefaction.A.qzv`
    - `ln -s results.B/Bacteria/03.1-alpha_rarefaction.qzv Analyze/03.1-alpha_rarefaction.B.qzv`


# [***$\not$<!-- @Hwrn -->*~~`\`~~**](README.md)
