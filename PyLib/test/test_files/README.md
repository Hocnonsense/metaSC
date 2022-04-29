<!--
 * @Date: 2022-04-29 11:02:00
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 11:10:18
 * @FilePath: /metaSC/PyLib/test/test_files/README.md
 * @Description:
-->

===

---
##
```log
(python39) [hwrn@localhost] -- 2022-04-29 11:07:52 -- (metaSC)
~/software/metaSC/PyLib/test/test_files$ python
Python 3.9.10 | packaged by conda-forge | (main, Feb  1 2022, 21:24:11)
[GCC 9.4.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> genome = "GCF_000215995"
>>> from PyLib.biotool.download import download_fna
>>> download_fna(genome)
2022-04-29 11:08:59  download.py : INFO  download to GCF_000215995.fna
2022-04-29 11:09:05  download.py : INFO  Downloading "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/215/995/GCF_000215995.1_ASM21599v1/GCF_000215995.1_ASM21599v1_genomic.fna.gz" to "GCF_000215995.fna.gz"

2022-04-29 11:09:06  download.py : INFO  Downloading file of size: 0.49 MB

2022-04-29 11:09:07  download.py : INFO

PosixPath('GCF_000215995.fna')
>>>
```

# [***$\not$<!-- @Hwrn -->*~~`\`~~**](README.md)
