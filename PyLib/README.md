<!--
 * @Date: 2020-10-02 20:40:15
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-19 19:42:45
 * @FilePath: /metaSC/PyLib/README.md
 * @Description:
-->
Python
===
Python files

- This package contains all python scripts, but not include any file or project that use python with other languages.
- This package is handled for 2 reasons:
    1.  Keep some useful functions.
        - These functions should be wrapped, with documents, typing, and tests
        - these functions include:
            - help me to handle other functions in this package (PyLib): [[PyLibTool/README]]
            - python API (glue water) for other program, e.g. output tsv formatting: reader
            - Python scripts that can get command parameters and called by bash
            - other useful functions: biotool
    2.  Demo to show example of other package
        - tool

- to test, run `pytest --cov --cov-report=html -m "not slowdownload"`

---
## [tmp](tmp/__init__.py)
- tmp files

## [seqPipe](seqPipe/__init__.py)
- scripts used in metagenomic analysis

## [setup](setup.py)
- clone: `git clone https://github.com/Hocnonsense/metaSC.git`
1.  install with repository and pip:
    - run `metaSC]$``python PyLib/setup.py develop`
    - to remove, just `pip uninstall metaSC`
2.  use PYTHONPATH
    - add this in bash environment:
        `export PYTHONPATH=~/<your-path-to>/metaSC`

## before setup
- please run `stubgen --output . PyLib/`

## functions
- [PyLibTool](PyLibTool/README.md)
- [Tool](tool/README.md)
- [BioTool](biotool/README.md)


## changelog
- 0.0.5
    - add tool.wildcards
- 0.0.4
    - update biotool.download
    - add retries_download
    - add RetriveUrl.download_genome
    - add RefSeqURL.download_gff
    - update RetriveUrl.download_fna to RetriveUrl.download_genome_to
    - update download_fna to download_genome

# [***$\not$<!-- @Hwrn -->*~~`\`~~**](../README.md)
