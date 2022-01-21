<!--
 * @Date: 2020-10-02 20:40:15
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-01-21 15:44:57
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

---
## [tmp](tmp/__init__.py)
- tmp files

## [seqPipe](seqPipe/__init__.py)
- scripts used in metagenomic analysis

## [setup](setup.py)
- clone: `git clone https://github.com/Hocnonsense/metaSC.git`
1.  install with repository and pip:
    - run `Hscripts]$``cd Python && python setup.py install`
    - to remove, just `pip uninstall Hscripts`
2.  use PYTHONPATH
    - add this in bash environment:
        `export PYTHONPATH=~/<your-path-to>/metaSC`

## functions
- [PyLibTool](PyLibTool/README.md)
- [Tool](tool/README.md)
- [BioTool](biotool/README.md)


# [***$\not$<!-- @Hwrn -->*~~`\`~~**](../README.md)
