<!--
 * @Date: 2022-01-21 11:26:36
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-01-21 15:20:47
 * @FilePath: /metaSC/PyLib/PyLibTool/README.md
 * @Description:
-->
PyLibTool
===

- Useful functions that tightly coupled with the structure and functions of this library

---
## file_info

- extract_doc: To extract metadata from `__doc__` of python file (or any `korofildheader-generated` fileheader)
- verbose_import:
    - Some files are updated unexpectedly, and scripts related to these file are affected.
    - This function helps to report the version of the file.
    - It use `logging.logger` to report -- and just return the logger

## tmpPkl
- use `pickle` to store temporary files

- TmpPkl: a wrapable class to avoid running a time/io comsuming function too often

## demo_exec
- a demo of executing python file


# [***$\not$<!-- @Hwrn -->*~~`\`~~**](README.md)
