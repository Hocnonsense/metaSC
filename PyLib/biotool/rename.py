# -*- coding: utf-8 -*-
"""
 * @Date: 2020-02-04 16:30:57
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-03 13:07:34
 * @FilePath: /HScripts/Python/mylib/biotool/rename.py
 * @Description: to rename sequence and store more messages if possible
        func:
            case '>':
                dict[i++] = readline()
                write '>'+i
            default:
                write line
            write dict in another file
"""


import re
# pylint: disable = anomalous-backslash-in-string


def save_name(path):
    """保存名字"""
    with path() as stdin, \
            path("namedict.txt", 'w') as outputname, \
            path("seqdumped.txt", 'w') as outputseq:
        dictnum = 0
        for line in stdin:
            if line.startswith(">"):
                outputname.write(str(dictnum) + ":" + line[1:])
                outputseq.write(">" + str(dictnum) + "\n")
                dictnum += 1
            else:
                outputseq.write(line)


def load_name(path, format_func):
    """读取名字"""
    with path("namedict.txt") as stdin:
        namedict = {}
        dictform = re.compile(r"(\d+):([^\n]+)")
        for i in stdin:
            nameitem = re.search(dictform, i)
            namedict[nameitem.group(1)] = format_func(nameitem.group(2))
    return namedict


def format_name(name, pattern=re.compile(r"\w+")):
    """名字正规化"""
    return "_".join(re.findall(pattern, name))


def rewrite_tree(path, filename, namedict):
    """重新绘制树"""
    savename = re.search(
        r"\\(\w+).\w+$", filename).group(1)   # "\\" in re is "\\\\", but for dir it is r"\"
    with path(filename) as stdin, path(savename + "_con.tre", 'w') as output:
        #
        treestr = stdin.readline()
        start = 0
        for i in re.finditer(r"[(,](\d+):", treestr):
            output.write(treestr[start:i.start(1)] + namedict[i.group(1)])
            start = i.end(1)
        output.write(treestr[start:])


if __name__ == "__main__":
    print("Warning: 缺少文件: seqdump.txt, 无法检验 `save_name` 函数")
    source_path = "demo\\degap"
    source_name = "seqdump.txt"
    file1 = __Path(source_name, source_path)
    # save_name(file1)
    formatnamedict = load_name(file1, format_name)
    contreefiles = [
        "Tree Inference\\191119081824\\degap_25.contree",
    ]
    for contree in contreefiles:
        rewrite_tree(file1, contree, formatnamedict)
