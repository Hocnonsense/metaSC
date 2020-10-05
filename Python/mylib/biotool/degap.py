# -*- coding: utf-8 -*-
# pylint: disable = anomalous-backslash-in-string
"""
 * @Date: 2020-02-04 15:41:48
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-03 13:04:57
 * @FilePath: /MyScripts/Python/mylib/biotool/degap.py
    @Description:
        读取 CLUSTAL O(1.2.4) multiple sequence alignment 文件
        按阈值删除碱基

        以以下文件格式为标准
__title__ \
\\s+
[(__name__\\s+[ATGC-]+\\s[\\d]+\\n)+ \
 [__name__\\s+[ATGC-]+\\n]? \
 \\s+]+
\\.

        经过讨论, 将ID作为key加入字典, 并append序列直至读完, 无需关心其他
        其他可变的位置:
            开头的'>', 序列后的数字
   @TODU: 接下来, 最好能加入检查处理后序列是否完全相同的部分
"""


import re
import numpy
from typing import Tuple


def report_gap(seqtable: numpy.ndarray):
    """
        @description:
            通过画图的方式显示位点 gap 数量的分布 (在同源的 n 条序列中)
            三个 `*` 的横坐标分别是 25%, 50%, 75%; 红线+蓝色条形图
    """
    import matplotlib.pyplot as plt
    seqnum = len(seqtable)
    length = len(seqtable[0])
    reporttable = [[] for i in range(seqnum + 1)]
    for i in range(length):
        row = seqtable[:, i]
        reporttable[sum(row != '-')].append(row)
    fig1 = [len(i) for i in reporttable]
    fig = [[0, fig1[0], fig1[0]]]
    for i in range(1, seqnum + 1):
        # if fig1[i] != 0:
        fig.append([i, fig1[i], fig[-1][2] + fig1[i]])
    fig = numpy.array(fig).T

    _, ax1 = plt.subplots()
    plt.axis([0, seqnum + 1, 0, length + 1])
    ax1.bar(fig[0], fig[1])
    middle = max(fig[1]) / 2
    ax1.scatter([seqnum * 0.25, seqnum * 0.5, seqnum * 0.75],
                [middle, middle, middle], marker="*")   # ruler of 25%, 50%, 75%
    ax1.set_xlabel("0 gap <-- " + str(seqnum // 2) + " gap --> " + str(seqnum) + " gap")
    ax1.set_ylabel("number of sites of same gap number")
    ax1.set_ylim(0, max(fig[1]) + 1)

    ax2 = ax1.twinx()
    ax2.plot(fig[0], fig[2], "-", color="r", linewidth=1)
    ax2.set_ylim(0, fig[2][-1] + 1)
    ax2.set_ylabel("sum of sites")

    plt.show()


def delete_gap(seqtable: numpy.ndarray, DEGAP_THRESHOLD: float) -> numpy.ndarray:
    """
        @description:
            删除节点
        @param DEGAP_THRESHOLD:
            需要去除的比例,
            若该值 > 1, 则保留比该值更多序列在此非 gap 的位点.
            如: delete_gap(array(
                    [["-", "A", "G"], ["-", "-", "G"], ["A", "A", "G"]], 2)
                -> 则仅删除第一个位点
        @return 删除的节点的比例 (输入值不一定 < 1), 处理后的序列
    """
    leaves = len(seqtable)
    if DEGAP_THRESHOLD > 1:
        DEGAP_THRESHOLD = DEGAP_THRESHOLD / leaves
    newseqtable = []
    for i in range(len(seqtable[0])):
        row = seqtable[:, i]
        if(sum(row != '-') / leaves > DEGAP_THRESHOLD):  # 计数
            newseqtable.append(row)
    return DEGAP_THRESHOLD, numpy.array(newseqtable).T   # 返回新的DNA序列


def clean_seq(sequence: list):
    return "".join([i if i in "ATCGatcg-" else "" for i in sequence])


def FORMAT(IOtext, *args: Tuple["namelist", "seqarray", "MAX_SEQUENCE_LEN"]) -> dict:
    """the format for read and write the sequence"""
    try:
        namelist, seqarray, MAX_SEQUENCE_LEN = args
    except ValueError:  # if no args, that is input
        seqdict = {}
        args = ""
        for line in IOtext:
            name, sequence = line
            seqdict[name] = seqdict.get(name, "") + sequence
        return seqdict
        # then format_input will change seqdict to (namedict:list, seqarray:numpy.ndarray)
    else:   # is for writing
        # for max sequence length in a line
        maxseqlength = len(seqarray[0]) // MAX_SEQUENCE_LEN
        for name, seq in zip(namelist, seqarray):
            IOtext.write(name)
            for i in range(maxseqlength):
                seqpart = seq[i * MAX_SEQUENCE_LEN:(i + 1) * MAX_SEQUENCE_LEN]
                # should change list to string, and go to the next line
                IOtext.write("".join(seqpart) + "\n")
            # the remain part if there are
            seqpart = seq[maxseqlength * MAX_SEQUENCE_LEN:]
            IOtext.write("".join(seqpart) + "\n")


def CLUSTAL(IOtext, *args: Tuple["namelist", "seqarray", "MAX_SEQUENCE_LEN"]) -> dict:
    """
        @description:
            本函数兼处理输入和输出
            read format like:
            __title__\\n
            \\s+[(__name__\\s+[ATGC-]+\\s[\\d]+\\n)+[__name__\\s+[ATGC-]+\\n]?\\s+]+
            \\.
        @param args:
            若无, 则读入文件, 输出 seqdict {name: sequence}
            否则写入文件
    """
    try:
        namelist, seqarray, MAX_SEQUENCE_LEN = args
    except ValueError:  # if no args, that is input
        seqdict = {}
        detect = re.compile(r'^(\w+)\s+([agtcAGTC-]+)\s')
        for line in IOtext:
            express = re.search(detect, line)
            if(express):
                name, sequence = express.group(1), express.group(2)
                seqdict[name] = seqdict.get(name, "") + sequence
        return seqdict
    else:
        maxseqlength = len(seqarray[0]) // MAX_SEQUENCE_LEN
        maxnamelength = len(max(*namelist, key=len)) + 4
        # make name the same length
        outputname = [i + " " * (maxnamelength - len(i)) for i in namelist]
        seqcount = len(outputname)
        for i in range(maxseqlength):   # cut sequence for MAX_SEQUENCE_LEN
            seqpart = seqarray[:, i * MAX_SEQUENCE_LEN:(i + 1) * MAX_SEQUENCE_LEN]
            for j in range(seqcount):
                IOtext.write(outputname[j] + "".join(seqpart[j]) + "\n")
            IOtext.write("\n")
        seqpart = seqarray[:, MAX_SEQUENCE_LEN * maxseqlength:]
        for j in range(seqcount):
            IOtext.write(outputname[j] + "".join(seqpart[j]) + "\n")
        IOtext.write("\n")


def FASTA(IOtext, *args: Tuple["namelist", "seqarray", "MAX_SEQUENCE_LEN"]) -> dict:
    try:
        namelist, seqarray, MAX_SEQUENCE_LEN = args
    except ValueError:
        seqdict = {}
        name = ""
        for line in IOtext:
            if line.startswith(">"):
                name = line
            else:
                seqdict[name] = seqdict.get(name, "") + clean_seq(line)
        return seqdict
    else:
        maxseqlength = len(seqarray[0]) // MAX_SEQUENCE_LEN
        for name, seq in zip(namelist, seqarray):
            IOtext.write(name)
            for i in range(maxseqlength):
                seqpart = seq[i * MAX_SEQUENCE_LEN:(i + 1) * MAX_SEQUENCE_LEN]
                IOtext.write("".join(seqpart) + "\n")
            seqpart = seq[maxseqlength * MAX_SEQUENCE_LEN:]
            IOtext.write("".join(seqpart) + "\n")


def format_input(format_func, path):
    """
        @description:
        @param format_func:
            the type of the file
            如 FASTA 等
        @param path:
            打开的文件 (__path) 或 str 格式
        @return namelist, seqarray
    """
    if isinstance(path, str):
        seqdict = format_func(path)
    else:
        with path() as _in:
            seqdict = format_func(_in)
    namelist = list(seqdict.keys())
    seqarray = numpy.array([list(seqdict[i]) for i in namelist])
    return namelist, seqarray


def format_output(format_func, path, DEGAP_THRESHOLD, namelist, seqarray, MAX_SEQUENCE_LEN, *args):
    with path('degap_' + str(int(DEGAP_THRESHOLD * 100)), 'w') as output:
        format_func(output, namelist, seqarray, MAX_SEQUENCE_LEN, *args)


# CLUSTAL 文件格式
DEGAP_THRESHOLD = 0.4   # 阈值
MAX_SEQUENCE_LEN = 60   # 一行最多的序列数
if __name__ == "__main__":
    source_path = "demo/degap"
    source_name = "alignment.txt"
    file1 = __Path(source_name, source_path)
    namelist, inputseqarray = format_input(FASTA, file1)
    report_gap(inputseqarray)
    for DEGAP_THRESHOLD in [0.25, ]:  # 0.4, 0.5, 0.75, 0.85, ]:
        DEGAP_THRESHOLD, outputseqarray = delete_gap(
            inputseqarray, DEGAP_THRESHOLD)
        format_output(FASTA, file1, DEGAP_THRESHOLD, namelist,
                      outputseqarray, MAX_SEQUENCE_LEN)
