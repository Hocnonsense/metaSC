# -*- coding: utf-8 -*-
"""
    @Date: 2020-01-20 23:12:59
    @LastEditors: Hwrn
    @LastEditTime: 2020-09-15 13:03:03
    @FilePath: /mylib/debug.2.py
    @Description:
"""


import pickle
from sys import stdout
import numpy as np
import matplotlib.pyplot as plt


with open(r"D:\Files\TODU\2020-09-MgAffact\Avg_folds.pickle", "rb") as fin:
    Avg_folds = pickle.load(fin)

# Draw pair picture for each data
# First, draw depth-length-[color=file]
# get length from .depth

# make ref dict for contigs
# now, we can only use the depth._cov file
NODE_length = np.zeros(len(list(Avg_folds.values())[0]), dtype=int)
ref_name = r"D:\Files\TODU\2020-09-MgAffact\origin-depth._cov"
# origin-02-length_cut.fa like ">{contig_name}\n{sequence}\n"
with open(ref_name) as ref:
    title = ref.readline()  # title
    line = ref.readline()  # loop starts
    while line:
        ID = line.strip().split()[0]
        # e.g. NODE_1_length_3290939_cov_10.181947
        msg = ID.split("_")
        NODE, length = int(msg[1]) - 1, int(msg[3])
        NODE_length[NODE] = length
        # LOOP
        line = ref.readline()

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
marker = "x"
samples = list(Avg_folds.keys())
# for i in range(9):
#     print(i, samples[i], colors[i], end=" ")
#     part = np.array([Avg_folds[samples[i]], NODE_length])
#     part = part[:, part[0] > 0]
#     print(len(part[0]))
#     plt.scatter(part[0], part[1], s=1, c=colors[i], marker=marker)
#
# plt.xlabel("depth")
# plt.ylabel("contig length")
# plt.axis([0, 100000, 500, 1.0e6])
# plt.legend(samples)
#
# plt.savefig('depths-length1.png', figsize=(40, 30),
#             dpi=500, bbox_inches='tight')

# (mylib) D:\Code\python\mylib> cd d:\Code\python\mylib && cmd /C "D:\Code\python\.Anaconda\envs\mylib\python.exe c:\Users\Hwrn\.vscode\extensions\ms-python.python-2020.8.108011\pythonFiles\lib\python\debugpy\launcher 64451 -- d:\Code\python\mylib\debug.py "
# 0 E001_04_map_bbmap.depth tab:blue 6169
# 1 E123_04_map_bbmap.depth tab:orange 50907
# 2 E152_04_map_bbmap.depth tab:green 977037
# 3 M001II_04_map_bbmap.depth tab:red 5228
# 4 M001_04_map_bbmap.depth tab:purple 1370062
# 5 M086II_04_map_bbmap.depth tab:brown 53463
# 6 M086_04_map_bbmap.depth tab:pink 1388588
# 7 M123_04_map_bbmap.depth tab:gray 81007
# 8 M152_04_map_bbmap.depth tab:olive 54913
# all right~

print(len(NODE_length), "scaffolds:", np.min(NODE_length), np.max(NODE_length))
clengths = np.log(NODE_length-499)  # minlength = 500
cm = plt.cm.get_cmap('jet')
plt.scatter(NODE_length, clengths, s=1, c=clengths, cmap=cm)
plt.xlabel("contig length")
plt.savefig("length-color.png")

for i in range(9):
    for j in range(i, 9):
        namei = samples[i].split("_")[0]
        namej = samples[j].split("_")[0]
        cross = np.array(
            [Avg_folds[samples[i]], Avg_folds[samples[j]], clengths])
        cross = cross[:, cross[0] > 0]
        cross = cross[:, cross[1] > 0]
        print(namei, namej, len(cross[0]))
        plt.clf()
        plt.scatter(cross[0], cross[1], s=1, c=cross[2], marker="x", cmap=cm)
        plt.xlabel(namei)
        plt.ylabel(namej)
        plt.title("cross map depth with length")
        plt.axis([0, 4000, 0, 4000])
        plt.savefig('{}-{}.png'.format(namei, namej),
                    figsize=(40, 40), dpi=500, bbox_inches='tight')

# (mylib) D:\Code\python\mylib> cd d:\Code\python\mylib && cmd /C "D:\Code\python\.Anaconda\envs\mylib\python.exe c:\Users\Hwrn\.vscode\extensions\ms-python.python-2020.8.108011\pythonFiles\lib\python\debugpy\launcher 62284 -- d:\Code\python\mylib\debug.2.py "
# 1589718 scaffolds: 500 3290939
# E001 E001 6169
# E001 E123 4688
# E001 E152 2912
# E001 M001II 3499
# E001 M001 1417
# E001 M086II 4356
# E001 M086 816
# E001 M123 4394
# E001 M152 1563
# E123 E123 50907
# E123 E152 35751
# E123 M001II 2684
# E123 M001 1907
# E123 M086II 5936
# E123 M086 1219
# E123 M123 23675
# E123 M152 11874
# E152 E152 977037
# E152 M001II 985
# E152 M001 904242
# E152 M086II 5585
# E152 M086 912248
# E152 M123 23757
# E152 M152 13381
# M001II M001II 5228
# M001II M001 1192
# M001II M086II 5000
# M001II M086 564
# M001II M123 2840
# M001II M152 824
# M001 M001 1370062
# M001 M086II 3917
# M001 M086 1363260
# M001 M123 2409
# M001 M152 4711
# M086II M086II 53463
# M086II M086 3585
# M086II M123 8314
# M086II M152 3478
# M086 M086 1388588
# M086 M123 2341
# M086 M152 4203
# M123 M123 81007
# M123 M152 15810
# M152 M152 54913
# all right~


stdout.write("all right~")
