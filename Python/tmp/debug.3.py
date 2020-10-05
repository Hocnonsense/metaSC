# -*- coding: utf-8 -*-
"""
    @Date: 2020-01-20 23:12:59
    @LastEditors: Hwrn
    @LastEditTime: 2020-09-15 18:19:03
    @FilePath: /mylib/debug.py
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

colors = ["#15A848", "#15A892", "#1574A8", "#152BA8", "#9215A8",
          "#A81574", "#A8152B", "#A84815", "#A89215", "#39DE1B",
          "#C11BDE", "#D24EE9",
          'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
          'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']
marker = "x"
samples = list(Avg_folds.keys())
for i in range(9):
    print(i, samples[i], colors[i], end=" ")
    part = np.array([Avg_folds[samples[i]], NODE_length])
    part = part[:, part[0] > 0]
    print(len(part[0]))
    plt.scatter(part[0], part[1], s=1, c=colors[i], marker=marker)

plt.xlabel("depth")
plt.ylabel("contig length")
# plt.axis([0, 100000, 500, 1.0e6])
plt.xscale("log")
plt.yscale("log")
plt.legend(samples, fontsize=7)
plt.savefig('depths-length.png', figsize=(40, 30),
            dpi=500, bbox_inches='tight')

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


stdout.write("all right~")
