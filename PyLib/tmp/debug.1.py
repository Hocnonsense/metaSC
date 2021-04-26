# -*- coding: utf-8 -*-
"""
    @Date: 2020-01-20 23:12:59
    @LastEditors: Hwrn
    @LastEditTime: 2020-09-14 09:04:46
    @FilePath: /mylib/debug.py
    @Description:
        plot by two depths
"""


from sys import stdout

import matplotlib.pyplot as plt
import numpy


file1, file2 = r"D:\Files\TODU\2020-09-MgAffact\M001_04_map_bbmap.depth", r"D:\Files\TODU\2020-09-MgAffact\M001II_04_map_bbmap.depth"

contig_name1 = {}
count1 = 0
with open(file1) as in1:
    title = in1.readline()
    line = in1.readline()
    while line:
        value = line.split("\t")
        ID, Avg_fold1 = value[0], float(value[1])
        if Avg_fold1:
            contig_name1[ID] = Avg_fold1  # #ID	Avg_fold	Length
        count1 += 1
        line = in1.readline()

contig_pair = {}
count2 = 0
with open(file2) as in2:
    title = in2.readline()
    line = in2.readline()
    while line:
        value = line.split("\t")
        ID, Avg_fold2 = value[0], float(value[1])
        if Avg_fold2 and ID in contig_name1:
            contig_pair[ID] = (contig_name1.pop(ID), Avg_fold2)
        count2 += 1
        line = in2.readline()

# for ID, Avg_fold1 in contig_name1.items():
#     contig_pair[ID] = (Avg_fold1, 0)

# Now, contig_pair is created.
# Then, to ... print it?

count3 = len(contig_pair)
report_table = numpy.array([[a1, a2]
                            for a1, a2 in contig_pair.values()]).T

print(count1, count2, count3)

plt.scatter(report_table[0], report_table[1])
plt.show()


stdout.write("done!")
