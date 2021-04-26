# -*- coding: utf-8 -*-
"""
    @Date: 2020-02-04 16:21:29
    @LastEditors: Hwrn
    @LastEditTime: 2020-08-18 16:58:24
    @FilePath: /mylib/mylib/biotool/vlookups.py
    @Description:
        查找, 对比两个表
        https://blog.csdn.net/weixin_44881875/article/details/102779975
"""

import pandas
import matplotlib.pyplot


if __name__ == "__main__":
    lines1 = pandas.read_csv('chr.COG_annotations.csv')
    lines2 = pandas.read_csv('D50.vs.D01.fil.csv')

    ndf1 = pandas.DataFrame(
        {'Gene_ID': lines1['Gene_ID'], 'COG_categories': lines1['COG_categories'], })

    df21 = lines2['Gene_ID'].dropna()
    good_n2 = ndf1[ndf1.Gene_ID.isin(df21)]  # Good!

    gn1 = ndf1.groupby('COG_categories')
    good_g1 = gn1.agg(lambda x: list(x)).reset_index()  # Good!
    good_g1['Gene_count_1'] = good_g1.apply(
        lambda x: len(x['Gene_ID']), axis=1)
    good_g1 = good_g1.rename(columns={'Gene_ID': 'Gene_ID_1'})

    good_g2 = good_n2.groupby('COG_categories').agg(
        lambda x: list(x)).reset_index()
    good_g2['Gene_count_2'] = good_g2.apply(
        lambda x: len(x['Gene_ID']), axis=1)

    good_gn = good_g1.merge(good_g2, on=['COG_categories'], how='left')
    print(good_gn)
    # gn1.apply(lambda x:x).values.tolist()
    good_gn.to_csv("result1.csv")

    goodgn = good_gn.set_index('COG_categories')
    goodgn.plot(kind='bar')
    matplotlib.pyplot.show()
