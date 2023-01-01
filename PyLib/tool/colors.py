# -*- coding: utf-8 -*-
"""
 * @Date: 2023-01-01 13:57:08
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-01-01 16:32:06
 * @FilePath: /2022_09-M_mem/home/hwrn/software/metaSC/PyLib/tool/colors.py
 * @Description:
"""
# """

import colorsys
import random


def get_n_random_hls_colors(num):
    # http://t.csdn.cn/ZIEZc
    hls_colors = []
    i = 0
    step = 360.0 / num
    while i < 360:
        h = i
        l = 30 + random.random() * 40  # median 50
        s = 80 + random.random() * 20  # max 100
        _hlsc = [h / 360.0, l / 100.0, s / 100.0]
        hls_colors.append(_hlsc)
        i += step

    return hls_colors


def get_n_hls_colors(num):
    l_range = (30, 70, 22)
    s_range = (60, 100, 22)

    hls_colors = []
    l = l_range[0]
    s = s_range[0]
    h_step = 360.0 / num
    i = 0
    while i < 360:
        h = i
        _hlsc = [h / 360.0, l / 100, s / 100]
        hls_colors.append(_hlsc)
        i += h_step
        if l == l_range[1]:
            l = l_range[0]
        else:
            l += l_range[2]
            if l > l_range[1]:
                l -= l_range[1] - l_range[0]
        if s == s_range[1]:
            s = s_range[0]
        else:
            s += s_range[2]
            if s > s_range[1]:
                s -= s_range[1] - s_range[0]
        if (l == l_range[1]) and (s == s_range[1]):
            l = l_range[0]
            s = s_range[0]

    return hls_colors


def ncolors(num, random=True):
    # http://t.csdn.cn/ZIEZc
    rgb_colors = []
    if num < 1:
        return rgb_colors
    hls_colors = (get_n_random_hls_colors if random else get_n_hls_colors)(num)
    for hlsc in hls_colors:
        _r, _g, _b = colorsys.hls_to_rgb(hlsc[0], hlsc[1], hlsc[2])
        r, g, b = [int(x * 255.0) for x in (_r, _g, _b)]
        rgb_colors.append([r, g, b])

    return rgb_colors


def rgb2str(r: int, g: int, b: int):
    # 原文链接：https://blog.csdn.net/qq_43707919/article/details/109320492
    hex1 = hex(r)[2:].upper()
    hex2 = hex(g)[2:].upper()
    hex3 = hex(b)[2:].upper()
    # 十进制转16进制时会出现缺省零的情况，用rjust函数可在字符串左侧填充0
    # 同理 ljust函数可在字符串的右侧填充0
    hex1 = hex1.rjust(2, "0")
    hex2 = hex2.rjust(2, "0")
    hex3 = hex3.rjust(2, "0")
    outputstr = "#"
    outputstr = outputstr + hex1 + hex2 + hex3
    return outputstr


def str2rgb(rgb: str):
    # 原文链接：https://blog.csdn.net/qq_43707919/article/details/109320492
    r = rgb[1:3]
    g = rgb[3:5]
    b = rgb[5:7]
    num1 = int("0x" + r, 16)
    num2 = int("0x" + g, 16)
    num3 = int("0x" + b, 16)
    return num1, num2, num3


if __name__ == "__main__":
    palletes1 = ncolors(40, random=False)
    [rgb2str(*i) for i in palletes1]
