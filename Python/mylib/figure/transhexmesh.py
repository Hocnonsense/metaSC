# -*- coding: utf-8 -*-
"""
 * @Date: 2020-01-20 23:12:59
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-05 15:34:31
 * @FilePath: /HScripts/Python/mylib/figure/transhexmesh.py
 * @Description: 画一个透明的六边形网格, 使之能够叠加在其他图案上.
    To draw a six-edge-axis
"""
import math
import PIL.Image as Image
import numpy as np
from matplotlib import pyplot as plt


def draw_transparent_hexagonal_mesh(row, column, unit_length):
    """
    画
    """
    width = 7
    slope = math.sqrt(3)/2
    color = [0, 0, 0, 255]
    row, col, collen = row, column, unit_length
    rowlen = collen*slope  # the width of the mesh
    imgarray = np.zeros(shape=(int(col * 3 / 2 * collen +
                                   collen / 2), int(row * rowlen * 2 + rowlen) + width, 4))
    for j in range(col // 2):
        _j = j * 3 * collen
        dj = int(_j)
        Aj = int(_j + collen / 2)
        Bj = int(_j + collen * 3 / 2)
        Cj = int(_j + collen * 2)
        Dj = int(_j + collen * 3)
        for i in range(row):
            # pylint: disable = anomalous-backslash-in-string
            """   // d \\
                A         a
               ||         |
                B         b
                  \\ C //   \\ c
                    ||        |
                     D        d
            """
            _i = i*2*rowlen
            ABi = int(_i)
            CDi = int(_i + rowlen)
            imgarray[Aj:Bj, ABi: ABi + width, :] = color  # A|B
            imgarray[Cj:Dj, CDi: CDi + width, :] = color  # C|D
            for x in range(dj, Aj):
                y1 = int(CDi - (x - dj) / (Aj - dj) * (CDi - ABi))
                y2 = int(CDi + (x - dj) / (Aj - dj) * (CDi - ABi))
                imgarray[x, y1: y1 + width, :] = color  # A/d
                imgarray[x, y2: y2 + width, :] = color  # d\a
            for x in range(Bj, Cj):
                y3 = int(ABi + (x - Bj) / (Cj - Bj) * (CDi - ABi))
                y4 = int(CDi + (Cj - x) / (Cj - Bj) * (CDi - ABi))
                imgarray[x, y3: y3 + width, :] = color  # B\C
                imgarray[x, y4: y4 + width, :] = color  # C/b

        imgarray[Aj:Bj, int(row * 2 * rowlen): int(row *
                                                   2 * rowlen) + width, :] = color  # last a|b
        imgarray[Cj:Dj, int(row * 2 * rowlen + rowlen
                            ): int(row * 2 * rowlen + rowlen) + width, :] = color  # last c|d
        for x in range(Dj, int(_j+collen*7/2)):
            y1 = int(row*2*rowlen+rowlen - (x-Dj)/(Aj-dj)*(CDi-ABi))
            imgarray[x, y1:y1+width, :] = color  # last A/d
        for x in range(Bj, Cj):
            y3 = int(row*2*rowlen + (x-Bj)/(Cj-Bj)*(CDi-ABi))
            imgarray[x, y3:y3+width, :] = color  # last b\c
    if col % 2 == 0:
        for i in range(row):
            _j = col // 2 * 3 * collen
            dj = int(_j)
            Aj = int(_j + collen / 2)
            Bj = int(_j+collen*3/2)
            Cj = int(_j+collen*2)
            Dj = int(_j+collen*3)
            _i = i*2*rowlen
            ABi = int(_i)
            CDi = int(_i+rowlen)
            for x in range(dj, Aj):
                y1 = int(CDi+rowlen*2 - (x-dj)/(Aj-dj)*(CDi-ABi))
                y2 = int(CDi + (x-dj)/(Aj-dj)*(CDi-ABi))
                imgarray[x, y1:y1+width, :] = color  # last
                imgarray[x, y2:y2+width, :] = color
    else:
        _j = col//2*3*collen
        dj = int(_j)
        Aj = int(_j + collen / 2)
        Bj = int(_j + collen * 3 / 2)
        Cj = int(_j + collen * 2)
        for i in range(row):
            _i = i * 2 * rowlen
            ABi = int(_i)
            CDi = int(_i + rowlen)
            imgarray[Aj:Bj, ABi: ABi + width, :] = color
            for x in range(dj, Aj):
                y1 = int(CDi - (x - dj) / (Aj - dj) * (CDi - ABi))
                y2 = int(CDi + (x - dj) / (Aj - dj) * (CDi - ABi))
                imgarray[x, y1: y1 + width, :] = color
                imgarray[x, y2: y2 + width, :] = color
            imgarray[Cj:Dj, CDi: CDi + width, :] = color
            for x in range(Bj, Cj):
                y3 = int(ABi + (x - Bj) / (Cj - Bj) * (CDi - ABi))
                y4 = int(CDi + (Cj - x) / (Cj - Bj) * (CDi - ABi))
                imgarray[x, y3: y3 + width, :] = color
                imgarray[x, y4: y4 + width, :] = color
        imgarray[Aj:Bj, int(row * 2 * rowlen): int(row *
                                                   2 * rowlen) + width, :] = color
    img = Image.fromarray(np.uint8(imgarray))
    plt.imshow(img)
    plt.show()
    img.save("a.png")


if __name__ == "__main__":
    draw_transparent_hexagonal_mesh(20, 20, 50)
