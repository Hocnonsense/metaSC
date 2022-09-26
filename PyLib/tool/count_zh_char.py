# -*- coding: utf-8 -*-
"""
 * @Date: 2022-09-25 16:01:17
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-09-25 22:04:41
 * @FilePath: /TODO/utils/count_zh_char.py
 * @Description:
    counting chinese characters in given file and given line range
"""

import sys
import re


def count_valid_md_zh_chars(line: str, skiphead=re.compile(r"^#+ ")):
    # ignore header in markdown file
    not_zh_chars = re.compile(r"[\u0000-\u0127（），。、≥≤]")
    if re.match(skiphead, line.strip()):
        return 0
    clean_zh = re.sub(not_zh_chars, "", line)
    return len(clean_zh)


if __name__ == "__main__":
    file, line_range = sys.argv[1:]
    splits = re.sub(r"[^\d]+", ",", line_range).split(",")
    drop_flag = True
    last_index = 0
    total_zh_chars = 0
    with open(file) as fi:
        for i in splits:
            index = int(i) - drop_flag
            for _, line in zip(range(last_index, index), fi):
                if not drop_flag:
                    total_zh_chars += count_valid_md_zh_chars(line)
                    # print(line)
                    # "drop" if drop_flag else "keep",
            drop_flag = not drop_flag
            last_index = index
        if not drop_flag:
            for line in fi:
                total_zh_chars += count_valid_md_zh_chars(line)
    print(total_zh_chars)

