# -*- coding: utf-8 -*-
"""
 * @Date: 2022-04-28 23:53:07
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 00:20:35
 * @FilePath: /metaSC/PyLib/test/tool/test_timer.py
 * @Description:
"""

from PyLib.tool.timer import Timer, TimeCode
from PyLib.PyLibTool.file_info import basicConfig


def test_Timer():
    t = Timer()
    t.make_checkpoint("test")
    timedelta = t.timedelta_to_checkpoint(t.timestamp(), checkpoint_key="test")
    print(
        t.format_time(
            timedelta,
            fmt="{days} days, {hours} hours, {seconds} seconds",
            zero_padding=0,
        )
    )
    print(t.time_elapsed())
    t.format_time(t.timedelta_to_checkpoint(t.timestamp()))
    t = Timer(3)  # 3 checkpoints expected until completion
    for _ in range(3):
        t.make_checkpoint()
        print("complete: %s" % t.complete)
        print(t.eta(fmt="ETA: {seconds} seconds"))


def test_format_time():
    """Run this and visually inspect its working"""
    from datetime import timedelta

    t = Timer()

    for exponent in range(1, 7):
        seconds = 10**exponent
        td = timedelta(seconds=seconds)

        print(f"TESTING {td}")
        fmts = [
            None,
            "SECONDS {seconds}",
            "MINUTES {minutes}",
            "HOURS   {hours}",
            "DAYS    {days}",
            "WEEKS   {weeks}",
            "MINUTES {minutes} SECONDS {seconds}",
            "SECONDS {seconds} MINUTES {minutes}",
            "HOURS   {hours}   MINUTES {minutes}",
            "DAYS    {days}    HOURS   {hours}",
            "WEEKS   {weeks}   DAYS    {days}",
            "WEEKS   {weeks}   HOURS   {hours}",
            "WEEKS   {weeks}   MINUTES {minutes}",
            "DAYS    {days}    MINUTES {minutes}",
            "HOURS   {hours}   SECONDS {seconds}",
            "DAYS    {days}    MINUTES {minutes} SECONDS {seconds}",
            "WEEKS   {weeks}   HOURS {hours}     DAYS    {days}    SECONDS {seconds} MINUTES {minutes}",
        ]
        for fmt in fmts:
            print(str(fmt), t.format_time(td, fmt=fmt))


def test_TimeCode():
    basicConfig(logger_level="DEBUG")
    import time

    # EXAMPLE 1
    with TimeCode() as t:
        time.sleep(5)

    # EXAMPLE 2
    try:
        with TimeCode() as t:
            time.sleep(5)
            raise Exception  # undefined variable
    except Exception:
        pass

    # EXAMPLE 3
    with TimeCode(quiet=True) as t:
        time.sleep(5)
    print(t.time)
