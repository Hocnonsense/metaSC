# -*- coding: utf-8 -*-
"""
 * @Date: 2020-02-05 17:15:25
 * @LastEditors: Hwrn
 * @LastEditTime: 2022-04-29 00:18:44
 * @FilePath: /metaSC/PyLib/tool/timer.py
 * @Description: Record time.
    update from anvio.terminal
"""

from datetime import datetime, timedelta
from collections import OrderedDict
import time
from PyLib.PyLibTool.file_info import verbose_import
from PyLib.tool.ttycolors import color_text
import pandas as pd

logger = verbose_import(__name__, __doc__)


class Timer:
    """Manages an ordered dictionary, where each key is a checkpoint name and value is a timestamp."""

    def __init__(
        self, required_completion_score=None, initial_checkpoint_key=0, score=0
    ):
        self.timer_start = self.timestamp()
        self.initial_checkpoint_key = initial_checkpoint_key
        self.last_checkpoint_key = self.initial_checkpoint_key
        self.checkpoints = OrderedDict([(initial_checkpoint_key, self.timer_start)])
        self.num_checkpoints = 0

        self.required_completion_score = required_completion_score
        self.score = score
        self.complete = False

        self.last_eta = None
        self.last_eta_timestamp = self.timer_start

        self.scores = {self.initial_checkpoint_key: self.score}

    def timestamp(self):
        return datetime.fromtimestamp(time.time())

    def timedelta_to_checkpoint(self, timestamp, checkpoint_key=None):
        timedelta = (
            timestamp - self.checkpoints[checkpoint_key or self.initial_checkpoint_key]
        )
        return timedelta

    def make_checkpoint(self, checkpoint_key=None, increment_to=None):
        if not checkpoint_key:
            checkpoint_key = self.num_checkpoints + 1

        if checkpoint_key in self.checkpoints:
            raise Exception(
                f"Timer.make_checkpoint :: {checkpoint_key} already exists as a checkpoint key. "
                f"All keys must be unique"
            )

        checkpoint = self.timestamp()

        self.checkpoints[checkpoint_key] = checkpoint
        self.last_checkpoint_key = checkpoint_key

        self.num_checkpoints += 1

        if increment_to:
            self.score = increment_to
        else:
            self.score += 1

        self.scores[checkpoint_key] = self.score

        if (
            self.required_completion_score
            and self.score >= self.required_completion_score
        ):
            self.complete = True

        return checkpoint

    def gen_report(self, title="Time Report"):
        checkpoint_last = self.initial_checkpoint_key

        logger.warning("", header=title, lc="yellow", nl_before=1, nl_after=0)

        for checkpoint_key, checkpoint in self.checkpoints.items():
            if checkpoint_key == self.initial_checkpoint_key:
                continue

            logger.info(
                str(checkpoint_key),
                "+%s"
                % self.timedelta_to_checkpoint(
                    checkpoint, checkpoint_key=checkpoint_last
                ),
            )
            checkpoint_last = checkpoint_key

        logger.info(
            f"Total elapsed="
            f"{self.timedelta_to_checkpoint(checkpoint, checkpoint_key=self.initial_checkpoint_key)}"
        )

    def gen_dataframe_report(self):
        """Returns a dataframe"""

        d = {"key": [], "time": [], "score": []}
        for checkpoint_key, checkpoint in self.checkpoints.items():
            d["key"].append(checkpoint_key)
            d["time"].append(checkpoint)
            d["score"].append(self.scores[checkpoint_key])

        return pd.DataFrame(d)

    def gen_file_report(self, filepath):
        """Writes to filepath, will overwrite"""

        self.gen_dataframe_report().to_csv(filepath, sep="\t", index=False)

    def calculate_time_remaining(self, infinite_default="∞:∞:∞"):
        if self.complete:
            return timedelta(seconds=0)
        if not self.required_completion_score:
            return None
        if not self.score:
            return infinite_default

        time_elapsed = self.checkpoints[self.last_checkpoint_key] - self.checkpoints[0]
        fraction_completed = self.score / self.required_completion_score
        time_remaining_estimate = time_elapsed / fraction_completed - time_elapsed

        return time_remaining_estimate

    def eta(self, fmt=None, zero_padding=0):
        # Calling format_time hundreds or thousands of times per second is expensive. Therefore if
        # eta was called within the last half second, the previous ETA is returned without further
        # calculation.
        eta_timestamp = self.timestamp()
        if (
            eta_timestamp - self.last_eta_timestamp < timedelta(seconds=0.5)
            and self.num_checkpoints > 0
        ):
            return self.last_eta

        eta = self.calculate_time_remaining()
        eta = (
            self.format_time(eta, fmt, zero_padding)
            if isinstance(eta, timedelta)
            else str(eta)
        )

        self.last_eta = eta
        self.last_eta_timestamp = eta_timestamp

        return eta

    def time_elapsed(self, fmt=None):
        return self.format_time(
            self.timedelta_to_checkpoint(self.timestamp(), checkpoint_key=0), fmt=fmt
        )

    def format_time(self, timedelta, fmt="{hours}:{minutes}:{seconds}", zero_padding=2):
        """Formats time

        Examples of `fmt`. Suppose the timedelta is seconds = 1, minutes = 1, hours = 1.

            {hours}h {minutes}m {seconds}s  --> 01h 01m 01s
            {seconds} seconds               --> 3661 seconds
            {weeks} weeks {minutes} minutes --> 0 weeks 61 minutes
            {hours}h {seconds}s             --> 1h 61s
        """

        unit_hierarchy = ["seconds", "minutes", "hours", "days", "weeks"]
        unit_denominations = {
            "weeks": 7,
            "days": 24,
            "hours": 60,
            "minutes": 60,
            "seconds": 1,
        }

        if not fmt:
            # use the highest two non-zero units, e.g. if it is 7200s, use {hours}h{minutes}m
            seconds = int(timedelta.total_seconds())
            if seconds < 60:
                fmt = "{seconds}s"
            else:
                m = 1
                for i, unit in enumerate(unit_hierarchy):
                    if not seconds // (m * unit_denominations[unit]) >= 1:
                        fmt = "{%s}%s{%s}%s" % (
                            unit_hierarchy[i - 1],
                            unit_hierarchy[i - 1][0],
                            unit_hierarchy[i - 2],
                            unit_hierarchy[i - 2][0],
                        )
                        break
                    elif unit == unit_hierarchy[-1]:
                        fmt = "{%s}%s{%s}%s" % (
                            unit_hierarchy[i],
                            unit_hierarchy[i][0],
                            unit_hierarchy[i - 1],
                            unit_hierarchy[i - 1][0],
                        )
                        break
                    else:
                        m *= unit_denominations[unit]

        # parse units present in fmt
        format_order = []
        for i, x in enumerate(fmt):
            if x == "{":
                for j, k in enumerate(fmt[i:]):
                    if k == "}":
                        unit = fmt[i + 1 : i + j]
                        format_order.append(unit)
                        break

        if not format_order:
            raise Exception(
                f"Timer.format_time :: fmt = '{fmt}' contains no time units."
            )

        for unit in format_order:
            if unit not in unit_hierarchy:
                raise Exception(
                    f"Timer.format_time :: '{unit}' is not a valid unit. Use any of {', '.join(unit_hierarchy)}."
                )

        # calculate the value for each unit (e.g. 'seconds', 'days', etc) found in fmt
        format_values_dict = {}
        smallest_unit = unit_hierarchy[
            [unit in format_order for unit in unit_hierarchy].index(True)
        ]
        units_less_than_or_equal_to_smallest_unit = unit_hierarchy[::-1][
            unit_hierarchy[::-1].index(smallest_unit) :
        ]
        seconds_in_base_unit = 1
        for a in [
            v
            for k, v in unit_denominations.items()
            if k in units_less_than_or_equal_to_smallest_unit
        ]:
            seconds_in_base_unit *= a
        r = int(timedelta.total_seconds()) // seconds_in_base_unit

        for i, lower_unit in enumerate(unit_hierarchy):
            if lower_unit in format_order:
                m = 1
                for upper_unit in unit_hierarchy[i + 1 :]:
                    m *= unit_denominations[upper_unit]
                    if upper_unit in format_order:
                        (
                            format_values_dict[upper_unit],
                            format_values_dict[lower_unit],
                        ) = divmod(r, m)
                        break
                else:
                    format_values_dict[lower_unit] = r
                    break
                r = format_values_dict[upper_unit]

        format_values = [format_values_dict[unit] for unit in format_order]

        style_str = "0" + str(zero_padding) if zero_padding else ""
        for unit in format_order:
            fmt = fmt.replace("{%s}" % unit, "%" + "%s" % (style_str) + "d")
        formatted_time = fmt % (*[format_value for format_value in format_values],)

        return formatted_time


class TimeCode(object):
    """Time a block of code.

    This context manager times blocks of code, and calls run.info afterwards to report
    the time (unless quiet = True). See also time_program()

    Parameters
    ==========
    sc: 'green'
        run info color with no runtime error
    success_msg: None
        If None, it is set to 'Code ran succesfully in'
    fc: 'green'
        run info color with runtime error
    failure_msg: None
        If None, it is set to 'Code failed within'
    run: Run()
        Provide a pre-existing Run instance if you want
    quiet: False,
        If True, run.info is not called and datetime object is stored
        as `time` (see examples)
    suppress_first: 0,
        Supress output if code finishes within this many seconds.
    """

    single_line_prefixes = {0: "✓ ", 1: "✖ "}

    def __init__(
        self,
        success_msg=None,
        sc="green",
        fc="red",
        failure_msg=None,
        quiet=False,
        suppress_first=0,
    ):

        self.quiet = quiet
        self.suppress_first = suppress_first
        self.sc, self.fc = sc, fc

        self.s_msg = success_msg or "Code finished after "
        self.f_msg = failure_msg or "Code encountered error after "

    def __enter__(self):
        self.timer = Timer()
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        self.time = self.timer.timedelta_to_checkpoint(self.timer.timestamp())

        if self.quiet or self.time <= timedelta(seconds=self.suppress_first):
            return

        return_code = 0 if exception_type is None else 1

        msg, color = (self.s_msg, self.sc) if not return_code else (self.f_msg, self.fc)
        message_line = color_text(
            f"\n{self.single_line_prefixes[return_code]}{msg}{self.time}\n", color
        )
        logger.info(message_line)
