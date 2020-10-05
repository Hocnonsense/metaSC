# -*- coding: utf-8 -*-
"""
 * @Date: 2020-02-05 17:15:25
 * @LastEditors: Hwrn
 * @LastEditTime: 2020-10-05 15:17:43
 * @FilePath: /HScripts/Python/mylib/tool/timer.py
 * @Description: Record time.
"""


from datetime import datetime, timedelta
from typing import Optional


def get_time(str_time: str()):
    """Get time"""
    return datetime.strptime(str_time, "%Y-%m-%d %H:%M:%S")


def timestamp(dt: Optional[datetime] = None) -> str:
    """ @description 返回时间戳
            Construct a datetime-timestamp from `dt` [default = `now()`], such as
            we use to timestamp a simulation output directory.
    """
    if not dt:
        dt = datetime.now()

    return dt.strftime('%Y%m%d.%H%M%S')


class Timer(object):
    """
    member:
        __init_at: datetime.datetime
            the time set the timer, useless
        begin_at: datetime.datetime | None
            if before now: timing on, since then
            elif after now: timing off, start since then
            elif is None: timing off, start when called
        alarm_delta: datetime.timedelta | None
            if is None: stop when called
            else: stop if time are all spend
    property:
        now: datetime.datetime.now()
        end_at: datetime.datetime | None

    method:
        is_time_on():
            if timer is counting: True
            else: False
        time_on(alarm_time, begin_at):
            if timer is on: return False
            else:
                change alarm_time to <class datetime.timedelta>
                set alarm_delta and begin_at
                return True
        time_off():
            if timer is off:
                return False
            else:
                close timer
                if alarm_delta last some time: update alarm_delta
                return begin_at, passed time, remain time
        check():
            if timer is off: return False
            else: return delta to end

    """
    def __init__(self):
        self.__init_at = datetime.now()
        self.begin_at: datetime = None
        self.alarm_delta: timedelta = None
        self.last_time_line = None, None, None

    def is_time_on(self):
        """
        Check timer state.
        If timer is counting: True
        else: False
        """
        if self.alarm_delta is None and self.begin_at is None:
            return False
        else:
            return True

    def time_on(self, alarm_time, begin_at: datetime = "now"):
        """
        Make timer start at time `begin_at` and stay awake durning `alarm_time`.
        @param alarm_time if alarm_time is datetime, it will end at that time.
                          elif alarm time is timedelta, this timer will end at that time.
        @param begin_at time it will start at

        If timer is on: return False
        else:
            change alarm_time to <class datetime.timedelta>
            set alarm_delta and begin_at
            return True
        """
        if self.is_time_on():
            return False
        if begin_at == "now":
            begin_at = datetime.now()
        if isinstance(alarm_time, datetime):
            alarm_time = alarm_time - begin_at
        elif not isinstance(alarm_time, timedelta):
            alarm_time = timedelta(seconds=alarm_time)
        self.begin_at = begin_at
        self.alarm_delta = alarm_time
        return self.check()

    def time_off(self):
        """
        Turn off timer if it is timing.
        if timer is off:
            return False
        else:
            close timer
            if alarm_delta last some time: update alarm_delta
            return begin_at, passed time, remain time
        """
        if not self.is_time_on():
            return False
        now = datetime.now()
        passed_delta = now - self.begin_at
        if self.alarm_delta:
            remain_delta = self.alarm_delta - passed_delta
            time_line = self.begin_at, passed_delta, remain_delta
            if timedelta(0) < remain_delta:
                self.alarm_delta = remain_delta
        else:
            time_line = self.begin_at, passed_delta, None
        self.begin_at = None
        self.alarm_delta = None
        self.last_time_line = time_line
        return time_line

    def check(self):
        """
        Look at the timer, report time.
        if timer is off: return False
        else: return delta to end
        """
        if not self.is_time_on():
            return False
        else:
            now = datetime.now()
            passed_delta = now - self.begin_at
            remain_delta = (self.alarm_delta or timedelta()) - passed_delta
            if self.alarm_delta and self.alarm_delta < passed_delta:
                self.time_off()
            return remain_delta

    @property
    def now(self):
        """Return datetime.now()."""
        return datetime.now()

    @property
    def end_at(self):
        """Return when the timer will be off."""
        if self.begin_at and self.alarm_delta:
            return self.begin_at + self.alarm_delta
        else:
            return None

    def regular_time(self, seconds: int = "self.alarm_delta"):
        """Return time in a regurlar format."""
        if seconds == "self.alarm_delta":
            seconds = self.alarm_delta.seconds
        day, hour, minute, sec = 0, 0, 0, int(seconds)
        minute, sec = sec // 60, sec % 60
        hour, minute = minute // 60, minute % 60
        day, hour = hour // 24, hour % 24
        return day, hour, minute, sec

    def show_time(self, seconds: int = "self.alarm_delta"):
        """Show time in a regurlar format."""
        if isinstance(seconds, datetime):
            return ":".join([str(i) for i in (seconds.hour, seconds.minute, seconds.second)])
        elif seconds:
            day, hour, minute, sec = self.regular_time(seconds)
            if day:
                return ":".join([str(i) for i in [day, hour, minute, sec]])
            if hour:
                return ":".join([str(i) for i in [hour, minute, sec]])
            return ":".join([str(i) for i in [minute, sec]])
        else:
            return ""


if __name__ == "__main__":
    timer_1 = Timer()
    timer_1.time_on(3)
    while timer_1.is_time_on():
        print("\r" + str(timer_1.check()), end="")
    print()
    print(timer_1.last_time_line)
    timer_1.time_on(4)
