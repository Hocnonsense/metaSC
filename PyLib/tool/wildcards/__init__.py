# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-26 11:57:34
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-02 19:07:09
 * @FilePath: /metaSC/PyLib/tool/wildcards/__init__.py
 * @Description:
"""

from . import common, formatter, pattern


def fmter(
    _pattern,
    nstep=1,
    keep_missing=False,
    fill_missing=False,
    missing_value=None,
    quote_all: formatter.Union[None, bool, formatter.string.Formatter] = False,
):
    p = formatter.StrPattern(_pattern).xformat(
        nstep=nstep + (1 if nstep else 0),
        keep_missing=keep_missing,
        fill_missing=fill_missing,
        missing_value=missing_value,
        quote_all=quote_all,
    )
    return p


def epder(
    *_patterns,
    keep_missing=False,
    fill_missing=False,
    fill_default=None,
    combinator=formatter.product,
):
    # check if remove missing is provided
    fps = [
        formatter.StrPattern(_pattern).xexpand(
            keep_missing=keep_missing,
            fill_missing=fill_missing,
            missing_value=fill_default,
            combinator=combinator,
        )
        for _pattern in _patterns
    ]

    def _epder(**kwargs):
        return [p for fp in fps for p in fp(**kwargs)]

    return _epder


def glob_path(
    _pattern,
    paths: pattern.Optional[pattern.Iterable] = None,
    followlinks=False,
    restriction: pattern.Optional[dict[str, str]] = None,
    omit_value: pattern.Optional[str] = None,
):
    return pattern.PathPattern(_pattern).glob_path(  # type: ignore
        paths=paths,
        followlinks=followlinks,
        restriction=restriction,
        omit_value=omit_value,
    )


def strip_constraints(_pattern):
    return pattern.PathPattern(_pattern).strip_constraints
