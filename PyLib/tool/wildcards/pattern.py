# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-26 11:03:36
 * @author: Johannes Köster
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-02 16:12:07
 * @FilePath: /metaSC/PyLib/tool/wildcards/pattern.py
 * @Description: extract regex methods from snakemake
"""


import collections
from collections.abc import Iterable
import os
import re
from itertools import chain
from typing import Iterable, Optional, overload

from .common import PathLike, Pattern, Unformattable

_wildcard_regex = re.compile(
    r"""
    \{
        (?=(   # This lookahead assertion emulates an 'atomic group'
               # which is required for performance
            \s*(?P<name>\w+)                    # wildcard name
            (\s*,\s*
                (?P<constraint>                 # an optional constraint
                    ([^{}]+ | \{\d+(,\d+)?\})*  # allow curly braces to nest one level
                )                               # ...  as in '{w,a{3,5}}'
            )?\s*
        ))\1
    \}
    """,
    re.VERBOSE,
)


class PathUnformattable(Unformattable):
    pass


class PathPattern(Pattern):
    def __init__(self, pattern, _type=None, default_constrains=".+"):
        super().__init__(pattern, _type)
        if _type is None:
            self.check_pattern()
        self.default_constrains = default_constrains

    def __repr__(self) -> str:
        return self.pattern

    def check_pattern(self):
        try:
            self.regex
        except re.error as e:
            e_msg = e.msg
        else:
            e_msg = ""
        finally:
            if e_msg:
                raise PathUnformattable(
                    "Number is not a valid value name when formatting\n"
                    + e_msg
                    + f"\nPlease check: {self.pattern}"
                )
        if "{}" in self.pattern:
            raise PathUnformattable(
                "Unnamed bracket is not allowed when formatting"
                + f"\nPlease check: {self.pattern}"
            )

    @property
    def names(self):
        names_ = [
            match.group("name") for match in _wildcard_regex.finditer(self.pattern)
        ]
        return sorted(set(names_), key=names_.index)

    @property
    def constraints(self):
        constraints: dict[str, re.Pattern] = {}
        for match in _wildcard_regex.finditer(self.pattern):
            if match.group("constraint"):
                constraints[match.group("name")] = re.compile(match.group("constraint"))
        return constraints

    def replace_constraint(
        self,
        constraints: dict[str, Optional[str]],
    ):
        """Update wildcard constraints

        Args:
        constraints (dict): dictionary of wildcard:constraint key-value pairs
        """

        def replace_constraint(match: re.Match):
            name = match.group("name")
            constraint = match.group("constraint")
            newconstraint = constraints.get(name)
            if name in examined_names:
                return match.group(0)
            examined_names.add(name)
            # Don't override if constraint already set
            if constraint is not None:
                return match.group(0)
            # Only update if a new constraint has actually been set
            elif newconstraint is not None:
                return f"{{{name},{newconstraint}}}"
            else:
                return match.group(0)

        examined_names: set[str] = set()
        updated = _wildcard_regex.sub(replace_constraint, self.pattern)

        # TODO: inherit flags as in snakemake
        # if isinstance(pattern, AnnotatedString):
        #    updated = AnnotatedString(updated)
        #    updated.flags = dict(pattern.flags)
        return type(self)(updated, self._type)

    @property
    def strip_constraints(self):
        """Return a string that does not contain any wildcard constraints."""
        # if is_callable(pattern):
        #    # do not apply on e.g. input functions
        #    return pattern

        def strip_constraint(match: re.Match):
            return "{{{}}}".format(match.group("name"))

        return type(self)(
            _wildcard_regex.sub(strip_constraint, self.pattern), self._type
        )

    @property
    def regex(self):
        f = []
        last = 0
        wildcards = set()
        for match in _wildcard_regex.finditer(self.pattern):
            f.append(re.escape(self.pattern[last : match.start()]))
            wildcard = match.group("name")
            if wildcard in wildcards:
                if match.group("constraint"):
                    raise ValueError(
                        "Constraint regex must be defined only in the first "
                        "occurence of the wildcard in a string."
                    )
                f.append(f"(?P={wildcard})")
            else:
                wildcards.add(wildcard)
                f.append(
                    "(?P<{}>{})".format(
                        wildcard,
                        match.group("constraint")
                        if match.group("constraint")
                        else ".+",
                    )
                )
            last = match.end()
        f.append(re.escape(self.pattern[last:]))
        f.append("$")  # ensure that the match spans the whole file
        return re.compile("".join(f))

    def glob(self, paths: Iterable[PathLike]):
        """
        Glob the values of the wildcards by matching the given pattern to the filesystem.
        Returns a named tuple with a list of values for each wildcard.
        """
        # an ordered dict can be more accurate
        if not self.names:
            return []
        Wildcards = collections.namedtuple("Wildcards", self.names)  # type: ignore

        for f in paths:
            match = re.match(self.regex, str(f))
            if match:
                yield PathLike, (Wildcards(**match.groupdict()))

    @overload
    def glob_path(self, paths: Iterable):
        pass

    @overload
    def glob_path(self, *, followlinks=False):
        pass

    @overload
    def glob_path(self, *, restriction=None, omit_value=None):
        pass

    def glob_path(
        self,
        paths: Optional[Iterable] = None,
        followlinks=False,
        restriction: Optional[dict[str, str]] = None,
        omit_value: Optional[str] = None,
    ):
        pattern_ = type(self)(os.path.normpath(self.pattern), str)
        first_wildcard = re.search("{[^{]", pattern_.pattern)
        if paths is None:
            dirname = (
                os.path.dirname(pattern_.pattern[: first_wildcard.start()])
                if first_wildcard
                else os.path.dirname(pattern_.pattern)
            ) or "."
            paths = (
                os.path.normpath(os.path.join(dirpath, f)) if dirpath != "." else f
                for dirpath, dirnames, filenames in os.walk(
                    dirname, followlinks=followlinks
                )
                for f in chain(filenames, dirnames)
            )
        for f, wildcards in pattern_.glob((os.path.normpath(p) for p in paths)):
            if restriction is not None:
                invalid = any(
                    omit_value not in v and v != wildcards[k]  # type: ignore
                    for k, v in restriction.items()
                )
                if not invalid:
                    yield f, wildcards
            else:
                yield f, wildcards
