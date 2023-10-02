# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-25 21:27:35
 * @author: Johannes Köster
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-10-02 18:49:18
 * @FilePath: /metaSC/PyLib/tool/wildcards/formatter.py
 * @Description:
"""

import shlex
import string
from itertools import product
from typing import Any, Union

from .common import Pattern, Unformattable, extract_args_kwargs, flatten, stepout


class SequenceFormatter(string.Formatter):
    """string.Formatter subclass with special behavior for sequences.

    This class delegates the formatting of individual elements to another
    formatter object. Non-list objects are formatted by calling the
    delegate formatter's "format_field" method. List-like objects
    (list, tuple, set, frozenset) are formatted by formatting each
    element of the list according to the specified format spec using
    the delegate formatter and then joining the resulting strings with
    a separator (space by default).

    """

    def __init__(
        self, separator=" ", element_formatter=string.Formatter(), *args, **kwargs
    ):
        self.separator = separator
        self.element_formatter = element_formatter

    def format_element(self, elem, format_spec):
        """Format a single element

        For sequences, this is called once for each element in a
        sequence. For anything else, it is called on the entire
        object. It is intended to be overridden in subclases.

        """
        return self.element_formatter.format_field(elem, format_spec)

    def format_field(self, value, format_spec):
        if isinstance(value, dict) or hasattr(value, "items"):
            return ",".join(
                f"{name}={value}"
                for name, value in sorted(value.items(), key=lambda item: item[0])
            )
        if isinstance(value, (list, tuple, set, frozenset)):
            return self.separator.join(
                self.format_element(v, format_spec) for v in value
            )
        return self.format_element(value, format_spec)


class QuotedFormatter(string.Formatter):
    """Subclass of string.Formatter that supports quoting.

    Using this formatter, any field can be quoted after formatting by
    appending "q" to its format string. By default, shell quoting is
    performed using "shlex.quote", but you can pass a different
    quote_func to the constructor. The quote_func simply has to take a
    string argument and return a new string representing the quoted
    form of the input string.

    Note that if an element after formatting is the empty string, it
    will not be quoted.

    """

    def __init__(self, quote_func=shlex.quote, *args, **kwargs):
        self.quote_func = quote_func
        super().__init__(*args, **kwargs)

    def format_field(self, value, format_spec):
        if format_spec.endswith("u"):
            format_spec = format_spec[:-1]
            do_quote = False
        else:
            do_quote = format_spec.endswith("q")
            if do_quote:
                format_spec = format_spec[:-1]
        formatted = super().format_field(value, format_spec)
        if do_quote and formatted != "":
            formatted = self.quote_func(formatted)
        return formatted


class AlwaysQuotedFormatter(QuotedFormatter):
    """Subclass of QuotedFormatter that always quotes.

    Usage is identical to QuotedFormatter, except that it *always*
    acts like "q" was appended to the format spec, unless u (for unquoted) is appended.

    """

    def format_field(self, value, format_spec: str):
        if not format_spec.endswith("q") and not format_spec.endswith("u"):
            format_spec += "q"
        return super().format_field(value, format_spec)


always_qfmt = SequenceFormatter(
    separator=" ", element_formatter=AlwaysQuotedFormatter(shlex.quote)
)
qfmt = SequenceFormatter(separator=" ", element_formatter=QuotedFormatter(shlex.quote))
nfmt = string.Formatter()


class GetKeyFormatter(string.Formatter):
    class T:
        def __getattribute__(self, __name: str) -> Any:
            return self

        def __getitem__(self, __name: str) -> Any:
            return self

    def __init__(self) -> None:
        self.names: list[str] = []
        self.args: set[int] = set()
        super().__init__()

    def get_value(self, key, args, kwargs):
        if isinstance(key, int):
            self.args.add(key)
            return self.T()
        if key not in self.names:
            self.names.append(key)
        return self.T()

    @classmethod
    def get_names(cls, _pattern: str):
        self = cls()
        self.format(str(_pattern))
        return self.args, self.names


class StrPattern(Pattern):
    def __init__(self, pattern, _type=None):
        super().__init__(pattern, _type)
        self._names: list[str] = None  # type: ignore
        self._arg: set[int] = None  # type: ignore

    @property
    def names(self) -> list[str]:
        if self._names is None:
            self._arg, self._names = GetKeyFormatter.get_names(self.pattern)
        return self._names

    def select_kwargs(
        self,
        keep_missing=False,
        fill_missing=False,
        missing_value=None,
    ):
        """
        Update wildcards in given filepattern.
        Return a function that update given kwargs:
            Any args not in given filepattern will be removed;
            Any args in given filepattern will be assigned with values
            Any args being None will be filled as empty string `''`

        Arguments
        keep_missing -- whether keep the missed values
            with their values wrapped as default, i.e., `{` + key + `}`
        fill_missing -- whether to fill the missed values with
            default missing value.
            If False, even the value in the dict was given as `fill_value`,
            an Unformattable exception will be raised
        missing_value -- the value assigned as the missed value.
        """

        def _check_kwargs(*args, **kwargs):
            checked = {}
            for name in self.names:
                if name in kwargs:
                    v = kwargs[name]
                elif keep_missing:
                    v = "{" + name + "}"
                else:
                    v = missing_value
                if not fill_missing and v == missing_value:
                    raise Unformattable(name)
                if callable(v):
                    raise Unformattable(
                        f"Callable/function {v} is passed as value for {name} in 'format' statement. "
                        "This is most likely not what you want, as expand takes iterables of values or single values for "
                        "its arguments. If you want to use a function to generate the values, you can wrap the entire "
                        "expand in another function that does the computation."
                    )
                checked[name] = v if v is not None else ""
            for name in self._arg:
                checked[name] = args[name]
            return checked

        return _check_kwargs

    def formatter(
        self,
        quote_all: Union[None, bool, string.Formatter] = False,
    ):
        if quote_all is None:
            fmt = nfmt
        elif quote_all is True:
            fmt = always_qfmt
        elif quote_all is False:
            fmt = qfmt
        else:
            fmt = quote_all

        def _format(*args, **kwargs):
            used_kwargs = kwargs

            try:
                return self._type(fmt.vformat(self.pattern, args, used_kwargs))
            except KeyError as ex:
                if (
                    "wildcards" in used_kwargs
                    and str(ex).strip("'") in used_kwargs["wildcards"].keys()
                ):
                    raise Unformattable(
                        "The name '{0}' is unknown in this context. "
                        "Did you mean 'wildcards.{0}'?".format(str(ex).strip("'"))
                    )
                raise Unformattable(
                    "The name {} is unknown in this context. Please "
                    "make sure that you defined that variable. "
                    "Also note that braces not used for variable access "
                    "have to be escaped by repeating them, "
                    "i.e. {{{{print $1}}}}".format(str(ex))
                )

        return _format

    def xformat(
        self,
        nstep=1,
        keep_missing=False,
        fill_missing=False,
        missing_value=None,
        quote_all: Union[None, bool, string.Formatter] = False,
    ):
        """Format a pattern in Snakemake style.

        This means that keywords embedded in braces are replaced by any variable
        values that are available in the current namespace.
        """

        variables = stepout(nstep + (1 if nstep else 0))

        def _format(*args, **kwargs):
            used_args, used_kwargs = extract_args_kwargs(
                self.select_kwargs(
                    keep_missing=keep_missing,
                    fill_missing=fill_missing,
                    missing_value=missing_value,
                )(*args, **(variables | kwargs))
            )
            return self.formatter(quote_all=quote_all)(*used_args, **used_kwargs)

        return _format

    def xexpand(
        self,
        keep_missing=False,
        fill_missing=False,
        missing_value=None,
        combinator=product,
    ):
        """
        Expand wildcards in given filepatterns.

        Arguments
        *args -- first arg: filepatterns as list or one single filepattern,
            second arg (optional): a function to combine wildcard values
            (itertools.product per default)
        **wildcards -- the wildcards as keyword arguments
            with their values as lists. If allow_missing=True is included
            wildcards in filepattern without values will stay unformatted.
        """

        def _format(*args, **kwargs):
            # remove unused wildcards to avoid duplicate filepatterns
            used_kwargs = self.select_kwargs(
                keep_missing=keep_missing,
                fill_missing=fill_missing,
                missing_value=missing_value,
            )(*args, **kwargs)
            return [
                self.formatter(quote_all=None)(*args, **kwargs)
                for args, kwargs in map(
                    extract_args_kwargs, map(dict, combinator(*flatten(used_kwargs)))
                )
            ]

        return _format
