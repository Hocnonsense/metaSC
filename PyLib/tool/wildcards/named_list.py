# -*- coding: utf-8 -*-
"""
 * @Date: 2023-09-26 10:50:32
 * @author: Johannes Köster
 * @LastEditors: Hwrn hwrn.aou@sjtu.edu.cn
 * @LastEditTime: 2023-09-26 11:00:58
 * @FilePath: /metaSC/PyLib/tool/named_list.py
 * @Description: NamedList from snakemake
"""


from typing import Any, Callable, Optional
import functools


class NamedList(list):
    """
    A list that additionally provides functions to name items. Further,
    it is hashable, however, the hash does not consider the item names.
    """

    def __init__(
        self,
        toclone: Optional["NamedList"] = None,
        fromdict: Optional[dict[str, Any]] = None,
        custom_map: Callable = None,
    ):
        """
        Create the object.

        Arguments
        toclone  -- another Namedlist that shall be cloned
        fromdict -- a dict that shall be converted to a
            Namedlist (keys become names)
        """
        list.__init__(self)
        self._names: dict[str, Any] = dict()

        self._allowed_overrides = ["index", "sort"]
        # white-list of attribute names that can be overridden in _set_name
        # default to throwing exception if called to prevent use as functions
        for name in self._allowed_overrides:
            setattr(self, name, functools.partial(self._used_attribute, _name=name))

        if toclone is not None:
            if custom_map is not None:
                self.extend(map(custom_map, toclone))
            else:
                self.extend(toclone)
            if isinstance(toclone, NamedList):
                self._take_names(toclone._get_names())
        if fromdict is not None:
            for key, item in fromdict.items():
                self.append(item)
                self._add_name(key)

    @staticmethod
    def _used_attribute(*args, _name, **kwargs):
        """
        Generic function that throws an `AttributeError`.

        Used as replacement for functions such as `index()` and `sort()`,
        which may be overridden by workflows, to signal to a user that
        these functions should not be used.
        """
        raise AttributeError(
            f"{_name}() cannot be used; attribute name reserved"
            " for use in some existing workflows"
        )

    def _add_name(self, name: str):
        """
        Add a name to the last item.

        Arguments
        name -- a name
        """
        self._set_name(name, len(self) - 1)

    def _set_name(self, name: str, index, end=None):
        """
        Set the name of an item.

        Arguments
        name  -- a name
        index -- the item index
        """
        if name not in self._allowed_overrides and hasattr(self.__class__, name):
            raise AttributeError(
                "invalid name: {name} is reserved for internal use".format(name=name)
            )

        self._names[name] = (index, end)
        if end is None:
            setattr(self, name, self[index])
        else:
            setattr(self, name, NamedList(toclone=self[index:end]))

    def _get_names(self):
        """
        Get the defined names as (name, index) pairs.
        """
        for name, index in self._names.items():
            yield name, index

    def _take_names(self, names):
        """
        Take over the given names.

        Arguments
        names -- the given names as (name, index) pairs
        """
        for name, (i, j) in names:
            self._set_name(name, i, end=j)

    def items(self):
        for name in self._names:
            yield name, getattr(self, name)

    def _allitems(self):
        next = 0
        for name, index in sorted(
            self._names.items(),
            key=lambda item: (
                item[1][0],
                item[1][0] + 1 if item[1][1] is None else item[1][1],
            ),
        ):
            start, end = index
            if end is None:
                end = start + 1
            if start > next:
                for item in self[next:start]:
                    yield None, item
            yield name, getattr(self, name)
            next = end
        for item in self[next:]:
            yield None, item

    def _insert_items(self, index, items):
        self[index : index + 1] = items
        add = len(items) - 1
        for name, (i, j) in self._names.items():
            if i > index:
                self._names[name] = (i + add, None if j is None else j + add)
            elif i == index:
                self._set_name(name, i, end=i + len(items))

    def keys(self):
        return self._names.keys()

    def _plainstrings(self):
        return self.__class__.__call__(toclone=self, plainstr=True)

    def _stripped_constraints(self):
        return self.__class__.__call__(toclone=self, strip_constraints=True)

    def _clone(self):
        return self.__class__.__call__(toclone=self)

    def get(self, key, default_value=None):
        return self.__dict__.get(key, default_value)

    def __getitem__(self, key):
        if isinstance(key, str):
            return getattr(self, key)
        else:
            return super().__getitem__(key)

    def __hash__(self):
        return hash(tuple(self))

    def __str__(self):
        return " ".join(map(str, self))
