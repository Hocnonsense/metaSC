# -*- coding: utf-8 -*-
"""
   @Description: 仅能在一定区域内运动的点的类, 上下, 左右联通

   @Author: Hwrn
   @Date: 2020-02-05 15:07:03
   @LastEditors: Hwrn
   @LastEditTime: 2020-02-05 17:25:39

   @TODU: to test the file
"""


class Boundary(object):
    """
    location of each point
    all Points move inside __XBoundary and __YBoundary
    """

    def __init__(self, xboundary, yboundary):
        self.__xboundary, self.__yboundary = xboundary, yboundary

    def __call__(self, x, y):
        # pylint: disable = invalid-name
        return Point(x, y, self)

    @property
    def x(self):
        """X length"""
        return self.__xboundary

    @property
    def y(self):
        """Y length"""
        return self.__yboundary


def istuple(tuple_x, _y=None):
    """
    to suit the situation if return a tuple
    """
    # pylint: disable = invalid-name
    try:
        x, y = tuple_x
    except TypeError:
        x, y = tuple, _y
    return (x, y)


class Point(object):
    """
    Point Class.
    """
    def __init__(self, x, y, boundary: Boundary = None):
        """
        @param boundary The boundary of Point.
                        Default None for endless ground
        """
        self.__x, self.__y = self.__ifTuple(x, y)
        self.boundary = boundary

    #X property
    def __getx(self):
        return self.__x

    def __setx(self, _x):
        self.__x = (_x % self.boundary.x) if self.boundary else _x
    x = property(__getx, __setx)

    # Y property
    def __gety(self):
        return self.__y

    def __sety(self, _y):
        self.__y = (_y % self.boundary.y) if self.boundary else _y
    y = property(__gety, __sety)

    def __str__(self):
        return "{X:" + "{:.0f}".format(self.__x) + \
               ",Y:" + "{:.0f}".format(self.__y) + "}" + \
               (" in [" + "{:.0f}".format(self.boundary.x) + \
                "," + "{:.0f}".format(self.boundary.y) + "]") if self.boundary else ""

    def __call__(self, d_x=0, d_y=0, absolute = False):
        """
        Go to where it would goto.  Return where it was
        @param absolute If location is absolutely.  Default use `move(d_x, d_y)`
                        `True` not use `move`
        """
        _x, _y = self.__x, self.__y
        if absolute:
            self.x, self.y = d_x, d_y
        else:
            self.__x, self.__y = self.move(d_x, d_y)
        return _x, _y

    def move(self, d_x=0, d_y=0):
        """
        Try to move to somewhere relative to location now.
        Return where it would be
        """
        return (self.x + d_x) % self.boundary.x if self.boundary else self.x + d_x, \
               (self.y + d_y) % self.boundary.y if self.boundary else self.y + d_y

    def around(self):
        """
        @Return things arrounding it.
        1 2 3
        4 * 5
        6 7 8
        """
        neighbors = list()
        for d_y in [-1, 0, 1]:
            for d_x in [-1, 0, 1]:
                _x, _y = self.move(d_x, d_y)
                neighbors.append((_x, _y))
        neighbors.remove((self.x, self.y))
        return neighbors
