import unittest
from itertools import product

from util.geometry import do_lines_intersect


class TestLineSegmentIntersect(unittest.TestCase):
    def test_regular_case(self):
        a = (1, 1)
        b = (3, 3)
        c = (1, 3)
        d = (3, 1)

        self.assertTrue(do_lines_intersect(a, b, c, d))

    def test_same_line(self):
        a = (1,1)
        b = (4,2)
        c = (1,1)
        d = (4,2)

        self.assertFalse(do_lines_intersect(a,b,c,d))

    def test_specific(self):
        a = (4,4)
        b = (8,6)
        c = (6,4)
        d = (6,7)

        self.assertTrue(do_lines_intersect(a,b,c,d))
        self.assertTrue(do_lines_intersect(b,a,c,d))
        self.assertTrue(do_lines_intersect(a,b,d,c))
        self.assertTrue(do_lines_intersect(b,a,d,c))

    def test_collinear(self):
        a = (1, 1)
        b = (4, 4)
        c = (2, 2)
        d = (3, 3)

        self.assertFalse(do_lines_intersect(a, b, c, d))

        a = (1, 1)
        b = (3, 3)
        c = (2, 2)
        d = (4, 4)

        self.assertFalse(do_lines_intersect(a, b, c, d))

    def test_point_on_line(self):
        a = (1,1)
        b = (4,4)
        c = (4,1)
        d = (3,3)

        self.assertFalse(do_lines_intersect(a,b,c,d))
        self.assertFalse(do_lines_intersect(d,c,a,b))

    def test_shared_point(self):
        a = (1, 1)
        b = (4, 4)
        c = (1, 1)
        d = (8, 6)

        self.assertFalse(do_lines_intersect(a, b, c, d))

    def test_parallel(self):
        a = (1, 1)
        b = (4, 1)
        c = (1, 2)
        d = (4, 2)

        self.assertFalse(do_lines_intersect(a, b, c, d))


if __name__ == "__main__":
    unittest.main()
