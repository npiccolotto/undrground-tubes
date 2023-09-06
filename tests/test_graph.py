import unittest
from itertools import product

from util.graph import are_port_edges_crossing
from util.enums import PortDirs


class TestCrossingPortEdges(unittest.TestCase):
    def test_commutative(self):
        u = {"port": "nw"}
        v = {"port": "s"}
        w = {"port": "n"}
        x = {"port": "e"}

        self.assertEqual(
            are_port_edges_crossing(u, v, w, x),
            are_port_edges_crossing(v, u, w, x),
            "first pair switched",
        )
        self.assertEqual(
            are_port_edges_crossing(u, v, x, w),
            are_port_edges_crossing(u, v, x, w),
            "second pair switched",
        )

    def test_no_crossing_no_overlap(self):
        u = {"port": "nw"}
        v = {"port": "s"}
        for w, x in [("n", "e"), ("n", "ne"), ("sw", "w")]:
            self.assertFalse(
                are_port_edges_crossing(u, v, {"port": w}, {"port": x}), (w, x)
            )

    def test_no_crossing_overlap(self):
        u = {"port": "nw"}
        v = {"port": "s"}
        for w, x in [("nw", "e"), ("n", "s"), ("s", "nw")]:
            self.assertFalse(
                are_port_edges_crossing(u, v, {"port": w}, {"port": x}), (w, x)
            )

    def test_crossing(self):
        u = {"port": "nw"}
        v = {"port": "s"}
        crossings = list(product(['sw', 'w'], ['n', 'ne', 'e', 'se']))
        for w, x in crossings:
            self.assertTrue(
                are_port_edges_crossing(u, v, {"port": w}, {"port": x}), (w, x)
            )


if __name__ == "__main__":
    unittest.main()
