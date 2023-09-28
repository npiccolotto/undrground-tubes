import unittest
from itertools import product

from util.graph import are_port_edges_crossing, get_relative_ports,remove_segment_in_circle
from util.enums import PortDirs


class TestCycleBreaker(unittest.TestCase):
    def test_circle_breaker(self):
        self.assertListEqual(
            remove_segment_in_circle(list(range(5)) + [0], 1, 3), [3, 4, 0, 1]
        )


class TestGetPorts(unittest.TestCase):
    def test_getrelativeports(self):
        self.assertEqual(get_relative_ports((1, 1), (2, 2)), ("nw", "se"))
        self.assertEqual(get_relative_ports((1, 1), (2, 1)), ("w", "e"))
        self.assertEqual(get_relative_ports((1, 1), (1, 2)), ("n", "s"))
        self.assertEqual(get_relative_ports((1, 2), (1, 1)), ("s", "n"))


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

    def test_crossing_overlap(self):
        self.assertTrue(
            are_port_edges_crossing(
                {"port": "e"}, {"port": "se"}, {"port": "e"}, {"port": "w"}
            )
        )

        u = {"port": "nw"}
        v = {"port": "s"}
        for w, x in [("nw", "e"), ("n", "s")]:
            self.assertTrue(
                are_port_edges_crossing(u, v, {"port": w}, {"port": x}), (w, x)
            )

        for w, x in [("nw", "e"), ("n", "s")]:
            self.assertFalse(
                are_port_edges_crossing(
                    u, v, {"port": w}, {"port": x}, cross_when_node_shared=False
                ),
                (u, v, w, x),
            )

    def test_same_edge(self):
        self.assertFalse(
            are_port_edges_crossing(
                {"port": "e"}, {"port": "se"}, {"port": "se"}, {"port": "e"}
            )
        )

    def test_crossing(self):
        u = {"port": "nw"}
        v = {"port": "s"}
        crossings = list(product(["sw", "w"], ["n", "ne", "e", "se"]))
        for w, x in crossings:
            self.assertTrue(
                are_port_edges_crossing(u, v, {"port": w}, {"port": x}), (w, x)
            )


if __name__ == "__main__":
    unittest.main()
