from enum import IntEnum, Enum

# Bast et al. say about edge cost:
# hop has to be >= 45 - 135
# 180 should go unpunished: 0
# 135 <= 90 <= 45
# center edges should be at least as expensive as 45 deg
# crossings edges are set to inf


class EdgePenalty(float, Enum):
    IN_SUPPORT = -1
    COMMON_MULTILAYER = -0.5

    # Bends
    ONE_EIGHTY = 0
    ONE_THIRTY_FIVE = 1
    NINETY = 2
    FORTY_FIVE = 3

    # To center
    TO_CENTER = 3

    # Using any edge between ports
    # TODO with zero cost there's no need to make paths short
    # i feel like we should set this to .1 or .2 or something...
    # so that an otherwise short path with one 135deg bend isn't more expensive than a very long straight line
    HOP = 2

    # crossings relative to node
    # diagonal edges between different nodes
    # there aren't too many in the grid graph of those so the effect is meh
    CROSSING_OUTSIDE = 1000

    # at port edge of an occupied node
    # would not recommend to set that because we often can't work around them anyways
    CROSSING_INSIDE_GLYPH = 0

    # at port edge of an unoccupied node
    # attention when setting this one
    # if port edges of the same node are penalized it tries to go over the center instead
    # so to avoid that, set the center penalty to the same value
    # attention #2: setting this also prevents parallel lines unless they connect exactly the same nodes
    # reason being that because of the former reason, a turn away from an otherwise parallel path is considered a crossing
    CROSSING_INSIDE_CELL = 0


class EdgeType(IntEnum):
    SUPPORT = -2  # edge in hypergraph support
    PHYSICAL = -1  # physical edge that is actually drawn
    CENTER = 0  # center to center edge


class NodeType(IntEnum):
    CENTER = 0
    PORT = 1


PortDirs = ["n", "ne", "e", "se", "s", "sw", "w", "nw"]
