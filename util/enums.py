from enum import IntEnum, Enum

# Bast et al. say about edge cost:
# hop has to be >= 45 - 135
# 180 should go unpunished: 0
# 135 <= 90 <= 45
# center edges should be at least as expensive as 45 deg
# crossings edges are set to inf


class EdgePenalty(float, Enum):
    IN_SUPPORT = -1

    # Bends
    ONE_EIGHTY = 0
    ONE_THIRTY_FIVE = 1
    NINETY = 1.5
    FORTY_FIVE = 2

    # To center
    TO_CENTER = 2

    # Using any edge between ports
    # TODO with zero cost there's no need to make paths short
    # i feel like we should set this to .1 or .2 or something...
    # so that an otherwise short path with one 135deg bend isn't more expensive than a very long straight line
    HOP = 2

    CROSSING = 1000


class EdgeType(IntEnum):
    DRAW = -3
    SUPPORT = -2  # edge in hypergraph support
    PHYSICAL = -1  # physical edge that is actually drawn
    CENTER = 0  # center to center edge


class NodeType(IntEnum):
    CENTER = 0
    PORT = 1

PortDirs = ["n", "ne", "e", "se", "s", "sw", "w", "nw"]
