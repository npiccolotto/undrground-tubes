import configparser
import numpy as np
from contextvars import ContextVar

config = configparser.ConfigParser()
config.read("config.ini")


DARK_MODE = ContextVar("DARK_MODE", default=config.getboolean("DRAW", "DarkMode"))
DRAW_GLYPHS = ContextVar("DRAW_GLYPHS", default=config.getboolean("DRAW", "DrawGlyphs"))
DRAW_GLYPHS_OVER_LINES = ContextVar(
    "DRAW_GLYPHS", default=config.getboolean("DRAW", "DrawGlyphsOverLines")
)
STROKE_WIDTH = ContextVar(
    "STROKE_WIDTH", default=config.getfloat("DRAW", "StrokeWidth")
)
LINE_GAP = ContextVar("LINE_GAP", default=config.getfloat("DRAW", "LineGap"))
CELL_SIZE_PX = ContextVar("CELL_SIZE_PX", default=config.getfloat("DRAW", "CellSizePx"))
GLYPH_SIZE_PX = ContextVar(
    "GLYPH_SIZE_PX", default=config.getfloat("DRAW", "GlyphSizePx")
)
MARGINS = ContextVar(
    "MARGINS",
    default=np.array(
        (
            config.getfloat("DRAW", "MarginHorizontal"),
            config.getfloat("DRAW", "MarginVertical"),
        )
    ),
)
NODE_CIRCLE_RADIUS = ContextVar(
    "NODE_CIRCLE_RADIUS", default=config.getfloat("DRAW", "NodeCircleRadius")
)
DRAW_DEG1_MARKS = ContextVar(
    "DRAW_DEG1_MARKS", default=config.getboolean("DRAW", "DrawDeg1Marks")
)
DRAW_DEG1_MARK_SIZE_PX = ContextVar(
    "DRAW_DEG1_MARK_SIZE_PX", default=config.getfloat("DRAW", "Deg1MarkSize")
)
DRAW_HUBS = ContextVar("DRAW_HUBS", default=config.getboolean("DRAW", "DrawHubs"))
DRAW_GRID_CELLS = ContextVar(
    "DRAW_GRID_CELLS", default=config.getboolean("DRAW", "DrawGridCells")
)
SET_COLORS = ContextVar(
    "SET_COLORS", default=config.get("DRAW", "SetColors").split(",")
)
GLYPH_TITLE = ContextVar(
    'GLYPH_TITLE', default = config.getboolean('DRAW', 'GlyphTitle')
)

SUB_SUPPORT_TYPE = ContextVar(
    "SUB_SUPPORT_TYPE", default=config.get("ROUTE", "SubSupportType")
)
SUB_SUPPORT_GROUPING = ContextVar(
    "SUB_SUPPORT_GROUPING", default=config.get("ROUTE", "SubSupportGrouping")
)

NUM_WEIGHTS = ContextVar("GRID_WIDTH", default=config.getint("GENERAL", "NumLayers"))
STRATEGY = ContextVar("STRATEGY", default=config.get("GENERAL", "Strategy"))
GRID_WIDTH = ContextVar("GRID_WIDTH", default=config.getint("GENERAL", "GridWidth"))
GRID_HEIGHT = ContextVar("GRID_HEIGHT", default=config.getint("GENERAL", "GridHeight"))

READ_DIR = ContextVar("READ_DIR", default=config.get("GENERAL", "ReadDir"))
WRITE_DIR = ContextVar("WRITE_DIR", default=config.get("GENERAL", "WriteDir"))


EDGE_SOFT_CONSTRAINT_WEIGHT = ContextVar(
    "EDGE_SOFT_CONSTRAINT_WEIGHT", default=config.getint("ROUTE", "EdgeLengthFactor")
)
EDGE_LAYER_SOFT_CONSTRAINT_WEIGHT = ContextVar(
    "EDGE_LAYER_SOFT_CONSTRAINT_WEIGHT",
    default=config.getint("ROUTE", "EdgeLayerLengthFactor"),
)
BEND_SOFT_CONSTRAINT_WEIGHT = ContextVar(
    "BEND_SOFT_CONSTRAINT_WEIGHT", default=config.getint("ROUTE", "BendLayerFactor")
)
ILP_MIP_GAP = ContextVar(
    "ILP_MIP_GAP", default=config.getfloat("ROUTE", "ILPMIPGap")
)
ILP_TIMEOUT =  ContextVar(
    "ILP_TIMEOUT", default=config.getfloat("ROUTE", "ILPTimeoutSecs")
)

LOOM_SOLVER = ContextVar(
    'LOOM_SOLVER', default=config.get('LOOM', 'Solver')
)
LOOM_TIMEOUT = ContextVar(
    'LOOM_TIMEOUT', default=config.getint('LOOM', 'TimeoutSecs')
)
