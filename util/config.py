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

SUB_SUPPORT_TYPE = ContextVar(
    'SUB_SUPPORT_TYPE', default = config.get('ROUTE', 'SubSupportType')
)
SUB_SUPPORT_GROUPING = ContextVar(
    'SUB_SUPPORT_GROUPING', default = config.get('ROUTE', 'SubSupportGrouping')
)

STRATEGY = ContextVar(
    'STRATEGY', default= config.get('GENERAL', 'Strategy')
)


GRID_WIDTH = ContextVar(
    'GRID_WIDTH', default= config.getint('GENERAL', 'GridWidth')
)

GRID_HEIGHT = ContextVar(
    'GRID_HEIGHT', default= config.getint('GENERAL', 'GridHeight')
)
