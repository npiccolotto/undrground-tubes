import math
import numpy as np
from collections import Counter
from pygame.math import Vector2


def dist_euclidean(p1, p2):
    # p1,p2 are tuples
    x1, y1 = p1
    x2, y2 = p2
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def get_side(a, b, c):
    "Given a line from a to b in 2D and a point c, determine if c is left, right or on the line."
    ax, ay = a
    bx, by = b
    cx, cy = c
    det = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    if det > 0:
        # left
        return -1
    if det < 0:
        # right
        return 1
    return 0


def get_closest_point(reference, options):
    dists = [dist_euclidean(reference, option) for option in options]
    closest = np.argmin(np.array(dists))
    return options[closest]


def centroid(points):
    x, y = zip(*points)
    n = len(x)
    return (sum(x) / n, sum(y) / n)


def are_faces_adjacent(face1, face2):
    # faces are a list of tuples (points)
    # faces are adjacent obv iff they share an edge, ie, two consecutive points
    # however, the way our grid is constructed, any two shared points must be a shared edge
    # thus we count the occurrences of points in both faces and if there are 2 points appearing 2 times, that's an edge
    c = Counter(face1 + face2)
    counter = Counter(list(c.values()))
    return counter[2] == 2


def triangular_area(a, b, c):
    aX, aY = a
    bX, bY = b
    cX, cY = c
    return (bX - aX) * (cY - aY) - (cX - aX) * (bY - aY)


# https://bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
def ccw(a, b, c):
    ax, ay = a
    bx, by = b
    cx, cy = c
    return (cy - ay) * (bx - ax) > (by - ay) * (cx - ax)


def do_lines_intersect(a, b, c, d):
    return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def get_angle(u, v):
    ux, uy = u
    vx, vy = v
    theta = math.atan2(uy - vy, ux - vx)
    if theta < 0:
        theta += 2 * math.pi
    return theta


def offset_point(p, angle, length):
    pn = np.array(p) + np.array(pol2cart(length, angle))
    x, y = pn
    return (x, y)


def offset_edge(edge, angle, length):
    s, t = edge
    return (offset_point(s, angle, length), offset_point(t, angle, length))


def is_point_inside_circle(point, circle):
    center, radius = circle
    d = dist_euclidean(point, center)
    return d <= radius


def get_segment_circle_intersection(segment, circle):
    # https://codereview.stackexchange.com/a/86428
    s1, s2 = segment
    ax, ay = s1
    bx, by = s2
    center, r = circle
    cx, cy = center

    Q = Vector2(cx, cy)  # Centre of circle
    P1 = Vector2(ax, ay)  # Start of line segment
    V = Vector2(bx, by) - P1  # Vector along line segment

    a = V.dot(V)
    b = 2 * V.dot(P1 - Q)
    c = P1.dot(P1) + Q.dot(Q) - 2 * P1.dot(Q) - r**2

    disc = b**2 - 4 * a * c
    if disc < 0:
        return None

    sqrt_disc = math.sqrt(disc)
    t1 = (-b + sqrt_disc) / (2 * a)
    t2 = (-b - sqrt_disc) / (2 * a)
    if not (0 <= t1 <= 1 or 0 <= t2 <= 1):
        return None

    # https://stackoverflow.com/a/1084899/490524
    if t1 >= 0 and t1 <= 1:
        i = P1 + t1 * V
        return i
    if t2 >= 0 and t2 <= 1:
        i = P1 + t2 * V
        return i

    return None


def simple_control_points(a, b, c, offset=2):
    "Given three points on polyline abc, return control points for middle point"
    # idea is following. determine angle between a and c
    # control points are offset from b in that angle and 2π-angle
    a = get_angle(a, c)
    c1 = offset_point(b, a, offset)
    a_ = 2 * math.pi - a
    if a_ < 0:
        a_ += 2 * math.pi
    c2 = offset_point(b, a_, offset)
    return (c1, c2)


# https://stackoverflow.com/a/50974391/490524
def define_circle(p1, p2, p3):
    """
    Returns the center and radius of the circle passing the given 3 points.
    In case the 3 points form a line, returns (None, infinity).
    """
    temp = p2[0] * p2[0] + p2[1] * p2[1]
    bc = (p1[0] * p1[0] + p1[1] * p1[1] - temp) / 2
    cd = (temp - p3[0] * p3[0] - p3[1] * p3[1]) / 2
    det = (p1[0] - p2[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p2[1])

    if abs(det) < 1.0e-6:
        return (None, np.inf)

    # Center of circle
    cx = (bc * (p2[1] - p3[1]) - cd * (p1[1] - p2[1])) / det
    cy = ((p1[0] - p2[0]) * cd - (p2[0] - p3[0]) * bc) / det

    radius = np.sqrt((cx - p1[0]) ** 2 + (cy - p1[1]) ** 2)
    return ((cx, cy), radius)
