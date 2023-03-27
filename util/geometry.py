import math
import numpy as np
from collections import Counter


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
    a, b = segment
    c, r = circle
    # crude, but should work:
    # if a or b is less than r away from c, we have an intersection
    d_ac = dist_euclidean(a, c)
    d_bc = dist_euclidean(b, c)
    if d_ac > r and d_bc > r:
        # in the general case both points could be outside the circle and yet the segment
        # still intersecting it, but due to our construction rules that can't really happen
        return None
    if d_ac <= r and d_bc <= r:
        raise Exception("somehow segment completely in circle")
    # the intersection is found by offsetting from the point inside the circle
    point_in_circle = a if d_ac < d_bc else b
    point_out_circle = b if point_in_circle == a else a
    offset_from_point_in_circle = min(d_ac, d_bc) - r
    angle = get_angle(point_in_circle, point_out_circle)
    return offset_point(point_in_circle, angle, offset_from_point_in_circle)
