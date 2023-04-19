import math
import numpy as np
from collections import Counter
from itertools import pairwise
from pygame.math import Vector2
import drawsvg as svg


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
    # Collinear lines count as intersecting or when they share a point
    return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)


def do_lines_intersect_strict(a, b, c, d):
    # collinear lines and idential points do NOT count as intersection

    abc = get_side(a, b, c)
    abd = get_side(a, b, d)
    cda = get_side(c, d, a)
    cdb = get_side(c, d, b)

    ab_separates_cd = abc != abd and abc != 0 and abd != 0
    cd_separates_ab = cda != cdb and cda != 0 and cdb != 0

    return ab_separates_cd and cd_separates_ab


def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return (rho, phi)


def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x, y)


def get_angle_vector(v):
    angle = math.atan2(v[1], v[0])

    if angle < 0:
        angle += 2 * math.pi

    return angle


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


def is_point_on_circle(point, circle, eps=10e-4):
    center, radius = circle
    d = dist_euclidean(point, center)
    return abs(d - radius) < eps


def get_line_segment_intersection(segment1, segment2):
    # http://paulbourke.net/geometry/pointlineplane/
    a, b = segment1
    c, d = segment2

    x1, y1 = a
    x2, y2 = b
    x3, y3 = c
    x4, y4 = d

    if (x1 == x2 and y1 == y2) or (x3 == x4 and y3 == y4):
        return [x2, y2]

    denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)

    # Lines are parallel
    if denominator == 0:
        return [x2, y2]

    ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator
    ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator

    x = x1 + ua * (x2 - x1)
    y = y1 + ua * (y2 - y1)

    return (x, y)


def get_segment_circle_intersection(segment, circle):
    # https://codereview.stackexchange.com/a/86428
    # TODO get rid of pygame vectors and convert to numpy
    # (it's probably numpy under the hood anyways)
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


def biarc_segment(p, q):
    pass


def are_vectors_parallel(a, b, eps=10e-4):
    ax, ay = a
    bx, by = b
    return abs(ax * by - bx * ay) < eps


def to_unit(v):
    n = np.linalg.norm(v)
    if n != 0:
        return v / n
    return v


def vector(a, b, unit=False):
    ax, ay = a
    bx, by = b
    v = np.array((bx - ax, by - ay))
    return to_unit(v) if unit else v


def normal_vector(a):
    x, y = a
    return np.array((-y, x))


def is_clockwise_order(a0, a1):
    # assumed both to be in [0, 2π]
    if a0 < a1:
        return a1 - a0 > math.pi

    return a1 + 2 * math.pi - a0 > math.pi


def get_svg_arc_params(
    startpoint, endpoint, center, radius, angle1, angle2, clockwise, eps=10e-4
):
    counterclockwise = not clockwise
    # we have a startpoint, from which the arc should originate
    # then a circle given by center and radius
    # the arc should be from angle1 to angle2 on that circle
    angle_diff = angle1 - angle2 if clockwise else angle2 - angle1

    # is the startpoint on the arc? if not we screwed up
    cx, cy = center
    sx, sy = startpoint
    sx_c = cx + radius * math.cos(angle1)
    sy_c = cy + radius * math.sin(angle1)

    if abs(sx - sx_c) > eps or abs(sy - sy_c) > eps:
        raise Exception("startpoint not on arc circle circumference")

    ex, ey = (cx + radius * math.cos(angle2), cy + radius * math.sin(angle2))
    if abs(endpoint[0] - ex) > eps or abs(endpoint[1] - ey) > eps:
        raise Exception("desired endpoint not the same as calculated endpoint'")

    large_arc_bit = angle_diff >= math.pi
    sweep_bit = counterclockwise

    return (ex, ey, radius, sweep_bit, large_arc_bit)


def interpolate_biarcs(points, **kwargs):
    if len(points) < 2:
        return None

    line = svg.Path(**kwargs)

    pairs = pairwise(points)
    for i, (p, q) in enumerate(pairs):
        seg_type, p = p
        _, q = q

        if seg_type == "L":
            line.M(p[0], p[1]).L(q[0], q[1])
            continue
        if seg_type != "B":
            raise Exception(f"unknown segment type '{seg_type}'")

        _, ps = points[i - 1]
        _, p0 = points[i]
        _, p4 = points[i + 1]
        _, pt = points[i + 2]

        # "Biarc approximation of NURBS curves", Piegl & Tiller, 2002
        tangent_s = vector(ps, p0, unit=True)
        tangent_e = vector(p4, pt, unit=True)

        if are_vectors_parallel(tangent_s, tangent_e) or (
            tangent_s.dot(p4 - p0) <= 0 and tangent_s.dot(tangent_e) <= 0
        ):
            # tangents are parallel (or some other bad case happened), but not necessarily collinear
            # so we draw a bezier curve
            d = dist_euclidean(p0, p4)
            cx1, cy1 = np.array(p0) + tangent_s * d / 4
            cx2, cy2 = np.array(p4) - tangent_e * d / 4

            line.C(cx1, cy1, cx2, cy2, p4[0], p4[1])
            continue

        v = p0 - p4
        a = 2 * (np.dot(tangent_s, tangent_e) - 1)
        b = 2 * np.dot(v, tangent_s + tangent_e)
        c = np.dot(v, v)

        alpha1 = (-b + math.sqrt(b**2 - 4 * a * c)) / (2 * a)
        alpha2 = (-b - math.sqrt(b**2 - 4 * a * c)) / (2 * a)

        alpha = alpha2 if alpha1 < 0 else alpha1

        p1 = p0 + tangent_s * alpha
        p3 = p4 + tangent_e * (-alpha)
        p2 = 1 / 2 * (p1 + p3)

        p1p0_norm = normal_vector(vector(p1, p0))
        p3p2_norm = normal_vector(vector(p3, p2))
        p3p4_norm = normal_vector(vector(p3, p4))

        pointOnP1P0Norm = p0 + p1p0_norm
        pointOnP3P2Norm = p2 + p3p2_norm
        pointOnP3P4Norm = p4 + p3p4_norm

        c1 = get_line_segment_intersection((p0, pointOnP1P0Norm), (p2, pointOnP3P2Norm))
        c2 = get_line_segment_intersection((p4, pointOnP3P4Norm), (p2, pointOnP3P2Norm))

        r1 = dist_euclidean(p0, c1)
        r2 = dist_euclidean(p4, c2)

        assert is_point_on_circle(p0, (c1, r1))
        assert is_point_on_circle(p2, (c1, r1))
        assert is_point_on_circle(p2, (c2, r2))
        assert is_point_on_circle(p4, (c2, r2))

        c1p0 = vector(c1, p0)
        c1p2 = vector(c1, p2)
        angle1a = get_angle_vector(c1p0)
        angle1b = get_angle_vector(c1p2)
        clockwise1 = is_clockwise_order(angle1a, angle1b)

        c2p2 = vector(c2, p2)
        c2p4 = vector(c2, p4)
        angle2a = get_angle_vector(c2p2)
        angle2b = get_angle_vector(c2p4)
        clockwise2 = is_clockwise_order(angle2a, angle2b)

        line.M(p0[0], p0[1])
        arc1_x, arc1_y, arc1_r, arc1_sweep, arc1_large = get_svg_arc_params(
            p0, p2, c1, r1, angle1a, angle1b, clockwise1
        )
        line.A(arc1_r, arc1_r, 0, arc1_large, arc1_sweep, arc1_x, arc1_y)
        line.M(arc1_x, arc1_y)
        arc2_x, arc2_y, arc2_r, arc2_sweep, arc2_large = get_svg_arc_params(
            p2, p4, c2, r2, angle2a, angle2b, clockwise2
        )
        line.A(arc2_r, arc2_r, 0, arc2_large, arc2_sweep, arc2_x, arc2_y)
        line.M(arc2_x, arc2_y)

    return line


def logical_coords_to_physical(x, y, lattice_type="sqr"):
    if lattice_type == "hex":
        if y % 2 == 1:
            return (x + 0.5, -y)
    return (x, -y)
