import copy
import math
import random
from typing import List, Tuple, Any

import matplotlib.pyplot as plt
import numpy as np
from shapely.affinity import rotate, translate
from shapely.geometry import Polygon, Point, LineString
from shapely.geometry.base import BaseGeometry


def intersection_of_rotated_polygon(angle: float, polygon_A: Polygon, polygon_B: Polygon,
                                    pivot_point: Point) -> float:
    """Calculates the intersection area as a function of the rotation angle.

    Args:
        angle: The rotation angle in radians.
        polygon_A: The fixed polygon.
        polygon_B: The rotating polygon.
        pivot_point: The pivot point of the rotating polygon.

    Returns:
        The intersection area.
    """

    rotated_B = rotate(polygon_B, angle, origin=pivot_point, use_radians=True)
    return polygon_A.intersection(rotated_B).area


def plot_intersection_area(angles: List[float], areas: List[float]) -> None:
    """Plots the intersection area as a function of the rotation angle.

    Args:
        angles: A list of rotation angles.
        areas: A list of corresponding intersection areas.
    """
    plt.plot(angles, areas)
    plt.xlabel("Rotation Angle (α)")
    plt.ylabel("Intersection Area")
    plt.show()


def generate_random_convex_polygon(num_vertices: int, bounds: Tuple[float, float, float, float]) -> Polygon:
    """
    Generates a random convex Shapely Polygon with a specified number of vertices within given bounds.

    Args:
        num_vertices (int): Number of vertices for the polygon (must be >= 3).
        bounds (Tuple[float, float, float, float]): The (min_x, max_x, min_y, max_y) bounds for the polygon vertices.

    Returns:
        Polygon: A randomly generated convex Shapely Polygon.
    """
    if num_vertices < 3:
        raise ValueError("A polygon must have at least 3 vertices.")

    min_x, max_x, min_y, max_y = bounds

    # Generate random points in x and y separately and sort them
    x_coords = sorted(random.uniform(min_x, max_x) for _ in range(num_vertices))
    y_coords = sorted(random.uniform(min_y, max_y) for _ in range(num_vertices))

    # Build the polygon using a method that guarantees convexity
    lower_hull = [(x, y) for x, y in zip(x_coords[:len(x_coords) // 2], y_coords[:len(y_coords) // 2])]
    upper_hull = [(x, y) for x, y in
                  zip(reversed(x_coords[len(x_coords) // 2:]), reversed(y_coords[len(y_coords) // 2:]))]

    # Combine hulls to form the full convex polygon
    points = lower_hull + upper_hull

    # Return a Shapely polygon
    return Polygon(points).convex_hull


def is_polygon_convex(polygon: Polygon) -> bool:
    """Check if a Shapely Polygon is convex."""
    coords = list(polygon.exterior.coords)[:-1]  # Exclude the closing coordinate

    def cross_product(o, a, b):
        """Compute the cross product of vectors OA and OB (O is the origin)."""
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    is_positive = None
    for i in range(len(coords)):
        o, a, b = coords[i], coords[(i + 1) % len(coords)], coords[(i + 2) % len(coords)]
        cross = cross_product(o, a, b)
        if cross != 0:
            current_sign = cross > 0
            if is_positive is None:
                is_positive = current_sign
            elif is_positive != current_sign:
                return False
    return True


def calculate_full_rotation_overlap(polygon: Polygon, pivot: Point, offset: Point) -> dict[float, float]:
    intersection_areas: dict[float, float] = {}
    fixed_polygon = copy.deepcopy(polygon)

    rotated_polygon = translate(polygon, offset.x, offset.y)

    for angle in np.linspace(0.0, 2 * np.pi, 360):
        intersection_area: float = intersection_of_rotated_polygon(angle, fixed_polygon, rotated_polygon,
                                                                   pivot_point=pivot)
        intersection_areas[angle] = intersection_area

    return intersection_areas


def find_min_max_points(values: list[Any]) -> list[tuple[int, float]]:
    if len(values) < 3:
        return []

    min_max_points = []
    for i in range(1, len(values) - 1):
        if values[i] > values[i - 1] and values[i] > values[i + 1]:
            min_max_points.append((i, values[i]))  # Local maximum
        elif values[i] < values[i - 1] and values[i] < values[i + 1]:
            min_max_points.append((i, values[i]))  # Local minimum

    return min_max_points


def find_extremes(polygon: Polygon, offset: Point, pivot:Point) -> list:

    overlaps = calculate_full_rotation_overlap(polygon, pivot=pivot, offset=offset)
    overlaps = np.round(list(overlaps.values()), 10)
    extremes = find_min_max_points(overlaps)
    return extremes


def find_long_diagonal(polygon: Polygon) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Finds the longest segment between any two vertices of the polygon.

    For a kite this is the major diagonal.

    Args:
        polygon: The polygon to inspect.

    Returns:
        The two vertices spanning the longest diagonal.
    """
    coords = list(polygon.exterior.coords)[:-1]

    best_pair = (coords[0], coords[1])
    best_length = -1.0
    for i in range(len(coords)):
        for j in range(i + 1, len(coords)):
            length = math.dist(coords[i], coords[j])
            if length > best_length:
                best_length = length
                best_pair = (coords[i], coords[j])

    return best_pair


def _cross_section_length(polygon: Polygon, u: np.ndarray, v: np.ndarray, s: float) -> float:
    """Length of the polygon's cross-section on the line perpendicular to `u` at offset `s`."""
    coords = np.array(polygon.exterior.coords)
    half_span = float(np.abs(coords @ v).max()) + 1.0

    base = s * u
    cut = LineString([base - half_span * v, base + half_span * v])
    return polygon.intersection(cut).length


def _slice_boundaries(polygon: Polygon, u: np.ndarray, v: np.ndarray, thickness: float) -> List[float]:
    """Computes the slice boundary offsets along `u`.

    Slicing starts at the tip of the long diagonal that is farthest from the widest
    cross-section, advances in steps of `thickness` until that cross-section, restarts
    exactly there and keeps going to the other tip. Therefore the widest cross-section is
    always a boundary, and only the last slice of each of the two runs may be thinner.
    """
    coords = np.array(polygon.exterior.coords)[:-1]
    projections = coords @ u
    s_min, s_max = float(projections.min()), float(projections.max())

    # The width of a convex polygon along `u` is piecewise linear in `s` with breakpoints at
    # the vertex projections, so the widest cross-section is attained at one of them.
    s_mid = max(projections, key=lambda s: _cross_section_length(polygon, u, v, float(s)))
    s_mid = float(s_mid)

    far_tip, near_tip = (s_min, s_max) if (s_mid - s_min) >= (s_max - s_mid) else (s_max, s_min)
    direction = 1.0 if near_tip > far_tip else -1.0

    boundaries = [far_tip]
    for start, stop in ((far_tip, s_mid), (s_mid, near_tip)):
        position = start
        while (stop - position) * direction > 1e-12:
            position = position + direction * thickness
            if (position - stop) * direction > 0:
                position = stop
            boundaries.append(position)

    return boundaries


def slice_polygon(polygon: Polygon, thickness: float | None = None,
                  default_slice_count: int = 20) -> List[Polygon]:
    """Cuts a polygon into thin slices perpendicular to its long diagonal.

    Args:
        polygon: The polygon to slice.
        thickness: The slice thickness along the long diagonal. Defaults to the polygon's
            extent along the long diagonal divided by `default_slice_count`.
        default_slice_count: Used to derive a default thickness when none is given.

    Returns:
        The slices, ordered starting from the tip farthest from the widest cross-section.
        Slices that end up empty are kept as empty polygons so that indices stay aligned
        between congruent polygons.
    """
    start, end = find_long_diagonal(polygon)
    axis = np.array([end[0] - start[0], end[1] - start[1]], dtype=float)
    u = axis / np.linalg.norm(axis)
    v = np.array([-u[1], u[0]])

    coords = np.array(polygon.exterior.coords)
    extent = float((coords @ u).max() - (coords @ u).min())
    if thickness is None:
        thickness = extent / default_slice_count
    if thickness <= 0:
        raise ValueError("Slice thickness must be positive.")

    boundaries = _slice_boundaries(polygon, u, v, thickness)
    half_span = float(np.abs(coords @ v).max()) + 1.0

    slices: List[Polygon] = []
    for lower, upper in zip(boundaries[:-1], boundaries[1:]):
        strip = Polygon([lower * u - half_span * v,
                         upper * u - half_span * v,
                         upper * u + half_span * v,
                         lower * u + half_span * v])
        piece = polygon.intersection(strip)
        if piece.is_empty:
            slices.append(Polygon())
        elif isinstance(piece, Polygon):
            slices.append(piece)
        else:
            slices.append(Polygon(piece.convex_hull))

    return slices


def neighbour_slice_intersections(slices_A: List[Polygon], slices_B: List[Polygon],
                                  neighbourhood: Tuple[int, ...] = (-1, 0, 1)) -> List[BaseGeometry]:
    """Collects the pairwise slice intersections summed by `neighbour_slice_overlap`.

    For every slice index `i` of A, the intersection with slices `i + k` of B is collected
    for each `k` in `neighbourhood`. Out-of-range indices are skipped.

    Args:
        slices_A: The slices of the fixed polygon.
        slices_B: The slices of the rotated polygon.
        neighbourhood: The relative slice indices to pair up.

    Returns:
        The non-empty intersection geometries, one per overlapping slice pair. The same
        region can appear more than once when a slice of A overlaps several slices of B,
        which mirrors how the areas are summed.
    """
    intersections: List[BaseGeometry] = []
    for i, slice_A in enumerate(slices_A):
        if slice_A.is_empty:
            continue
        min_x, min_y, max_x, max_y = slice_A.bounds
        for k in neighbourhood:
            j = i + k
            if not 0 <= j < len(slices_B):
                continue
            slice_B = slices_B[j]
            if slice_B.is_empty:
                continue
            other_min_x, other_min_y, other_max_x, other_max_y = slice_B.bounds
            if other_min_x > max_x or other_max_x < min_x or other_min_y > max_y or other_max_y < min_y:
                continue
            piece = slice_A.intersection(slice_B)
            if not piece.is_empty:
                intersections.append(piece)

    return intersections


def neighbour_slice_overlap(slices_A: List[Polygon], slices_B: List[Polygon],
                            neighbourhood: Tuple[int, ...] = (-1, 0, 1)) -> float:
    """Sums the overlap of each slice of A with the neighbouring slices of B.

    For every slice index `i` of A, the overlap with slices `i + k` of B is accumulated for
    each `k` in `neighbourhood`. Out-of-range indices are skipped.

    Args:
        slices_A: The slices of the fixed polygon.
        slices_B: The slices of the rotated polygon.
        neighbourhood: The relative slice indices to pair up.

    Returns:
        The total overlap area.
    """
    return sum(piece.area for piece in neighbour_slice_intersections(slices_A, slices_B, neighbourhood))


def calculate_full_rotation_slice_overlap(polygon: Polygon, pivot: Point, offset: Point,
                                          thickness: float | None = None) -> dict[float, float]:
    """Calculates the neighbouring-slice overlap sum for a full rotation.

    The slices of the rotating polygon are cut once in its own frame and then transformed
    rigidly, so they stay perpendicular to its long diagonal at every angle.

    Args:
        polygon: The fixed polygon.
        pivot: The pivot point of the rotating polygon.
        offset: The translation applied to the rotating polygon.
        thickness: The slice thickness along the long diagonal.

    Returns:
        A mapping of rotation angle to the neighbouring-slice overlap sum.
    """
    slices_A = slice_polygon(polygon, thickness=thickness)
    base_slices_B = [translate(piece, offset.x, offset.y) for piece in slice_polygon(polygon, thickness=thickness)]

    overlaps: dict[float, float] = {}
    for angle in np.linspace(0.0, 2 * np.pi, 360):
        slices_B = [rotate(piece, angle, origin=pivot, use_radians=True) for piece in base_slices_B]
        overlaps[angle] = neighbour_slice_overlap(slices_A, slices_B)

    return overlaps


def generate_kite(half_angle: float, diagonal: float, side_length: float) -> Polygon:
    """
    Generates a kite-shaped polygon using the given half-angle (in degrees), diagonal length, and side length.

    :param half_angle: Half the angle (in degrees) between the diagonal and a side of the kite.
    :param diagonal: The length of the major diagonal.
    :param side_length: The length of each of the kite's sides.
    :return: A Shapely Polygon representing the kite.
    """
    # Convert angle to radians
    half_angle_rad = math.radians(half_angle)

    # Calculate the perpendicular distance from the diagonal's center to the side
    perpendicular_distance = side_length * math.sin(half_angle_rad)

    # Calculate the horizontal offset for the side vertices
    half_minor_diagonal = side_length * math.cos(half_angle_rad)

    # Define the points of the kite
    # (0, 0) is the center of the kite
    top = (0, diagonal / 2)
    bottom = (0, -diagonal / 2)
    left = (-half_minor_diagonal, perpendicular_distance - (diagonal / 2))
    right = (half_minor_diagonal, perpendicular_distance - (diagonal / 2))

    # Create the kite polygon
    kite = Polygon([top, right, bottom, left, top])
    return kite
