import copy
import random
from typing import List, Tuple, Any

import matplotlib.pyplot as plt
import numpy as np
from shapely.affinity import rotate, translate
from shapely.geometry import Polygon, Point


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
