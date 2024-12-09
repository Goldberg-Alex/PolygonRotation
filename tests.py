import random
import unittest

from shapely.affinity import translate
from shapely.geometry import Polygon
from shapely.geometry.point import Point

from Plotting import interactive_plot
from main import generate_random_convex_polygon, is_polygon_convex, find_extremes, generate_kite


class TestGenerateRandomConvexPolygon(unittest.TestCase):

    def test_polygon_vertex_count(self):
        """Test if the generated polygon has the correct number of vertices."""
        num_vertices = 5
        bounds = (0, 10, 0, 10)
        polygon = generate_random_convex_polygon(num_vertices, bounds)

        self.assertEqual(len(polygon.exterior.coords) - 1, num_vertices)

    def test_polygon_with_minimum_vertices(self):
        """Test polygon generation with the minimum number of vertices."""
        num_vertices = 3
        bounds = (0, 10, 0, 10)
        polygon = generate_random_convex_polygon(num_vertices, bounds)
        self.assertEqual(len(polygon.exterior.coords) - 1, num_vertices)

    def test_polygon_is_valid(self):
        """Test if the generated polygon is valid."""
        num_vertices = 8
        bounds = (0, 10, 0, 10)
        polygon = generate_random_convex_polygon(num_vertices, bounds)
        self.assertTrue(polygon.is_valid)


class TestIsPolygonConvex(unittest.TestCase):
    def test_convex_polygon(self):
        """Test that a convex polygon is identified correctly."""
        convex_polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])  # Square (Convex)
        self.assertTrue(is_polygon_convex(convex_polygon), "Convex polygon incorrectly identified as concave.")

    def test_concave_polygon(self):
        """Test that a concave polygon is identified correctly."""
        concave_polygon = Polygon([(0, 0), (2, 0), (1, 1), (2, 2), (0, 2)])  # Star-like (Concave)
        self.assertFalse(is_polygon_convex(concave_polygon), "Concave polygon incorrectly identified as convex.")

    def test_triangle(self):
        """Test that a triangle is always convex."""
        triangle = Polygon([(0, 0), (2, 0), (1, 1)])  # Triangle (Convex)
        self.assertTrue(is_polygon_convex(triangle), "Triangle incorrectly identified as concave.")

    def test_line_polygon(self):
        """Test that a degenerate polygon (line) is not convex."""
        line_polygon = Polygon([(0, 0), (1, 1), (2, 2)])  # Degenerate polygon (Line)
        self.assertTrue(is_polygon_convex(line_polygon), "Degenerate line incorrectly identified as concave.")

    def test_collinear_polygon(self):
        """Test that a polygon with collinear points is considered convex."""
        collinear_polygon = Polygon([(0, 0), (2, 0), (4, 0)])  # Collinear points (Convex)
        self.assertTrue(is_polygon_convex(collinear_polygon), "Collinear polygon incorrectly identified as concave.")


class TestFullRotations(unittest.TestCase):
    def test_rotate_square(self):
        polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])  # Square (Convex)
        pivot = Point([0, 2])

        offset = Point([0, 0])
        self.assertTrue(find_extremes(polygon=polygon, offset=offset))
        interactive_plot(polygon, pivot, offset)

    def test_full_hypothesis(self):
        bound = 10
        random.seed(0)
        for _ in range(100):
            polygon = generate_random_convex_polygon(5, (-bound, bound, -bound, bound))
            offset = Point((random.randrange(-bound, bound), random.randrange(-bound, bound)))
            for pivot in translate(polygon, offset.x, offset.y).exterior.coords:
                extremes = find_extremes(polygon=polygon, offset=offset, pivot=Point(pivot))
                if len(extremes) > 2:
                    print(f"found example!"
                          f"\npolygon = {polygon.exterior.coords.xy},"
                          f"\noffset = {offset.coords.xy}"
                          f"\npivot = {pivot}"
                          f"{extremes=}")

                    interactive_plot(polygon=polygon, offset=offset, pivot=Point(pivot))

    def test_full_hypothesis_kite(self):
        """check the hypothesis specifically for kite polygons"""
        bound = 10
        random.seed(0)
        for _ in range(100):
            polygon = generate_kite(15, 8, side_length=3)
            offset = Point((random.randrange(-bound, bound), random.randrange(-bound, bound)))
            for pivot in translate(polygon, offset.x, offset.y).exterior.coords:
                extremes = find_extremes(polygon=polygon, offset=offset, pivot=Point(pivot))
                if len(extremes) > 2:
                    print(f"found example!"
                          f"\npolygon = {polygon.exterior.coords.xy},"
                          f"\noffset = {offset.coords.xy}"
                          f"\npivot = {pivot}"
                          f"{extremes=}")

                    interactive_plot(polygon=polygon, offset=offset, pivot=Point(pivot))

if __name__ == '__main__':
    unittest.main()
