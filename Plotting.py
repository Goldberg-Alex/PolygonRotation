import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
from shapely.affinity import rotate, translate
from shapely.geometry.point import Point
from shapely.geometry.polygon import Polygon

from main import (calculate_full_rotation_overlap, calculate_full_rotation_slice_overlap,
                  neighbour_slice_intersections, slice_polygon)


def _draw_slices(ax, slices, color: str) -> None:
    """Draws the outline of every slice of a polygon."""
    for piece in slices:
        if piece.is_empty:
            continue
        x, y = piece.exterior.xy
        ax.plot(x, y, color=color, linewidth=0.5, alpha=0.4)


def _draw_slice_intersections(ax, intersections, color: str, label: str) -> None:
    """Fills the slice-pair intersections that make up the slice overlap sum."""
    labelled = False
    for piece in intersections:
        for part in getattr(piece, "geoms", [piece]):
            if part.is_empty or not isinstance(part, Polygon):
                continue
            x, y = part.exterior.xy
            ax.fill(x, y, color=color, alpha=0.3, linewidth=0,
                    label=None if labelled else label)
            labelled = True


def interactive_plot(polygon: Polygon, pivot: Point, offset: Point, slice_thickness: float | None = None):
    # The full-rotation curves and the slices do not depend on the slider, so they are
    # computed once instead of on every slider move.
    intersection_areas = calculate_full_rotation_overlap(polygon=polygon, pivot=pivot, offset=offset)
    slice_overlaps = calculate_full_rotation_slice_overlap(polygon=polygon, pivot=pivot, offset=offset,
                                                           thickness=slice_thickness)
    slices_A = slice_polygon(polygon, thickness=slice_thickness)
    base_slices_B = [translate(piece, offset.x, offset.y)
                     for piece in slice_polygon(polygon, thickness=slice_thickness)]

    def update(angle: float):

        # Update the top plot
        top_ax.clear()
        x1, y1 = polygon.exterior.xy
        top_ax.plot(x1, y1, label="Original Polygon")
        rotated_polygon = translate(polygon, offset.x, offset.y)
        x2, y2 = rotate(rotated_polygon, angle, origin=pivot, use_radians=True).exterior.xy
        top_ax.plot(x2, y2, label="Rotated Polygon")

        _draw_slices(top_ax, slices_A, color="tab:blue")
        rotated_slices_B = [rotate(piece, angle, origin=pivot, use_radians=True) for piece in base_slices_B]
        _draw_slices(top_ax, rotated_slices_B, color="tab:orange")
        _draw_slice_intersections(top_ax, neighbour_slice_intersections(slices_A, rotated_slices_B),
                                  color="tab:green", label="Slice overlap")

        top_ax.set_title("Polygons")
        top_ax.legend()
        top_ax.axis('equal')

        area_ax.clear()
        area_ax.plot(intersection_areas.keys(), intersection_areas.values())

        area_ax.plot([angle, angle], [0, max(intersection_areas.values())], color="black", linestyle="--")
        area_ax.set_title("Area per angle")

        slice_ax.clear()
        slice_ax.plot(slice_overlaps.keys(), slice_overlaps.values(), color="tab:green")
        slice_ax.plot([angle, angle], [0, max(slice_overlaps.values())], color="black", linestyle="--")
        slice_ax.set_title("Neighbouring slice overlap per angle")

        fig.canvas.draw_idle()

    # Create the figure and axes
    fig, (top_ax, area_ax, slice_ax) = plt.subplots(3, 1, figsize=(6, 10))

    # Create the slider
    ax_slider = plt.axes([0.01, 0.3, 0.03, 0.6], facecolor='lightgrey')  # x, y, width, height
    slider = Slider(ax_slider, 'Angle (rad)', 0, 2 * np.pi, valinit=0,orientation='vertical')
    fig.tight_layout(rect=(0, 0.2, 0.8, 1))  # Reserve space for the slider at the bottom

    # Connect the slider to the update function
    slider.on_changed(update)

    # Initial plot
    update(0)

    # plt.tight_layout()
    plt.get_current_fig_manager().window.state('zoomed')
    plt.show()


def plot_polygon_rotation(intersection_areas:dict[float,float]):
    plt.plot(intersection_areas.values(), intersection_areas.keys())
    plt.grid(True)
    plt.show()

def plot_polygon(polygon: Polygon):
    """Utility function to plot a Shapely polygon."""
    x, y = polygon.exterior.xy
    plt.figure()
    plt.fill(x, y, alpha=0.5, edgecolor='black', facecolor='blue')
    plt.title("Generated Convex Polygon")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.grid(True)
    plt.show()



