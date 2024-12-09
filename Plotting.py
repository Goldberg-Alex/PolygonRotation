import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider
from shapely.affinity import rotate, translate
from shapely.geometry.point import Point
from shapely.geometry.polygon import Polygon

from main import calculate_full_rotation_overlap


def interactive_plot(polygon:Polygon, pivot:Point,offset:Point ):
    def update(angle:float):

        # Update the top plot
        top_ax.clear()
        x1, y1 = polygon.exterior.xy
        top_ax.plot(x1, y1, label="Original Polygon")
        rotated_polygon = translate(polygon, offset.x, offset.y)
        x2, y2 = rotate(rotated_polygon, angle, origin=pivot, use_radians=True).exterior.xy
        top_ax.plot(x2, y2, label="Rotated Polygon")
        top_ax.set_title("Polygons")
        top_ax.legend()
        top_ax.axis('equal')

        intersection_areas = calculate_full_rotation_overlap(polygon=polygon, pivot=pivot,offset=offset)
        bottom_ax.clear()
        bottom_ax.plot(intersection_areas.keys(), intersection_areas.values())

        bottom_ax.plot([angle, angle], [0, max(intersection_areas.values())], color="black", linestyle="--")
        bottom_ax.set_title("Area per angle")

        fig.canvas.draw_idle()

    # Create the figure and axes
    fig, (top_ax, bottom_ax) = plt.subplots(2, 1, figsize=(6, 8))

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



