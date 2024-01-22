import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Arrow, Circle, FancyArrowPatch, PathPatch
from matplotlib.transforms import Affine2D
import numpy as np

def draw_curved_arrow(ax, center, radius, start_hour, end_hour, width=0.2, color='black'):
    arrow_length = 0.5
    arrow_width = 0.1

    # Convert hours to degrees
    start_angle = start_hour * 360 / 12
    end_angle = end_hour * 360 / 12

    # Draw the arc
    arc = Arc(xy=center, width=2 * radius, height=2 * radius,
              theta1=start_angle, theta2=end_angle, color=color)
    ax.add_patch(arc)

    # Calculate arrowhead position
    arrow_x = center[0] + radius * np.cos(np.radians(end_angle))
    arrow_y = center[1] + radius * np.sin(np.radians(end_angle))

    # Draw the arrowhead
    arrow_head = FancyArrowPatch((center[0], center[1]), (arrow_x, arrow_y),
                                 mutation_scale=15, color=color, arrowstyle='->',
                                 linewidth=arrow_width)
    ax.add_patch(arrow_head)

def plot_plasmid_map(features):
    fig, ax = plt.subplots(figsize=(8, 8))

    # Set the aspect ratio to be equal
    ax.set_aspect('equal', adjustable='box')

    # Plot a circle representing the plasmid
    circle = plt.Circle((0, 0), 1, edgecolor='black', facecolor='none')
    ax.add_patch(circle)

    # Plot features as curved arrows
    for feature in features:
        start_hour = feature['start_hour']
        end_hour = feature['end_hour']
        color = feature.get('color', 'black')

        draw_curved_arrow(ax, center=(0, 0), radius=1, start_hour=start_hour, end_hour=end_hour, color=color)

    # Set axis limits and remove ticks
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.show()

# Example usage
features = [
    {'start_hour': 2, 'end_hour': 8, 'color': 'red'},
    {'start_hour': 6, 'end_hour': 9, 'color': 'blue'},
    {'start_hour': 10, 'end_hour': 2, 'color': 'green'},
]

plot_plasmid_map(features)