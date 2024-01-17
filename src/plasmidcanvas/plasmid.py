
from typing import MutableSequence
import matplotlib.pyplot as plt
from math import pi
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Wedge

from feature import Feature, RectangleFeature, ArrowFeature

class Plasmid:

    DEFAULT_CIRCLE_RADIUS: float = 1000
    DEFAULT_CIRCLE_CENTER: tuple[float, float] = (0,0)
    DEFAULT_PLASMID_LINE_WIDTH: float = DEFAULT_CIRCLE_RADIUS * 0.05
    DEFAULT_PLASMID_NAME: str = "Untitled Plasmid"

    name: str = DEFAULT_PLASMID_NAME
    base_pairs: int = 0

    features: MutableSequence[Feature] = []

    radius: float

    def __init__(self, name: str, base_pairs: int) -> None:
        self.set_base_pairs(base_pairs)
        self.set_name(name)

    def plot(self) -> None:
        fig: Figure
        ax: Axes
        # Create plot
        fig, ax = plt.subplots()

        # Set x,y scaling to be equal
        ax.set_aspect('equal')

        # Place the plasmid circle onto the figure
        self.render(ax)

        ax.set_xlim((-self.DEFAULT_CIRCLE_RADIUS * 1.5, self.DEFAULT_CIRCLE_RADIUS * 1.5))
        ax.set_ylim((-self.DEFAULT_CIRCLE_RADIUS * 1.5, self.DEFAULT_CIRCLE_RADIUS * 1.5))

        

        # Add all features to the plasmid map by running their render() method
        for feature in self.get_features():
            feature.render(ax, self.get_base_pairs(), (0,0), self.DEFAULT_CIRCLE_RADIUS, self.DEFAULT_PLASMID_LINE_WIDTH)

    def save_to_file(self, fig: Figure, filename: str) -> None:
        pass


    def render(self, ax: Axes) -> None:
        center: tuple[float, float] = self.DEFAULT_CIRCLE_CENTER
        
        plasmid_circle=Wedge(center=center, r=self.DEFAULT_CIRCLE_RADIUS, theta1=0, theta2=360, width=self.DEFAULT_PLASMID_LINE_WIDTH)

        # TODO - Read this from PlasmidStyle

        # TODO - Figure out what colours a user can provide
        plasmid_circle.set_color("orange")

        # Make a horizontally aligned label in the center of the circle
        # Create bold plasmid name line and non-bold base pair line, using \bf math text
        plasmid_label_text = r"$\bf{" + self.name + "}$\n" + f"{self.base_pairs}bp"
        ax.annotate(plasmid_label_text, xy=center, ha="center", va="center")
        ax.add_patch(plasmid_circle)


    def add_feature(self, feature: Feature) -> None:
        # TODO - Add checks to ensure a multi-span feature cannot be created that is "out of bounds"
        self.features.append(feature)

    def add_restriction_site(self) -> None:
        pass

    # ====== Getters + Setters =======
    
    def get_features(self) -> MutableSequence[Feature]:
        return self.features;

    def get_base_pairs(self) -> int:
        return self.base_pairs

    def set_base_pairs(self, base_pairs) -> None:
        self.base_pairs = base_pairs

    def get_name(self) -> str:
        return self.name

    def set_name(self, name: str) -> None:
        self.name = name

class PlasmidStyle:

    base_pair_label_interval: int
    #default_rectangle_style: FeatureStyle
    #default_circle_style: FeatureStyle
    #default_arrow_style: FeatureStyle
    #default_label_style: LabelStyle

    def __init__(self) -> None:
        pass


plasmid = Plasmid("MyPlasmid", 500)

first_rect = RectangleFeature(100,200)
first_arrow = ArrowFeature(300,350)
second_arrow = ArrowFeature(400,500)

plasmid.add_feature(first_rect)
plasmid.add_feature(first_arrow)
plasmid.add_feature(second_arrow)

plasmid.plot()
