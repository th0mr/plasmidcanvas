
from typing import MutableSequence
import matplotlib.pyplot as plt
from math import pi
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Wedge
import numpy as np

from feature import Feature, RectangleFeature, ArrowFeature, RestrictionSite, SinglePairLabel



class Plasmid:

    DEFAULT_CIRCLE_RADIUS: float = 1000
    DEFAULT_CIRCLE_CENTER = (0,0)
    DEFAULT_PLASMID_LINE_WIDTH: float = DEFAULT_CIRCLE_RADIUS * 0.10
    DEFAULT_PLASMID_NAME: str = "Untitled Plasmid"
    DEFAULT_PLASMID_NUMBER_OF_MARKERS: int = 16
    DEFAULT_MARKER_DISTANCE_SF = 1.10

    name: str = DEFAULT_PLASMID_NAME
    base_pairs: int = 0

    features: MutableSequence[Feature] = []

    radius: float

    def __init__(self, name: str, base_pairs: int) -> None:
        self.set_base_pairs(base_pairs)
        self.set_name(name)

    def _degrees_to_basepair(self, degree: float) -> int:
        return round((degree / 360) * self.base_pairs)

    def _basepair_to_degrees(self, basepair: int) -> float:
        return (basepair / self.base_pairs) * 360

    def plot(self) -> None:
        fig: Figure
        ax: Axes
        # Create plot
        fig, ax = plt.subplots()

        # Change figure height and width
        fig.set_figheight(6)
        fig.set_figwidth(6)
        # Set dpi to print quality
        fig.set_dpi(300)

        # Set x,y scaling to be equal
        ax.set_aspect('equal')

        # Turn axis off
        plt.axis('off')

        # Place the plasmid circle onto the figure
        self.render(ax)

        XY_SCALING_FACTOR = 1.6

        ax.set_xlim((-self.DEFAULT_CIRCLE_RADIUS * XY_SCALING_FACTOR, self.DEFAULT_CIRCLE_RADIUS * XY_SCALING_FACTOR))
        ax.set_ylim((-self.DEFAULT_CIRCLE_RADIUS * XY_SCALING_FACTOR, self.DEFAULT_CIRCLE_RADIUS * XY_SCALING_FACTOR))
        
        # Assemble a set of points to plot label markers around
        degrees_to_place_markers = np.linspace(0, 360, self.DEFAULT_PLASMID_NUMBER_OF_MARKERS, endpoint=False)
        for degree in degrees_to_place_markers:
            x,y = (self.DEFAULT_CIRCLE_CENTER[0] + self.DEFAULT_CIRCLE_RADIUS * self.DEFAULT_MARKER_DISTANCE_SF * np.sin(np.deg2rad(degree)), 
                   self.DEFAULT_CIRCLE_CENTER[1] + self.DEFAULT_CIRCLE_RADIUS * self.DEFAULT_MARKER_DISTANCE_SF * np.cos(np.deg2rad(degree)))
            ax.text(x,y, s=f"{self._degrees_to_basepair(degree)}", horizontalalignment='center', fontstyle='italic', alpha=0.5, fontsize=7)

        # Add all features to the plasmid map by running their render() method
        for feature in self.get_features():
            feature.render(ax, self.get_base_pairs(), self.DEFAULT_CIRCLE_CENTER, self.DEFAULT_CIRCLE_RADIUS, self.DEFAULT_PLASMID_LINE_WIDTH)

        

    def save_to_file(self, fig: Figure, filename: str) -> None:
        pass

    def render(self, ax: Axes) -> None:
        center: tuple[float, float] = self.DEFAULT_CIRCLE_CENTER
        
        plasmid_circle=Wedge(center=center, r=self.DEFAULT_CIRCLE_RADIUS, theta1=0, theta2=360, width=self.DEFAULT_PLASMID_LINE_WIDTH)

        # TODO - Read this from PlasmidStyle

        # TODO - Figure out what colours a user can provide
        plasmid_circle.set_color("grey")
        plasmid_circle.set_alpha(0.5)

        # Make a horizontally aligned label in the center of the circle
        # Create bold plasmid name line and non-bold base pair line, using \bf math text
        plasmid_label_text = r"$\bf{" + self.name + "}$\n" + f"{self.base_pairs}bp"
        ax.annotate(plasmid_label_text, xy=center, ha="center", va="center")
        ax.add_patch(plasmid_circle)


    def add_feature(self, feature: Feature) -> None:
        # TODO - Add checks to ensure a multi-span feature cannot be created that is "out of bounds"
        self.features.append(feature)

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


## TESTING
    

# Define a plasmid of X base pairs long, with a name
plasmid = Plasmid("pBR322", 4361)

# Adding an arrow
# for pBR322 this is TcR
tcr = ArrowFeature("TcR", 86,1276)
# # Customise the thinkness of the line relative to the thickness of the plasmid circle
# tcr.set_line_width_scale_factor(1.0)
plasmid.add_feature(tcr)

# Add rop protein for pBR322
rop = ArrowFeature("rop", 1915,2106)
plasmid.add_feature(rop)

# Add a rectangle, base of mobility for pBR322
bom = RectangleFeature("bom", 2208,2348)
plasmid.add_feature(bom)

# Add ori
# TODO - Change direction of this arrow to counter clockwise
# TODO - Should this be done by f.setDirection or by providing a counter clockwise order to the constructor
#        e.g. ArrowFeature(10, 1)
ori = ArrowFeature("ori", 2534, 3122)
ori.color = "orange"
plasmid.add_feature(ori)

# Add ampr - technically this arrow should have a portion segmented for its signal sequence
# TODO - This should be counter clockwise
ampr = ArrowFeature("ampr", 3293, 4153)
ampr.color = "red"
plasmid.add_feature(ampr)

# Add ampr promoter as an arrow
# TODO - This should be counter clockwise
ampr_promoter = ArrowFeature("ampr promoter", 4154, 4258)
ampr_promoter.color = "darkred"
plasmid.add_feature(ampr_promoter)

restriction_site_1 = RestrictionSite("BamHI", 375)
restriction_site_2 = RestrictionSite("BfuAI - BspMI", 1054)
restriction_site_3 = RestrictionSite("Bpu10I", 1581)
restriction_site_4 = RestrictionSite("AflIII - PciI", 2473)
restriction_site_5 = RestrictionSite("AhdI", 3366)

# Add the sites to the plasmid
plasmid.add_feature(restriction_site_1)
plasmid.add_feature(restriction_site_2)
plasmid.add_feature(restriction_site_3)
plasmid.add_feature(restriction_site_4)
plasmid.add_feature(restriction_site_5)

# Plot the plasmid
plasmid.plot()
