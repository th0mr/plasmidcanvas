
from typing import MutableSequence
import matplotlib.pyplot as plt
from math import pi
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Wedge
import numpy as np

from feature import CurvedMultiPairLabel, Feature, MultiPairFeature, RectangleFeature, ArrowFeature, RestrictionSite, SinglePairLabel



class Plasmid:

    DEFAULT_CIRCLE_RADIUS: float = 1000
    DEFAULT_CIRCLE_CENTER: tuple[float, float] = (0,0)
    DEFAULT_PLASMID_LINE_WIDTH: float = DEFAULT_CIRCLE_RADIUS * 0.10
    DEFAULT_PLASMID_NAME: str = "Untitled Plasmid"
    
    SUPPORTED_MARKER_STYLES = ["auto", "n_markers", "none"]
    DEFAULT_MARKER_STYLE = "auto"
    DEFAULT_MARKER_DISTANCE_SF = 1.10
    # Used for n_markers style
    DEFAULT_NUMBER_OF_MARKERS: int = 16

    name: str = DEFAULT_PLASMID_NAME
    base_pairs: int
    _center: tuple[float, float] = DEFAULT_CIRCLE_CENTER
    _radius: float = DEFAULT_CIRCLE_RADIUS
    _features: MutableSequence[Feature] = []
    
    _plasmid_line_width: float = DEFAULT_PLASMID_LINE_WIDTH
    _marker_style: str = DEFAULT_MARKER_STYLE
    _marker_distance_sf: float = DEFAULT_MARKER_DISTANCE_SF
    _number_of_markers: int = DEFAULT_NUMBER_OF_MARKERS

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

        ax.set_xlim((-self._radius * XY_SCALING_FACTOR, self._radius * XY_SCALING_FACTOR))
        ax.set_ylim((-self._radius * XY_SCALING_FACTOR, self._radius * XY_SCALING_FACTOR))
        
        # Place numbered tick markers around the plasmid to indicate increments of basepairs
        self._place_markers_at_degrees(ax, self._get_markers())

        # Add all features to the plasmid map by running their render() method

        multi_pair_features = [feature for feature in self.get_features() if issubclass(feature.__class__, MultiPairFeature)]

        orbit = 0

        for feature in self.get_features():

            pre_placement_orbit = orbit

            # Deal with multi-pair feature overlaps
            if issubclass(feature.__class__, MultiPairFeature):
                for potential_overlap in multi_pair_features:
                    # Check if start_pair of feature lies between the start and end of the potential overlap feature
                    if potential_overlap.get_start_pair() < feature.get_start_pair() < potential_overlap.get_end_pair(): 
                        orbit += 1

            # If we have not incremented orbit during this feature placement, reset back to orbit 0
            if pre_placement_orbit == orbit:
                orbit = 0

            print(f"orbit for feature {feature.get_name()} = {orbit}")
            feature.render(ax, self.get_base_pairs(), self.DEFAULT_CIRCLE_CENTER, self._radius, self.DEFAULT_PLASMID_LINE_WIDTH)
            
            

    def _get_markers(self) -> list[float]:

        if self._marker_style == "auto":
            # Increment based on number of digits in the base pairs
            # e.g. if bp = 9    then increment is 10^(1-1) = 10^0 = 1
            #      if bp = 4050 then increment is 10^(4-1) = 10^3 = 1000
            marker_increment = 10**(len(str(self.base_pairs)) - 1)
            
            # If this marker increment leads to less or equal to 2 non-zero markers then half the marker increment 
            while (self.base_pairs / marker_increment <= 2):
                marker_increment = int(marker_increment / 2)

            marker_basepairs = list(range(0, self.base_pairs, max(1,marker_increment)))

            # Remove the last marker if it is going to be close to the inital zero marker
            # i.e. if the final marker is within 15% of a marker increment away from the final base pair
            if marker_basepairs[-1] >= self.base_pairs - (marker_increment * 0.15):
                marker_basepairs.pop(-1)

            # Turn base pairs into degrees
            degrees_to_place_markers = map(self._basepair_to_degrees, marker_basepairs)
            return degrees_to_place_markers
        
        elif self._marker_style == "n_markers":
            # Assemble a set of points to plot label markers around
            degrees_to_place_markers = np.linspace(0, 360, self._number_of_markers, endpoint=False)
            return degrees_to_place_markers
                
        elif self._marker_style == "none":
            return []

    def _place_markers_at_degrees(self, ax, degrees_to_place_markers: list[float]) -> None:
        for degree in degrees_to_place_markers:
            x,y = (self._center[0] + self._radius * self._marker_distance_sf * np.sin(np.deg2rad(degree)), 
            self._center[1] + self._radius * self._marker_distance_sf * np.cos(np.deg2rad(degree)))
            ax.text(x,y, s=f"{self._degrees_to_basepair(degree)}", horizontalalignment='center', fontstyle='italic', alpha=0.5, fontsize=7)


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
        self._features.append(feature)

    # ====== Getters + Setters =======
    
    def get_features(self) -> MutableSequence[Feature]:
        return self._features

    def get_base_pairs(self) -> int:
        return self.base_pairs

    def set_base_pairs(self, base_pairs) -> None:
        self.base_pairs = base_pairs

    def get_name(self) -> str:
        return self.name

    def set_name(self, name: str) -> None:
        self.name = name

    def get_center(self) -> tuple[float, float]:
        return self._center
    
    def set_center(self, center: tuple[float, float]) -> None:
        self._center = center

    def get_plasmid_line_width(self) -> float:
        return self._plasmid_line_width

    def set_plasmid_line_width(self, plasmid_line_width: float) -> None:
        self.plasmid_line_Width = plasmid_line_width

    def get_marker_style(self) -> str:
        return self._marker_style

    def set_marker_style(self, marker_style: str) -> None:
        if marker_style not in self.SUPPORTED_MARKER_STYLES:
            raise ValueError(f"'{marker_style}' is not a supported marker style. The following are supported {self.SUPPORTED_MARKER_STYLES}")
        self._marker_style = marker_style

    def get_marker_distance_sf(self) -> str:
        return self._marker_distance_sf
    
    def set_marker_distance_sf(self, marker_distance_sf) -> None:
        self._marker_distance_sf = marker_distance_sf

    def get_number_of_markers(self) -> int:
        return self._number_of_markers
    
    def set_number_of_markers(self, number_of_markers: int) -> None:
        self._number_of_markers = number_of_markers

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
plasmid.set_marker_style("auto")

# Adding an arrow
# for pBR322 this is TcR
tcr = ArrowFeature("TcR", 86,1276)
# # # Customise the thinkness of the line relative to the thickness of the plasmid circle
# # tcr.set_line_width_scale_factor(1.0)
plasmid.add_feature(tcr)

# # Add rop protein for pBR322
rop = ArrowFeature("rop", 1915,2106)
plasmid.add_feature(rop)

# # Add a rectangle, base of mobility for pBR322
bom = RectangleFeature("bom", 2208,2348)
plasmid.add_feature(bom)

# # Add ori
ori = ArrowFeature("ori", 2534, 3122, -1)
ori.color = "orange"
plasmid.add_feature(ori)

# # Add ampr - technically this arrow should have a portion segmented for its signal sequence
ampr = ArrowFeature("ampr", 3293, 4153, -1)
ampr.color = "red"
plasmid.add_feature(ampr)

overlapping = ArrowFeature("overlapping feature", 3500, 4300)
overlapping.color = "pink"
plasmid.add_feature(overlapping)

# # Add ampr promoter as an arrow
ampr_promoter = ArrowFeature("ampr promoter", 4154, 4258, -1)
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

# plasmid.add_feature(CurvedMultiPairLabel("TcR LOOOOOOOOOONG", 86, 1210))
# plasmid.add_feature(CurvedMultiPairLabel("rop", 1915, 2040))
# plasmid.add_feature(CurvedMultiPairLabel("bom", 1915, 2040))
# plasmid.add_feature(CurvedMultiPairLabel("ori", 2599, 3122))
# plasmid.add_feature(CurvedMultiPairLabel("ampr", 3358, 4153))
# plasmid.add_feature(CurvedMultiPairLabel("word", 4000, 4361))
# plasmid.add_feature(CurvedMultiPairLabel("word", 4300, 100))
# plasmid.add_feature(CurvedMultiPairLabel("ampr", 400, 3000))

# Plot the plasmid
plasmid.plot()
