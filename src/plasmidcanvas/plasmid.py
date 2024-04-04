from typing import MutableSequence
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Wedge
import numpy as np

from plasmidcanvas.feature import Feature, LabelBase, MultiPairFeature, RectangleFeature, ArrowFeature, RestrictionSite, SinglePairFeature, SinglePairLabel
from plasmidcanvas._utils import DEFAULT_LABEL_FONT_SIZE


class Plasmid:
    """
    Circular object representing a plasmid object on which to plot Features onto. 
    """
    _DEFAULT_CIRCLE_RADIUS: float = 1000
    _DEFAULT_CIRCLE_CENTER: tuple[float, float] = (0, 0)
    _DEFAULT_PLASMID_LINE_WIDTH: float = _DEFAULT_CIRCLE_RADIUS * 0.10
    _DEFAULT_PLASMID_LINE_WIDTH_SF: float = 1
    _DEFAULT_PLASMID_NAME: str = "Untitled Plasmid"
    _DEFAULT_PLASMID_COLOR: str = "grey"

    SUPPORTED_MARKER_STYLES = ["auto", "n_markers", "none"]
    _DEFAULT_MARKER_STYLE = "auto"
    _DEFAULT_MARKER_DISTANCE_SF = 1.10
    # Used for n_markers style
    _DEFAULT_NUMBER_OF_MARKERS: int = 16

    SUPPORTED_TICK_STYLES = ["auto", "none"]
    _DEFAULT_TICK_STYLE = "auto"

    name: str = _DEFAULT_PLASMID_NAME
    base_pairs: int
    _center: tuple[float, float] = _DEFAULT_CIRCLE_CENTER
    _radius: float = _DEFAULT_CIRCLE_RADIUS
    _features: MutableSequence[Feature] = []

    _plasmid_line_width: float = _DEFAULT_PLASMID_LINE_WIDTH
    _plasmid_line_width_sf: float = _DEFAULT_PLASMID_LINE_WIDTH_SF
    _marker_style: str = _DEFAULT_MARKER_STYLE
    _marker_distance_sf: float = _DEFAULT_MARKER_DISTANCE_SF
    _number_of_markers: int = _DEFAULT_NUMBER_OF_MARKERS
    _tick_style: str = _DEFAULT_TICK_STYLE
    _plasmid_color: str = _DEFAULT_PLASMID_COLOR
    _tick_color: str = _DEFAULT_PLASMID_COLOR
    _feature_label_font_size: int = DEFAULT_LABEL_FONT_SIZE

    def __init__(self, name: str, base_pairs: int) -> None:
        """
        Creates a Plasmid object with a given name and number of base pairs

        Parameters
        ----------
        name : str
            The name of the plasmid e.g. pBR322
        base_pairs : int
            The number of base pairs in the plasmid

        Raises
        ------
        ValueError
            If the number of base pairs is negative

        Examples
        --------
        Creating a Plasmid::

            from plasmidcanvas.plasmid import Plasmid
            myplasmid = Plasmid("pBR322", 4361)
        """
        if base_pairs <= 0:
            raise ValueError("The number of base pairs for a Plasmid must be greater than zero")
        self.set_base_pairs(int(base_pairs))
        self.set_name(name)
        self._features = []

    def _degrees_to_basepair(self, degree: float) -> int:
        """Takes a degree value, returns a base pair that is closest to that degree"""        
        return round((degree / 360) * self.base_pairs)

    def _basepair_to_degrees(self, basepair: int) -> float:
        """Takes a base pair value, returns the degree on the plasmid circle that it is located at"""    
        return (basepair / self.base_pairs) * 360

    def plot(self) -> Figure:
        """
        Plots all features added to the Plasmid object onto a matplotlib Figure.

        Note
        ----
        Unless you are working in an interactive environment, e.g. Jupyter, it is recommended 
        to use save_to_file() to view your plasmid instead.

        Returns
        -------
        figure : Figure
            A matplotlib figure object with the plasmid and its features plotted onto it

        Example
        -------
        Obtaining a plasmid Figure::

            from plasmidcanvas.plasmid import Plasmid
            myplasmid = Plasmid("pBR322", 4361)
            figure = myplasmid.plot()
            proceed to view or work with figure...
        """
    
        fig: Figure
        ax: Axes
        # Create plot
        fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

        # Set x,y scaling to be equal
        ax.set_aspect('equal')

        # Turn axis off
        plt.axis('off')

        # Place the plasmid circle onto the figure
        self._render(ax)

        XY_SCALING_FACTOR = 1.6

        ax.set_xlim((-self._radius * XY_SCALING_FACTOR, self._radius * XY_SCALING_FACTOR))
        ax.set_ylim((-self._radius * XY_SCALING_FACTOR, self._radius * XY_SCALING_FACTOR))

        # Place numbered tick markers around the plasmid to indicate increments of basepairs
        self._place_markers_at_degrees(ax, self._get_markers())
        self._place_ticks_at_degrees(ax, self._get_ticks())

        # ========== Placing features ============
        # Add all features to the plasmid map by running their _render() method

        placed_features = []
        # All features are kept in orbits, decrementing the radius that the features are placed
        # in based on the orbit if the feature is going to overlap. higher number = closer to the center
        orbit = 0

        multi_pair_features = [feature for feature in self.get_features() if issubclass(feature.__class__, MultiPairFeature)]

        # Sorting gives the proper cascading effect when placing multipair features in orbit
        # Also provides non-determinism, as otherwise the order in which features are added by the user
        # impacts the layout.
        sorted_multi_pair_features = sorted(multi_pair_features, key=lambda feature: feature.start_pair)
        non_multi_pair_features = [feature for feature in self.get_features() if not issubclass(feature.__class__, MultiPairFeature)]

        for feature in sorted_multi_pair_features:

            # A quick check to see if the feature's font size has changed from the default, if not, apply Plasmid
            # object's feature_label_font_size to the feature prior to rendering
            # Otherwise leave the font size as is, as it means the user has manually adjusted it

            if feature.get_label_font_size() == DEFAULT_LABEL_FONT_SIZE:
                feature.set_label_font_size(self.get_feature_label_font_size())

            pre_placement_orbit = orbit

            # Get all current placed multi pair features
            sorted_placed_multi_pair_features = sorted([feature for feature in placed_features if issubclass(feature.__class__, MultiPairFeature)],
                                                        key=lambda feature: feature.start_pair)
            for potential_overlap in sorted_placed_multi_pair_features:
                # Check if start_pair of feature lies between the start and end of the potential overlap feature on the same orbit
                if (potential_overlap.get_start_pair() <= feature.get_start_pair() <= potential_overlap.get_end_pair()
                    and potential_overlap.get_orbit() == orbit):
                    orbit += 1

            if pre_placement_orbit == orbit:
                # Reset orbit to 0
                orbit = 0

            ORBIT_GAP_SF = 1.25
            orbit_radius = self._radius - (self.get_plasmid_line_width() * ORBIT_GAP_SF * orbit)

            feature.set_orbit(orbit)

            feature._render(ax, self.get_base_pairs(), self._DEFAULT_CIRCLE_CENTER, orbit_radius, self.get_plasmid_line_width())

            placed_features.append(feature)

        for feature in non_multi_pair_features:
            feature._render(ax, self.get_base_pairs(), self._DEFAULT_CIRCLE_CENTER, self._radius, self.get_plasmid_line_width())
            placed_features.append(feature)

        fig.tight_layout()
        return fig

    def _get_markers(self) -> list[float]:
        """
        Returns a list of degrees to place base pair markers.

        Returns
        -------
        degrees_to_place_markers : list[float]
            A list of degrees to place markers at. Locations / number of markers is dependant on self._marker_style.
        """

        if self._marker_style == "auto":
            # Increment based on number of digits in the base pairs
            # e.g. if bp = 9    then increment is 10^(1-1) = 10^0 = 1
            #      if bp = 4050 then increment is 10^(4-1) = 10^3 = 1000
            marker_increment = 10**(len(str(self.base_pairs)) - 1)

            # If this marker increment leads to less or equal to 2 non-zero markers then half the marker increment 
            while (self.base_pairs / marker_increment <= 2):
                marker_increment = int(marker_increment / 2)

            marker_basepairs = list(range(0, self.base_pairs, max(1, marker_increment)))

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

        return []

    def _get_ticks(self) -> list[float]:
        """
        Returns a list of degrees to place base pair ticks.

        Returns
        -------
        degrees_to_place_markers : list[float]
            A list of degrees to place ticks at. Locations / number of markers is dependant on self_tick_style.
        """

        if self._tick_style == "auto":
            # Increment based on number of digits in the base pairs
            # e.g. if bp = 9    then increment is 10^(1-1) = 10^0 = 1
            #      if bp = 4050 then increment is 10^(4-1) = 10^3 = 1000
            marker_increment = 10**(len(str(self.base_pairs)) - 1) / 10
            # makes 4 digit bp plasmids have ticks of 100, 3 digit to 10 etc.

            marker_basepairs = list(range(0, self.base_pairs, max(1, round(marker_increment))))

            # Turn base pairs into degrees
            degrees_to_place_markers = map(self._basepair_to_degrees, marker_basepairs)
            return degrees_to_place_markers

        return []

    def _place_markers_at_degrees(self, ax, degrees_to_place_markers: list[float]) -> None:
        """
        Places base pair markers at a set of given degrees.

        Parameters
        -------
        ax : Axes
            A matplotlib Axes object on which to plot the marker text
        degrees_to_place_markers : list[float]
            A list of degrees to place markers at
        """

        for degree in degrees_to_place_markers:
            x, y = ((self._center[0] + self._radius * self._marker_distance_sf * np.sin(np.deg2rad(degree))),
                    (self._center[1] + self._radius * self._marker_distance_sf * np.cos(np.deg2rad(degree))))
            ax.text(x, y, s=f"{self._degrees_to_basepair(degree)}", horizontalalignment='center', fontstyle='italic', alpha=0.5, fontsize=7)

    def _place_ticks_at_degrees(self, ax: Axes, degrees_to_place_ticks: list[float]) -> None:
        """
        Places base pair ticks (small lines around circle) at a set of given degrees. 

        Parameters
        -------
        ax : Axes
            A matplotlib Axes object on which to plot the marker text
        degrees_to_place_markers : list[float]
            A list of degrees to place ticks at
        """

        for degree in degrees_to_place_ticks:

            radians = np.deg2rad(degree)

            start_xy = (((self._center[0] + self._radius) * np.sin(radians)),
                        ((self._center[1] + self._radius) * np.cos(radians)))

            TICK_DISTANCE_SF = 1.03

            end_xy = (((self._center[0] + self._radius) * np.sin(radians) * TICK_DISTANCE_SF),
                      ((self._center[1] + self._radius) * np.cos(radians) * TICK_DISTANCE_SF))

            ax.plot([start_xy[0], end_xy[0]], [start_xy[1], end_xy[1]], color=self._plasmid_color, alpha=0.3, zorder=0)

    def save_to_file(self, filename: str) -> None:
        """
        Plots the plasmid by calling Plasmid.plot() and saves the figure to an image with a given filename.

        Parameters
        ----------
        filename : str
            filename or path to save the figure to. e.g. 
                "myplasmid.png" would save it to the working directory. 
                "build/myplasmid.png" would save the same image inside the "build" folder.
            A file extention should be included in the filename. Any matplotlib supported file extention is supported.
            e.g. png, pdf, ps, eps and svg.

        Examples
        --------
        To save a plasmid to a png::

            from plasmidcanvas.plasmid import Plasmid
            myplasmid = Plasmid("pBR322", 4361)
            myplasmid.savefig("figure.png")
        """
        # TODO - Path validation
        fig = self.plot()
        fig.savefig(fname=filename, bbox_inches="tight")
        print(f"Plasmid map saved to {filename}")

    def _render(self, ax: Axes) -> None:
        """ 
        Renders a blank plasmid circle onto a given matplotlib Axes object.

        Parameters
        ----------
        ax : Axes
            matplotlib Axes object to plot the plasmid circle onto
        """
        center: tuple[float, float] = self._DEFAULT_CIRCLE_CENTER

        plasmid_circle = Wedge(center=center, r=self._DEFAULT_CIRCLE_RADIUS, theta1=0, theta2=360, width=self.get_plasmid_line_width())

        # TODO - Read this from PlasmidStyle

        # TODO - Figure out what colours a user can provide
        plasmid_circle.set_color(self._plasmid_color)
        plasmid_circle.set_alpha(0.5)

        # Make a horizontally aligned label in the center of the circle
        # Create bold plasmid name line and non-bold base pair line, using \bf math text
        plasmid_label_text = r"$\bf{" + self.name + "}$\n" + f"{self.base_pairs}bp"
        ax.annotate(plasmid_label_text, xy=center, ha="center", va="center")
        ax.add_patch(plasmid_circle)

    def add_feature(self, feature: Feature) -> None:
        """ 
        Adds a feature to the Plasmid object. All features to be plotted must be created and then added to the plasmid this way.
        See the example below. 
        
        Parameters
        ----------
        feature : Feature
            A Feature object to add to the plasmid. e.g RectangleFeature, ArrowFeature, SinglePairLabel, RestrictionSite etc.
            See the documentation for plasmidcanvas.feature to see the availble Feature types and their usage.

        Raises
        ------
        ValueError
            If the feature lies out of bounds of the Plasmid's base pair range
            
        Note
        ----
        This does not render the feature, that only happens when Plasmid.plot() or Plasmid.save_to_file() is ran.
        Therefore, plasmid wide and feature specific customisations can be made after the feature is added if you wish.

        Examples
        --------
        Adding an Arrow Featue::

            from plasmidcanvas.plasmid import Plasmid
            from plasmidcanvas.feature import ArrowFeature, RectangleFeature, RestrictionSite
            myplasmid = Plasmid("pBR322", 4361)
            myplasmid.savefig("figure.png")

            arrow = ArrowFeature(1000, 2000)
            # Optional feature customisation here
            myplasmid.add_feature(arrow)

            # Save plasmid out
            myplasmid.save_to_file("figure.png")
        """
        # Check if the feature to be added is out of bounds
        if issubclass(feature.__class__, MultiPairFeature):
            if feature.get_end_pair() > self.get_base_pairs():
                raise ValueError(f"""Feature {feature.get_name()} with end pair = {feature.get_end_pair()} 
                                   is out of bounds for a Plasmid of size {self.get_base_pairs()}""")
        if issubclass(feature.__class__, SinglePairFeature):
            if feature.get_base_pair() > self.get_base_pairs():
                raise ValueError(f"""Feature {feature.get_name()} with base pair = {feature.get_base_pair()} 
                                    is out of bounds for a Plasmid of size {self.get_base_pairs()}""")
        
        self._features.append(feature)

    # ====== Getters + Setters =======

    def get_features(self) -> MutableSequence[Feature]:
        """
        Get the features associated with the plasmid.

        Returns
        -------
        features : MutableSequence[Feature]
            The features associated with the plasmid.
        """
        return self._features

    def get_base_pairs(self) -> int:
        """
        Get the number of base pairs in the plasmid.

        Returns
        -------
        base_pairs : int
            The number of base pairs in the plasmid.
        """
        return self.base_pairs

    def set_base_pairs(self, base_pairs) -> None:
        """
        Set the number of base pairs in the plasmid.

        Parameters
        ----------
        base_pairs : int
            The number of base pairs to set the plasmid to.
        """
        self.base_pairs = base_pairs

    def get_name(self) -> str:
        """
        Get the name of the plasmid.

        Returns
        -------
        name : str
            The name of the plasmid.
        """
        return self.name

    def set_name(self, name: str) -> None:
        """
        Set the name of the plasmid.

        Parameters
        ----------
        name : str
            The name to set for the plasmid.
        """
        self.name = name

    def get_color(self) -> str:
        """
        Get the color of the plasmid circle.

        Returns
        -------
        color : str
            The color of the plasmid circle.
        """
        return self._plasmid_color

    def set_color(self, color: str) -> None:
        """
        Set the color of the plasmid circle.

        Parameters
        ----------
        color : str
            The color to set for the plasmid circle. Use words e.g "red" or hex values e.g. "#FFFFFF"
        """
        self._plasmid_color = color

    def get_center(self) -> tuple[float, float]:
        """
        Get the center coordinates of the plasmid.

        Returns
        -------
        center : tuple[float, float]
            The (x, y) coordinates of the center of the plasmid.
        """
        return self._center

    def set_center(self, center: tuple[float, float]) -> None:
        """
        Set the center coordinates of the plasmid.

        Parameters
        ----------
        center : tuple[float, float]
            The (x, y) coordinates to set as the center of the plasmid.
        """
        self._center = center

    def get_plasmid_line_width_sf(self) -> float:
        """
        Get the scale factor for the plasmid circle line width.

        Returns
        -------
        line_width_sf : float
            The scale factor for the plasmid circle line width.
        """
        return self._plasmid_line_width_sf

    def set_plasmid_line_width_sf(self, line_width_sf: float) -> None:
        """
        Set the scale factor for the plasmid circle line width.

        Parameters
        ----------
        line_width_sf : float
            The scale factor to set for the plasmid circle line width.
            e.g. 1.5 makes the plasmid circle 1.5 times as thick as the default
        """
        self._plasmid_line_width_sf = line_width_sf

    def get_plasmid_line_width(self) -> float:
        """
        Get the plasmid line width.

        Returns
        -------
        plasmid_line_width : float
            The plasmid line width.
        """
        return self._plasmid_line_width * self.get_plasmid_line_width_sf()

    def set_plasmid_line_width(self, plasmid_line_width: float) -> None:
        """
        Set the plasmid line width.

        Parameters
        ----------
        plasmid_line_width : float
            The plasmid line width to set.

        Note
        ----
        Either use set_plasmid_line_width or set_plasmid_line_width_sf. 
        Increasing both values could create a very thick plasmid circle!
        """
        self._plasmid_line_width = plasmid_line_width

    def get_marker_style(self) -> str:
        """
        Get the style of markers for the plasmid.

        Returns
        -------
        marker_style : str
            The style of markers used for the plasmid.
        """
        return self._marker_style

    def set_marker_style(self, marker_style: str) -> None:
        """
        Set the style of markers for the plasmid.
        
        Parameters
        ----------
        marker_style : str
            The style of markers to set for the plasmid.
            Currently supported styles

                "auto" (default) - Automatically apply markers at a reasonable marker interval based on the plasmid size

                "n_markers" - Place n equidistant markers around the circle. The default is 16, unless changed with Plasmid.set_number_of_markers()

                "none" - No markers added to the circle 

        Raises
        ------
        ValueError
            If the given marker_style is not in Plasmid.SUPPORTED_MARKER_STYLES, as listed above.

        """
        if marker_style not in self.SUPPORTED_MARKER_STYLES:
            raise ValueError(f"'{marker_style}' is not a supported marker style. The following are supported {self.SUPPORTED_MARKER_STYLES}")
        self._marker_style = marker_style

    def get_marker_distance_sf(self) -> float:
        """
        Get the scale factor for the distance between the circle and the marker text.

        Returns
        -------
        marker_distance_sf : float
            The scale factor for the marker distance.
        """
        return self._marker_distance_sf

    def set_marker_distance_sf(self, marker_distance_sf) -> None:
        """
        Set the scale factor for the distance between the circle and the marker text.

        Parameters
        ----------
        marker_distance_sf : float
            The scale factor to set for the marker distance.

        Examples
        --------
        This value is small by default (1.03). It is advised to alter this by applying a to the existing scale factor
        to achieve reliable increases.::

            myplasmid.set_marker_distance(myplasmid.get_marker_distance_sf() * 1.5)
        """
        self._marker_distance_sf = marker_distance_sf

    def get_number_of_markers(self) -> int:
        """
        Get the number of markers on the plasmid. Only used when the plasmid's marker style is "n_markers".

        Returns
        -------
        number_of_markers : int
            The number of markers on the plasmid.
        """
        return self._number_of_markers

    def set_number_of_markers(self, number_of_markers: int) -> None:
        """
        Set the number of markers to produce around the plasmid when the plasmid's marker style is "n_markers".

        Parameters
        ----------
        number_of_markers : int
            The number of markers to place around the plasmid when using n_marker style. 
        """
        self._number_of_markers = number_of_markers

    def get_tick_style(self) -> str:
        """
        Get the style of tick placement on the plasmid.

        Returns
        -------
        tick_style : str
            The style of ticks used for the plasmid.
        """
        return self._tick_style

    def set_tick_style(self, tick_style: str) -> None:
        """
        Set the style of tick placement on the plasmid.

        Parameters
        ----------
        tick_style : str
            The style of ticks to set for the plasmid.

            Currently supported:

                "auto" (default) - Automatically draw on ticks at a reasonable marker interval based on the plasmid size

                "none" - No ticks are drawn around the circle

        Raises
        ------
        ValueError
            If the given tick_style is not in Plasmid.SUPPORTED_TICK_STYLES, as listed above.

        """
        if tick_style not in self.SUPPORTED_TICK_STYLES:
            raise ValueError(f"'{tick_style}' is not a supported tick style. The following are supported {self.SUPPORTED_TICK_STYLES}")
        self._tick_style = tick_style

    def get_tick_color(self) -> str:
        """
        Get the color of ticks for the plasmid.

        Returns
        -------
        tick_color : str
            The color of ticks used for the plasmid.
        """
        return self._tick_color

    def set_tick_color(self, tick_color: str) -> None:
        """
        Set the color of ticks for the plasmid.

        Parameters
        ----------
        tick_color : str
            The color of ticks to set for the plasmid. Use words e.g "red" or hex values e.g. "#FFFFFF"
        """
        self._tick_color = tick_color

    def get_feature_label_font_size(self) -> int:
        """
        Get the override font size of feature labels for the plasmid.
        This font size will be applied to any labels / features that have not already had their font size changed manually.

        Returns
        -------
        feature_label_font_size : int
            The font size of feature labels used for the plasmid. Given as a pt value.
        """
        return self._feature_label_font_size

    def set_feature_label_font_size(self, feature_label_font_size: int) -> None:
        """
        Set the font size of feature labels for the plasmid.
        This font size will be applied to any labels / features that have not already had their font size changed manually.

        Parameters
        ----------
        feature_label_font_size : int
            The font size of feature labels to set for the plasmid. Given as a pt value.
        """
        self._feature_label_font_size = feature_label_font_size

