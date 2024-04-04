from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.patches import FancyArrowPatch, Polygon, RegularPolygon, Wedge
import numpy as np

from plasmidcanvas._utils import circular_length, circular_midpoint, to_counter_clockwise
from plasmidcanvas._curvedtext import CurvedText
from plasmidcanvas._utils import DEFAULT_LABEL_FONT_SIZE

# ==================================
# Abstract Feature Classes


class Feature:
    """
    Abstract base class representing a feature
    """
    _DEFAULT_COLOR: str = "#069AF3"

    name: str
    color: str = _DEFAULT_COLOR

    def __init__(self, name: str) -> None:
        """
        Initializes the Feature object.

        Parameters
        ----------
        name : str
            The name of the feature.
        """
        self.set_name(name)

    def _render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        """
        Abstract method to render the feature on a matplotlib Axes object.

        Parameters
        ----------
        ax : Axes
            The matplotlib Axes object to render the feature on.
        p_total_base_pairs : int
            The total number of base pairs on the plasmid.
        p_center : tuple(float, float)
            The center coordinates of the plasmid.
        p_radius : float
            The radius length of the plasmid.
        p_line_width : float
            The line width of the plasmid circle.
        """
        pass

    def get_name(self) -> str:
        """
        Get the name of the feature.

        Returns
        -------
        str
            The name of the feature.
        """
        return self.name

    def set_name(self, name: str):
        """
        Set the name of the feature.

        Parameters
        ----------
        name : str
            The name to set for the feature.
        """
        self.name = name

    def get_color(self) -> str:
        """
        Get the color of the feature.

        Returns
        -------
        str
            The color of the feature.
        """
        return self.color

    def set_color(self, color: str) -> None:
        """
        Set the color of the feature.

        Parameters
        ----------
        color : str
            The color to set for the feature. Use words e.g "red" or hex values e.g. "#FFFFFF"
        """
        self.color = color


class MultiPairFeature(Feature):
    """
    An abstract Feauture type that spans multiple base pairs
    """
    _label_font_size: int = DEFAULT_LABEL_FONT_SIZE

    start_pair: int
    end_pair: int

    SUPPORTED_LABEL_STYLES = ["on-circle", "off-circle"]
    label_style = ["off-circle"]

    _orbit: int = 0

    def __init__(self, name: str, start_pair: int, end_pair: int) -> None:
        """
        Initializes the MultiPairFeature object.

        Parameters
        ----------
        name : str
            The name of the feature.
        start_pair : int
            The starting base pair of the feature.
        end_pair : int
            The ending base pair of the feature.
        """
        super().__init__(name)
        self.set_start_pairs(start_pair)
        self.set_end_pairs(end_pair)

    def length(self) -> int:
        """
        Calculates and returns the length of the feature.
        i.e. end_pair - start_pair

        Returns
        -------
        int
            The length of the feature.
        """
        return self.get_end_pair() - self.get_start_pair()

    def get_start_pair(self) -> int:
        """
        Get the starting base pair of the feature.

        Returns
        -------
        int
            The starting base pair.
        """
        return self.start_pair

    def set_start_pairs(self, start_pair: int) -> None:
        """
        Set the starting base pair of the feature.

        Parameters
        ----------
        start_pair
            The starting base pair to set.
        """
        self.start_pair = start_pair

    def get_end_pair(self) -> int:
        """
        Get the ending base pair of the feature.

        Returns
        -------
        int
            The ending base pair.
        """
        return self.end_pair

    def set_end_pairs(self, end_pair) -> None:
        """
        Set the ending base pair of the feature.

        Parameters
        ----------
        end_pair
            The ending base pair to set.
        """
        self.end_pair = end_pair

    def get_orbit(self) -> int:
        """
        Get the orbit of the feature.

        Returns
        -------
        int
            The orbit value of the feature.
        """
        return self._orbit

    def set_orbit(self, orbit: int) -> None:
        """
        Set the orbit of the feature. Mostly used internally for overlapping features, but can be used to force features into an orbit.
        The higher the orbit, the closer to the center of the plasmid the feature is placed. i.e orbit = 0 is the default, orbit = 4 circles four times the width of a feature closer to the center.

        Parameters
        ----------
        orbit : int
            The orbit value to set.
        """
        self._orbit = orbit

    def get_label_font_size(self) -> int:
        """
        Get the font size of the labels associated with the feature.

        Returns
        -------
        int
            The font size of the labels.
        """
        return self._label_font_size

    def set_label_font_size(self, font_size: int) -> None:
        """
        Set the font size of the labels associated with the feature.

        Parameters
        ----------
        font_size : int
            The font size to set.
        """
        self._label_font_size = font_size

    def get_label_styles(self) -> list[str]:
        """
        Get the list of label styles associated with the feature.

        Returns
        -------
        label_style : list[str]
            The list of label styles.
        """
        return self.label_style

    def set_label_styles(self, label_styles: list[str]):
        """
        Set the list of label styles associated with the feature.

        Parameters
        ----------
        label_styles : list[str]
            The list of label styles to set.
            Currently supported styles

                "on-circle" - The label is applied onto the feature using curved text. (The base pair location is not included in this labelling)
                
                "off-circle" - A label is placed off the circle with the name and bp range, connected by a small line

        Note
        ----
        A plasmid can have both, one or none of the label styles applied at the same time

        Raises
        ------
        ValueError
            If the given label_style is not in MultiPairFeature.SUPPORTED_LABEL_STYLES as listed above.
        """
        for style in label_styles:
            if style not in self.SUPPORTED_LABEL_STYLES:
                raise ValueError(
                    f"The label style '{style}' is not supported. Only the following are supported {self.SUPPORTED_LABEL_STYLES}")
        # If all styles are valid, replace the style array with the given one
        self.label_style = label_styles

    def _get_feature_labels(self, p_total_base_pairs: int, start_pair: int = None, end_pair: int = None,
                            styles: list[str] = None) -> None:
        """
        Get the labels associated with the feature based off the value of label_style

        Parameters
        ----------
        p_total_base_pairs : int
            The total number of base pairs on the plasmid.
        start_pair : int, optional
            The starting base pair.
        end_pair : int, optional
            The ending base pair.
        styles : list[str], optional
            The list of label styles. Used to force the creation of a specific type of label if needed.

        Returns
        -------
        list[Label]
            The list of labels associated with the feature.
        """

        # Passing in values for start_pair and end_pair that are not None allows us to
        # create features wherever we need, otherwise the default start_pair and end_pair will be used
        start_pair = self.get_start_pair() if start_pair is None else start_pair
        end_pair = self.get_end_pair() if end_pair is None else end_pair

        # Allows us to force styles if needed
        styles = list(dict.fromkeys(self.label_style)) if styles is None else styles

        labels = []
        for style in styles:
            if style not in self.SUPPORTED_LABEL_STYLES:
                raise ValueError(
                    f"{style} is not in the list of supported styles {self.SUPPORTED_LABEL_STYLES}")

            if style == "on-circle":
                # Adds a curved label
                label_text = f"{self.get_name()}"
                label = CurvedMultiPairLabel(label_text, start_pair, end_pair)
                label.set_orbit(self.get_orbit())
                label.set_font_size(self._label_font_size)
                labels.append(label)

            elif style == "off-circle":
                # Adds an off circle label
                # Calculate the midpoint, accounting for features that pass through 0 using circular_midpoint
                label_base_pair_location: int = round(
                    circular_midpoint(start_pair, end_pair, p_total_base_pairs))
                label_text = f"{self.get_name()} ({start_pair} - {end_pair})"
                label = SinglePairLabel(label_text, label_base_pair_location)
                label._orbit_ofset = self.get_orbit()
                label.set_line_color(self.color)
                label.set_font_size(self._label_font_size)
                label.set_font_color(self.color)
                labels.append(label)

        return labels


class SinglePairFeature(Feature):
    """
    An abstract Feature class representing a Feature that relates to only a single base pair
    """
    base_pair: int

    def __init__(self, name: str, base_pair: int) -> None:
        """
        Initializes the SinglePairFeature object.

        Parameters
        ----------
        name : str
            The name of the feature.
        base_pair : int
            The base pair associated with the feature.
        """
        super().__init__(name)
        self.set_base_pair(base_pair)

    def get_base_pair(self) -> int:
        """
        Get the base pair associated with the feature.

        Returns
        -------
        int
            The base pair.
        """
        return self.base_pair

    def set_base_pair(self, base_pair) -> None:
        """
        Set the base pair associated with the feature.

        Parameters
        ----------
        base_pair
            The base pair to set.
        """
        self.base_pair = base_pair


class LabelBase():
    """
    A base class for labels.
    """
    
    DEFAULT_FONT_SIZE: int = DEFAULT_LABEL_FONT_SIZE
    DEFAULT_FONT_COLOR: str = "black"

    label_text: str = "UntitledLabel"
    font_color: str = DEFAULT_FONT_COLOR
    font_size: int = DEFAULT_FONT_SIZE

    def __init__(self, name: str) -> None:
        """
        Initializes the LabelBase object.

        Parameters
        ----------
        name : str
            The name of the label.
        """
        self.set_label_text(name)

    def get_font_color(self) -> str:
        """
        Get the color of the label font.

        Returns
        -------
        str
            The color of the label font.
        """
        return self.font_color

    def set_font_color(self, font_color: str) -> None:
        """
        Set the color of the label font.

        Parameters
        ----------
        font_color : str
            The color of the label font.
        """
        self.font_color = font_color

    def set_label_text(self, label_text: str) -> None:
        """
        Set the text of the label.

        Parameters
        ----------
        label_text : str
            The text of the label.
        """
        self.label_text = label_text

    def get_label_text(self) -> str:
        """
        Get the text of the label.

        Returns
        -------
        label_text : str
            The text of the label.
        """
        return self.label_text

    def get_font_size(self) -> int:
        """
        Get the size of the label font.

        Returns
        -------
        font_size : int
            The size of the label font.
        """
        return self.font_size

    def set_font_size(self, font_size: int) -> None:
        """
        Set the size of the label font.

        Parameters
        ----------
        font_size : int
            The size of the label font.
        """
        if font_size > 0:
            self.font_size = font_size
        else:
            raise ValueError(f"Font size cannot be negative")
# ==================================
# Concrete Feature Classes


class SinglePairLabel(SinglePairFeature, LabelBase):
    """
    A class representing a label associated with a single base pair feature on a plasmid. Inheriting from both SinglePairLabel and LabelBase.
    """
    _DEFAULT_LINE_LENGTH_SF: float = 1
    _DEFAULT_LINE_ALPHA: float = 0.3
    _DEFAULT_LINE_COLOR: str = "black"

    line_length_sf: float = _DEFAULT_LINE_LENGTH_SF
    line_color: str = _DEFAULT_LINE_COLOR

    # Used internally for making label lines longer for internal orbits
    _orbit_ofset: int = 0

    def __init__(self, name: str, base_pair: int) -> None:
        """
        Initializes the SinglePairLabel object.

        Parameters
        ----------
        name : str
            The name of the label.
        base_pair : int
            The base pair associated with the label.
        """
        super().__init__(name=name, base_pair=base_pair)
        LabelBase.__init__(self, name)

    def _render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        """
        Renders the label on the provided Axes object.

        Parameters
        ----------
        ax : Axes
            The matplotlib Axes object to render the label on.
        p_total_base_pairs : int
            The total number of base pairs in the plasmid.
        p_center : tuple
            The coordinates of the center of the plasmid.
        p_radius : float
            The radius of the plasmid.
        p_line_width : float
            The width of the lines connecting the label to the features.
        """
        degrees = (self.get_base_pair() / p_total_base_pairs) * 360
        radians = np.deg2rad(degrees)

        # Calculate the x and y of the inner coordinates of the line
        start_xy = ((p_center[0] + p_radius) * np.sin(radians),
                    (p_center[1] + p_radius) * np.cos(radians))

        # To make the line for the label stick out farther we times by a scale factor to make the radius larger
        end_xy = (((p_center[0] + p_radius + (p_line_width * (2 + self._orbit_ofset))) * self.line_length_sf) * np.sin(
            radians),
                  ((p_center[1] + p_radius + (p_line_width * (2 + self._orbit_ofset)) * self.line_length_sf) * np.cos(
                      radians)))

        align = "left" if (degrees <= 180) else "right"

        # Plot the line and text for the label
        # Using zorder=3 to bring the line and label forward so they dont get covered up by other overlapping features

        ax.plot([start_xy[0], end_xy[0]], [start_xy[1], end_xy[1]], color=self.get_line_color(),
                alpha=self._DEFAULT_LINE_ALPHA, zorder=3)
        ax.text(x=end_xy[0], y=end_xy[1], s=self.get_label_text(), fontsize=self.get_font_size(),
                color=self.font_color, horizontalalignment=align, verticalalignment="center", zorder=3)

    def get_line_color(self) -> str:
        """
        Get the color of the line connecting the label to the feature.

        Returns
        -------
        str
            The color of the line.
        """
        return self.line_color

    def set_line_color(self, line_color: str) -> None:
        """
        Set the color of the line connecting the label to the feature.

        Parameters
        ----------
        line_color : str
            The color of the line.
        """
        self.line_color = line_color

    def get_line_length_sf(self) -> float:
        """
        Get the scale factor for the length of the line connecting the label to the feature.

        Returns
        -------
        float
            The scale factor for the line length.
        """
        return self.line_length_sf

    def set_line_length_sf(self, line_length_sf: float) -> None:
        """
        Set the scale factor for the length of the line connecting the label to the feature.

        Parameters
        ----------
        line_length_sf : float
            The scale factor for the line length.
        """
        self.line_length_sf = line_length_sf


# Currently just an alias for a SinglePairLabel
class RestrictionSite(SinglePairLabel):
    """
    A class representing a restriction site on a plasmid as an extention of SinglePairLabel.
    """
    def __init__(self, text: str, base_pair: int) -> None:
        """
        Initializes the RestrictionSite object. Creating a label in the format "<text> (<base_pair>)". 
        For more control over the label text, use SinglePairLabel.

        Parameters
        ----------
        text : str
            The text to display on the restiction site. e.g. "BspMI", "BfuAI, BspMI", "Multiple Cloning Site"
        base_pair : int
            The base pair position of the restriction site.
        """
        super().__init__(text, base_pair)
        self.set_label_text(f"{self.get_name()} ({self.get_base_pair()})")


class CurvedMultiPairLabel(MultiPairFeature, LabelBase):
    """
    A class representing a curved label associated with multiple base pairs on a plasmid.
    """

    SUPPORTED_CURVE_ALIGNMENTS: list[str] = ["bottom", "top"]
    # Set to either 'bottom' or 'top' to snap charecters to the bottom or top of the curve
    _alignment_to_curve: str = "bottom"

    def __init__(self, name: str, start_pair: int, end_pair: int) -> None:
        """
        Initializes the CurvedMultiPairLabel object.

        Parameters
        ----------
        name : str
            The name of the curved label.
        start_pair : int
            The starting base pair position of the label.
        end_pair : int
            The ending base pair position of the label.
        """
        super().__init__(name, start_pair, end_pair)
        LabelBase.__init__(self, name)

    def _render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        """
        Renders the curved label on the plot.

        Parameters
        ----------
        ax : Axes
            The matplotlib Axes object.
        p_total_base_pairs : int
            The total number of base pairs.
        p_center : tuple
            The coordinates of the center of the plot.
        p_radius : float
            The radius of the plot.
        p_line_width : float
            The width of the line.
        """

        start_radians = np.deg2rad((self.get_start_pair() / p_total_base_pairs) * 360)
        end_radians = np.deg2rad((self.get_end_pair() / p_total_base_pairs) * 360)

        label_center_base_pair = round(
            circular_midpoint(start=self.get_start_pair(), end=self.get_end_pair(), total_count=p_total_base_pairs))
        center_bp_angle = (label_center_base_pair / p_total_base_pairs) * 360

        # We do a slightly different approach depending on whether the center of the feature falls in the
        # top or the bottom of the circle to ensure the text is readable
        if 90 < center_bp_angle < 270:
            # We use 'top' text to curve aligment for these curved labels to place the letters on the top of the curve
            self.set_curve_alignment("top")
            # The curve must be "backwards" i.e. created from the end to the start of the feature
            curve = [np.sin(np.linspace(end_radians, start_radians)),
                     np.cos(np.linspace(end_radians, start_radians))]
        else:
            # The default style for the text to curve alignment is 'bottom'
            # meaning the letters are split up and placed under the curve, with the top of the letter closest to the curve
            # We make the curve go from start to end here to achieve the right results

            # We also account for curves that pass through zero by assisting in assembling a split curve and stitching it back together
            if end_radians < start_radians:
                curve_to_two_pi = [np.sin(np.linspace(start_radians, (2 * np.pi))),
                                   np.cos(np.linspace(start_radians, (2 * np.pi)))]
                curve_from_zero_to_end = [np.sin(np.linspace(0, end_radians)),
                                          np.cos(np.linspace(0, end_radians))]

                # TODO - COMBINE THESE TWO TO CREATE THE RIGHT EFFECT
                curve = [np.concatenate((curve_to_two_pi[0], curve_from_zero_to_end[0])),
                         np.concatenate((curve_to_two_pi[1], curve_from_zero_to_end[1]))]

            else:
                curve = [np.sin(np.linspace(start_radians, end_radians)),
                         np.cos(np.linspace(start_radians, end_radians))]

        # ============== Centering the text on the feature ===============

        # To center the text in the middle of the line we need to know the font size height
        # We need a renderer to do this.
        # Get a throwaway figure, were only using it to get a renderer object
        fig = plt.figure(0, figsize=(6, 6), dpi=300)
        renderer = fig.canvas.get_renderer()
        # Make some text, we wont be displaying this
        useless_text = plt.text(0, 0, 'useless text', fontsize=self.get_font_size())
        useless_text.set_alpha(0.0)
        # Get the bounding box for the rendered text
        bb = useless_text.get_window_extent(renderer=renderer)
        text_height = bb.height
        plt.close(0)

        ORBIT_GAP_SF = 1.25

        # Scale the curve out to the radius of the circle minus the line width, adjusted for text height
        scaled_curve = [[x * (p_center[0] + (p_radius - (p_line_width / 2) - (text_height))) for x in curve[0]],
                        [y * (p_center[1] + (p_radius - (p_line_width / 2) - (text_height))) for y in curve[1]]]

        # ================================================================

        CurvedText(
            x=scaled_curve[0],
            y=scaled_curve[1],
            text=self.get_label_text(),
            va=self._alignment_to_curve,
            axes=ax,  # This calls ax.add_artist, no need to add it
            fontsize=self.get_font_size()
        )

    def set_curve_alignment(self, curve_align: str) -> None:
        """
        Sets the alignment of the curved label. i.e. does the text follow the bottom or top of the curve?
        (use bottom for curves in the bottom half of the circle and top for the top half)

        Parameters
        ----------
        curve_align : str
            The alignment of the curved label ('bottom' or 'top').
        """
        self._alignment_to_curve = curve_align

    def get_curve_alignment(self) -> str:
        """
        Returns the alignment of the curved label.

        Returns
        -------
        str
            The alignment of the curved label.
        """
        return self._alignment_to_curve

class RectangleFeature(MultiPairFeature):
    """
    A class representing a curved rectangular feature associated with multiple base pairs on a plasmid.
    """
    # Center and radius for plotting the curved rectangle against plasmid circle
    radius: float

    line_width_scale_factor: float = 1

    def _render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        """
        Renders the rectangular feature on the plot.

        Parameters
        ----------
        ax : Axes
            The matplotlib Axes object.
        p_total_base_pairs : int
            The total number of base pairs.
        p_center : tuple
            The coordinates of the center of the plot.
        p_radius : float
            The radius of the plot.
        p_line_width : float
            The width of the line.
        """

        # Creating a rectangle
        # =============================================

        # 1 - Work out angles from 12 o clock
        # theta1 is the angle from the start of the plasmid (0th bp) to the start of the plasmid
        start_angle = (self.get_start_pair() / p_total_base_pairs) * 360
        # theta2 is the angle from the start of the feature to the end of the feature
        start_to_end_angle = start_angle + ((self.length() / p_total_base_pairs) * 360)

        # 2 - Convert angles starting from 12oclock clockwise into ones starting from 3oclock counterclockwise
        # Span from end to start to reverse the normal counterclockwise behaviour of matplotlib when
        theta_1 = to_counter_clockwise(start_to_end_angle)
        theta_2 = to_counter_clockwise(start_angle)

        # TODO - Support on / off the circle placement by moving radius

        # 3 - Place the wedge

        # Amount to stretch the arrow outwards / inwards to achieve the line width scale effect
        adjust_amount = ((p_line_width * self.get_line_width_scale_factor()) - p_line_width) / 2
        line_width_adjusted_radius = p_radius + adjust_amount
        rectangle = Wedge(center=p_center, r=line_width_adjusted_radius, theta1=theta_1, theta2=theta_2,
                          width=(p_line_width * self.get_line_width_scale_factor()), color=self.color)
        ax.add_patch(rectangle)

        # Labelling
        # ==================================================

        if "none" not in list(dict.fromkeys(self.label_style)):
            labels = self._get_feature_labels(p_total_base_pairs)
            for label in labels:
                label._render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)

    def get_line_width_scale_factor(self) -> float:
        """
        Returns the scale factor for the width of the rectangle.

        Returns
        -------
        float
            The scale factor for the width of the rectangle.
        """
        return self.line_width_scale_factor

    def set_line_width_scale_factor(self, sf: float) -> None:
        """
        Sets the scale factor for the width of the rectangle.

        Parameters
        ----------
        sf : float
            The scale factor for the width of the rectangle.
        """
        self.line_width_scale_factor = sf

class DirectionalMultiPairFeature(MultiPairFeature):
    """
    An abstract class representing a directional feature associated with multiple base pairs on a plasmid.
    """
    direction: int = 1

    def __init__(self, name: str, start_pair: int, end_pair: int, direction: int = 1):
        """
        Initializes a directional multi-pair feature.

        Parameters
        ----------
        name : str
            The name of the feature.
        start_pair : int
            The start base pair.
        end_pair : int
            The end base pair.
        direction : int, optional
            The direction of the feature, either 1 for clockwise or -1 for anti-clockwise (default is 1).
        """
        super().__init__(name, start_pair, end_pair)
        self.set_direction(direction)

    def get_direction(self) -> int:
        """
        Returns the direction of the feature.

        Returns
        -------
        int
            The direction of the feature. (1 for clockwise or -1 for anti-clockwise)
        """
        return self.direction

    def set_direction(self, direction: int) -> None:
        """
        Sets the direction of the feature.

        Parameters
        ----------
        direction : int
            The direction of the feature, either 1 for clockwise or -1 for anti-clockwise.
        
        Raises
        ------
        ValueError
            If direction is neither 1 or -1.
        """
        if direction not in [1, -1]:
            raise ValueError(
                f"{direction} is an invalid direction value. Direction can only be 1 (clockwise) or -1 (anti-clockwise)")
        self.direction = direction

class ArrowFeature(DirectionalMultiPairFeature):
    """
    A class representing a curved arrow feature associated with multiple base pairs on a plasmid.
    """
    line_width_scale_factor: float = 1

    def _render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        """
        Renders the arrow feature on the provided axes.

        Parameters
        ----------
        ax : Axes
            The axes object on which to render the arrow.
        p_total_base_pairs : int
            The total number of base pairs on the plasmid.
        p_center : tuple
            The coordinates of the center of the plasmid.
        p_radius : float
            The radius of the plasmid.
        p_line_width : float
            The width of the line representing the plasmid.
        """

        # NOTE - To clarify the mentions of "start" and "end" of an arrow
        # This side of the arrow is the "end" <| This side is the "start"

        # ======== Determine the length of the arrow head ======

        # The arrow head takes up 50% of the feature, up to a maximum of 1.5% total size of the plasmid
        # This avoids distortion if the arrows are too long, but allows smaller features to still have
        # a visable arrow
        MAX_PERCENT_AS_ARROW_HEAD = 0.50
        ARROW_LENGTH_CAP = p_total_base_pairs * 0.015
        max_arrow_length = circular_length(self.get_start_pair(), self.get_end_pair(), p_total_base_pairs) * MAX_PERCENT_AS_ARROW_HEAD
        length_of_arrow_head = max_arrow_length if max_arrow_length < ARROW_LENGTH_CAP else ARROW_LENGTH_CAP

        # ======= Determining the arrow start and end ========

        # If direction is clockwise (1) then set the arrowhead to be made on the end of the feature
        if self.get_direction() == 1:
            arrow_start_base_pair = int((self.get_end_pair() - length_of_arrow_head) % p_total_base_pairs)
            arrow_end_base_pair = self.get_end_pair()
        # If direction is counter-clockwise (-1) then set the arrowhead to be made at the start of the feature
        elif self.get_direction() == -1:
            arrow_start_base_pair = int((self.get_start_pair() + length_of_arrow_head) % p_total_base_pairs)
            arrow_end_base_pair = self.get_start_pair()
        else:
            raise ValueError(f"direction={self.get_direction()} is invalid")

        # Can be a float value
        angle_of_arrow_start = (arrow_start_base_pair / p_total_base_pairs) * 360
        angle_of_arrow_start_radians = np.deg2rad(angle_of_arrow_start)
        angle_of_end_of_arrow = (arrow_end_base_pair / p_total_base_pairs) * 360
        angle_of_arrow_end_radians = np.deg2rad(angle_of_end_of_arrow)

        # ======== Work out three points to assemble the traingle =======

        triangle_front_xy = ((p_center[0] + p_radius - (p_line_width / 2)) * np.sin(angle_of_arrow_end_radians),
                             (p_center[1] + p_radius - (p_line_width / 2)) * np.cos(angle_of_arrow_end_radians))

        # Amount to stretch the arrow outwards / inwards to achieve the line width scale effect
        adjust_amount = ((p_line_width * self.get_line_width_scale_factor()) - p_line_width) / 2

        triangle_point_inner_circle_xy = (((p_center[0] + p_radius - p_line_width - adjust_amount) * np.sin(
            angle_of_arrow_start_radians)),
                                           ((p_center[1] + p_radius - p_line_width - adjust_amount) * np.cos(
                                               angle_of_arrow_start_radians)))

        triangle_point_outer_circle_xy = (((p_center[0] + p_radius + adjust_amount) * np.sin(
            angle_of_arrow_start_radians)),
                                           ((p_center[1] + p_radius + adjust_amount) * np.cos(
                                               angle_of_arrow_start_radians)))

        triangle_points = [triangle_front_xy, triangle_point_inner_circle_xy, triangle_point_outer_circle_xy]

        # ======== Create the triangle and rectangle to make a curved arrow =======

        ax.add_patch(                    #Create triangle as arrow head
            Polygon(
                xy=triangle_points,            # (x,y) values
                color=self.get_color()
            )
        )

        if self.get_direction() == 1:
            rectangle_start_pair = self.get_start_pair()
            rectangle_end_pair = arrow_start_base_pair

        elif self.get_direction() == -1:
            rectangle_start_pair = arrow_start_base_pair
            rectangle_end_pair = self.get_end_pair()

        rectangle = RectangleFeature(self.name, rectangle_start_pair, rectangle_end_pair)

        rectangle.set_line_width_scale_factor(self.line_width_scale_factor)
        rectangle.set_color(self.get_color())

        # Prevent labelling inside the rectangles render, as rectangle_end_pair is not
        # the actual end_pair of the whole feature, so prevent labelling and control it below
        rectangle.label_style = ["none"]

        # Label the arrow feature
        if "none" not in list(dict.fromkeys(self.label_style)):
            labels = self._get_feature_labels(p_total_base_pairs, self.get_start_pair(), self.get_end_pair())   
            for label in labels:
                    label._render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)

        rectangle._render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)

    def get_line_width_scale_factor(self) -> float:
        """
        Returns the scale factor for the width of the arrow.

        Returns
        -------
        float
            The scale factor for the width of the arrow.
        """

        return self.line_width_scale_factor

    def set_line_width_scale_factor(self, sf: float) -> None:
        """
        Sets the scale factor for the width of the arrow.

        Parameters
        ----------
        sf : float
            The scale factor to set.
        """
        self.line_width_scale_factor = sf
