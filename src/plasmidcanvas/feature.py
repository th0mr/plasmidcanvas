from matplotlib.axes import Axes
from matplotlib.patches import FancyArrowPatch, Polygon, RegularPolygon, Wedge
import numpy as np

from curvedtext import CurvedText

#from plasmidcanvas.curvedtext import CurvedText

# ==================================
# Abstract Feature Classes

# Takes an angle read clockwise from 12oclock, where the plasmid logically starts
# and converts it to an angle that matplotlib wedges can use, i.e. counterclockwise
# starting from 3
def to_counter_clockwise(clockwise_angle):
    counterclockwise_angle = (90 - clockwise_angle) % 360
    return counterclockwise_angle

class Feature:

    DEFAULT_COLOR: str = "#069AF3"

    name: str
    color: str = DEFAULT_COLOR

    def __init__(self, name: str) -> None:
        self.set_name(name)

    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        pass

    def get_name(self) -> str:
        return self.name
    
    def set_name(self, name: str):
        self.name = name

    def get_color(self) -> str:
        return self.color

    def set_color(self, color: str) -> None:
        self.color = color
    

class MultiPairFeature(Feature):

    start_pair: int
    end_pair: int

    SUPPORTED_LABEL_STYLES = ["on-circle", "off-circle", "inside-circle"]
    label_style = ["off-circle"]

    _orbit: int = 0

    def __init__(self, name: str, start_pair: int, end_pair: int)  -> None:
        super().__init__(name)
        self.set_start_pairs(start_pair)
        self.set_end_pairs(end_pair)

    def length(self) -> int:
        return self.get_end_pair() - self.get_start_pair()

    def get_start_pair(self) -> int:
        return self.start_pair

    def set_start_pairs(self, start_pair) -> None:
        self.start_pair = start_pair

    def get_end_pair(self) -> int:
        return self.end_pair

    def set_end_pairs(self, end_pair) -> None:
        self.end_pair = end_pair

    def get_orbit(self) -> int:
        return self._orbit
    
    def set_orbit(self, orbit: int) -> None:
        self._orbit = orbit

    # TODO - Rethink how labelling works
    def _get_feature_labels(self, start_pair=None, end_pair=None) -> None:

        # Passing in values for start_pair and end_pair that are not None allows us to
        # create features wherever we need, otherwise the default start_pair and end_pair will be used
        start_pair = self.get_start_pair() if start_pair is None else start_pair
        end_pair = self.get_end_pair() if end_pair is None else end_pair

        labels = []
        for style in  list(dict.fromkeys(self.label_style)):
            if style not in self.SUPPORTED_LABEL_STYLES:
                raise ValueError(f"{style} is not in the list of supported styles {self.SUPPORTED_LABEL_STYLES}")
            
            if style == "on-circle":
                # Text must curve
                curved_text = "TODO DO SOMETHING"
            
            elif style == "off-circle":
                # Add a label
                
                label_base_pair_location: int = round((start_pair + end_pair) / 2)
                label_text = f"{self.get_name()} ({start_pair} - {end_pair})"
                print(label_text)
                label = SinglePairLabel(label_text, label_base_pair_location)
                label._orbit_ofset = self.get_orbit()
                # TODO - MAKING LABEL LINE LENGTH SF TIMES BY 0.5 MAKES THEM GO IN
                #label.line_length_sf=label.line_length_sf * 0.5
                # Set line color to the same as the feature colour
                # TODO - MAKE THIS CHANGE
                label.set_line_color(self.color)
                label.set_font_color(self.color)
                labels.append(label)

        return labels
                

class SinglePairFeature(Feature):

    base_pair: int

    def __init__(self, name: str, base_pair: int) -> None:
        super().__init__(name)
        self.set_base_pair(base_pair)
        
    def get_base_pair(self) -> int:
        return self.base_pair

    def set_base_pair(self, base_pair) -> None:
        self.base_pair = base_pair

# ==================================
# Concrete Feature Classes

class LabelBase():
    DEFAULT_FONT_SIZE: int = 7
    DEFAULT_FONT_COLOR: str = "black"

    label_text: str = "UntitledLabel"
    font_color: str = DEFAULT_FONT_COLOR
    font_size: int = DEFAULT_FONT_SIZE

    def __init__(self, name: str) -> None:
        self.set_label_text(self.get_name())

    def get_font_color(self) -> str:
        return self.font_color

    def set_font_color(self, font_color: str) -> None:
        self.font_color = font_color

    def set_label_text(self, label_text:str) -> None:
        self.label_text = label_text

    def get_label_text(self) -> str:
        return self.label_text
    
    def get_font_size(self) -> int:
        return self.font_size
    
    def set_font_size(self, font_size: int) -> None:
        if font_size > 0:
            self.font_size = font_size
        else:
            raise ValueError(f"Font size cannot be negative")


class SinglePairLabel(SinglePairFeature, LabelBase):
    
    DEFAULT_LINE_LENGTH_SF: float = 1
    DEFAULT_LINE_ALPHA: float = 0.3
    DEFAULT_LINE_COLOR: str = "black"

    line_length_sf: float = DEFAULT_LINE_LENGTH_SF
    line_color: str = DEFAULT_LINE_COLOR

    # Used internally for making label lines longer for internal orbits
    _orbit_ofset: int = 0
    
    def __init__(self, name: str, base_pair: int) -> None:
        super().__init__(name=name, base_pair=base_pair)
        LabelBase.__init__(self, name)
        
    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        
        degrees = (self.get_base_pair() / p_total_base_pairs) * 360
        radians = np.deg2rad(degrees)
        
        # Calculate the x and y of the inner coordinates of the line
        start_xy = ((p_center[0] + p_radius) * np.sin(radians),
                    (p_center[1] + p_radius) * np.cos(radians))


        # To make the line for the label stick out farther we times by a scale factor to make the radius larger
        end_xy = (((p_center[0] + p_radius + (p_line_width*(2 + self._orbit_ofset))) * self.line_length_sf) * np.sin(radians),
                  ((p_center[1] + p_radius + (p_line_width*(2 + self._orbit_ofset)) * self.line_length_sf) * np.cos(radians)))

        align = "left" if (degrees <= 180) else "right"

        # Plot the line and text for the label
        # Using zorder=2 to bring the line and label forward so they dont get covered up by other overlapping features

        ax.plot([start_xy[0], end_xy[0]], [start_xy[1], end_xy[1]], color=self.get_line_color(), alpha=self.DEFAULT_LINE_ALPHA, zorder=3)
        ax.text(x=end_xy[0], y=end_xy[1], s=self.get_label_text(), fontsize=self.DEFAULT_FONT_SIZE,
                color=self.font_color, horizontalalignment=align, verticalalignment="center", zorder=3)

    def get_line_color(self) -> str:
        return self.line_color

    def set_line_color(self, line_color: str) -> None:
        self.line_color = line_color

    def get_line_length_sf(self) -> float:
        return self.line_length_sf
    
    def set_line_length_sf(self, line_length_sf: float) -> None:
        self.line_length_sf = line_length_sf


# Currently just an alias for a SinglePairLabel
class RestrictionSite(SinglePairLabel):
    
    def __init__(self, name:str , base_pair: int) -> None:
        super().__init__(name, base_pair)
        self.set_label_text(f"{self.get_name()} ({self.get_base_pair()})")

class CurvedMultiPairLabel(MultiPairFeature, LabelBase):

    def __init__(self, name: str, start_pair: int, end_pair: int) -> None:
        super().__init__(name, start_pair, end_pair)
        LabelBase.__init__(self, name)

    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:

        start_radians = np.deg2rad((self.get_start_pair() / p_total_base_pairs) * 360)
        end_radians   = np.deg2rad((self.get_end_pair() / p_total_base_pairs) * 360)

        print(f"start_radians = {start_radians}")
        print(f"end_radians = {end_radians}")

        # SWAP SIN AND COS TO MAKE IT WORK FOR TOP RIGHT QUADRANT

        curve = [np.cos(np.linspace(start_radians, end_radians)),
                 np.sin(np.linspace(start_radians, end_radians))]

        scaled_curve = [[x * (p_center[0] + p_radius) for x in curve[0]],
                        [y * (p_center[1] + p_radius) for y in curve[1]]]

        ax.plot(scaled_curve, color='b')

        CurvedText(
            x = scaled_curve[0],
            y = scaled_curve[1],
            text = self.get_label_text(),
            va = 'bottom',
            axes = ax,   # This calls ax.add_artist, no need to add it
            fontsize = self.get_font_size()
        )

class RectangleFeature(MultiPairFeature):

    # Center and radius for plotting the curved rectangle against plasmid circle
    #center
    radius: float

    line_width_scale_factor: float = 1

    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        
        print(f"endpair = {self.get_end_pair()}")
        # Creating a rectangle
        # =============================================
        
        # 1 - Work out angles from 12 o clock
        # theta1 is the angle from the start of the plasmid (0th bp) to the start of the plasmid
        start_angle = (self.get_start_pair() / p_total_base_pairs) * 360
        # theta2 is the angle from the start of the feature to the end of the feature
        start_to_end_angle = start_angle + ((self.length() / p_total_base_pairs) * 360)

        # 2 - Convert angles starting from 12oclock clockwise into ones starting from 3oclock counterclockwise
        # Span from end to start to reverse the normal counterclockwise behaviour of matplotlib when
        theta_1=to_counter_clockwise(start_to_end_angle)
        theta_2=to_counter_clockwise(start_angle)

        # 3 - TODO - Adjust placement of the rectangle
        # For the rectangle to sit on top of the plasmid ring, the radius must be beyond 
        # the radius of the plasmid + half the plasmid circle width
        r_radius: float = p_radius

        # TODO - Support on / off the circle placement by moving radius

        # 4 - Place the wedge
        rectangle = Wedge(center=p_center, r=r_radius, theta1=theta_1, theta2=theta_2, width=(p_line_width * self.line_width_scale_factor), label="blah testing", color=self.color)
        ax.add_patch(rectangle)

        # Labelling 
        # ==================================================
        
        if "none" not in list(dict.fromkeys(self.label_style)):
            labels = self._get_feature_labels()   
            for label in labels:
                label.render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)     
    
    def set_line_width_scale_factor(self, sf: float) -> None:
        self.line_width_scale_factor = sf

class DirectionalMultiPairFeature(MultiPairFeature):

    direction: int = 1

    def __init__(self, name: str, start_pair: int, end_pair: int, direction: int = 1):
        super().__init__(name, start_pair, end_pair)
        self.set_direction(direction)

    def get_direction(self) -> int:
        return self.direction
    
    def set_direction(self, direction: int) -> None:
        if direction not in [1, -1]:
            raise ValueError(f"{direction} is an invalid direction value. Direction can only be 1 (clockwise) or -1 (anti-clockwise)")
        self.direction = direction

class ArrowFeature(DirectionalMultiPairFeature):    

    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
              
        print(f"pre-rectangle endpair = {self.get_end_pair()}")
        
        # NOTE - To clarify the mentions of "start" and "end" of an arrow
        # This side of the arrow is the "end" <| This side is the "start"

        #======== Determine the length of the arrow head ======
        
        # The arrow head takes up 50% of the feature, up to a maximum of 1.5% total size of the plasmid
        # This avoids distortion if the arrows are too long, but allows smaller features to still have
        # a visable arrow
        MAX_PERCENT_AS_ARROW_HEAD = 0.50
        ARROW_LENGTH_CAP = p_total_base_pairs * 0.015
        max_arrow_length = self.length() * MAX_PERCENT_AS_ARROW_HEAD
        length_of_arrow_head = max_arrow_length if max_arrow_length < ARROW_LENGTH_CAP else ARROW_LENGTH_CAP        


        #======= Determining the arrow start and end ========

        # If direction is clockwise (1) then set the arrowhead to be made on the end of the feature
        if self.get_direction() == 1:
            arrow_start_base_pair = self.get_end_pair() - length_of_arrow_head
            arrow_end_base_pair = self.get_end_pair()
        # If direction is counter-clockwise (-1) then set the arrowhead to be made at the start of the feature
        elif self.get_direction() == -1:
            arrow_start_base_pair = self.get_start_pair() + length_of_arrow_head
            arrow_end_base_pair = self.get_start_pair()
        else:
            raise ValueError(f"direction={self.get_direction()} is invalid")
    
        # Can be a float value
        angle_of_arrow_start = (arrow_start_base_pair / p_total_base_pairs) * 360
        angle_of_arrow_start_radians = np.deg2rad(angle_of_arrow_start)
        angle_of_end_of_arrow = (arrow_end_base_pair / p_total_base_pairs) * 360
        angle_of_arrow_end_radians = np.deg2rad(angle_of_end_of_arrow)

        #======== Work out three points to assemble the traingle =======

        triangle_front_xy = ((p_center[0]+p_radius-(p_line_width/2))*np.sin(angle_of_arrow_end_radians),
                             (p_center[1]+p_radius-(p_line_width/2))*np.cos(angle_of_arrow_end_radians))
        
        triangle_point_inner_circle_xy = ((p_center[0]+p_radius-p_line_width)*np.sin(angle_of_arrow_start_radians),
                                          (p_center[1]+p_radius-p_line_width)*np.cos(angle_of_arrow_start_radians))
        
        triangle_point_outer_circle_xy = ((p_center[0]+p_radius)*np.sin(angle_of_arrow_start_radians),
                                          (p_center[1]+p_radius)*np.cos(angle_of_arrow_start_radians))

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
            rectangle_end_pair = int(arrow_start_base_pair)
        elif self.get_direction() == -1:
            rectangle_start_pair = int (arrow_start_base_pair)
            rectangle_end_pair = self.get_end_pair()
            
        rectangle = RectangleFeature(self.name, rectangle_start_pair, rectangle_end_pair)
        rectangle.set_color(self.get_color())
        # Prevent labelling inside the rectangles render, as rectangle_end_pair is not
        # the actual end_pair of the whole feature, so prevent labelling and control it below
        rectangle.label_style = ["none"]
        
        # Label the arrow feature
        if "none" not in list(dict.fromkeys(self.label_style)):
            labels = self._get_feature_labels()   
            for label in labels:
                label.render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)

        rectangle.render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)


class CircleFeature(MultiPairFeature):

    def __init__(self, start_pair: int, end_pair: int) -> None:
        super().__init__(start_pair, end_pair)
