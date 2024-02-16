from matplotlib.axes import Axes
from matplotlib.patches import FancyArrowPatch, Polygon, RegularPolygon, Wedge
import numpy as np

from plasmidcanvas.curvedtext import CurvedText

# ==================================
# Abstract Feature Classes

# Takes an angle read clockwise from 12oclock, where the plasmid logically starts
# and converts it to an angle that matplotlib wedges can use, i.e. counterclockwise
# starting from 3
def to_counter_clockwise(clockwise_angle):
    counterclockwise_angle = (90 - clockwise_angle) % 360
    return counterclockwise_angle

class Feature:

    name: str

    def __init__(self, name: str) -> None:
        self.set_name(name)

    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        pass

    def get_name(self) -> str:
        return self.name
    
    def set_name(self, name: str):
        self.name = name
    

class MultiPairFeature(Feature):

    start_pair: int
    end_pair: int

    def __init__(self, name: str, start_pair: int, end_pair: int)  -> None:
        super().__init__(name)
        self.set_start_pairs(start_pair)
        self.set_end_pairs(end_pair)

    def length(self) -> int:
        return self.get_end_pairs() - self.get_start_pairs()

    def get_start_pairs(self) -> int:
        return self.start_pair

    def set_start_pairs(self, start_pair) -> None:
        self.start_pair = start_pair

    def get_end_pairs(self) -> int:
        return self.end_pair

    def set_end_pairs(self, end_pair) -> None:
        self.end_pair = end_pair

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

class SinglePairLabel(SinglePairFeature):
    
    DEFAULT_LINE_LENGTH_SF: float = 1.2
    DEFAULT_FONT_SIZE: int = 7
    DEFAULT_LINE_ALPHA: float = 0.2
    DEFAULT_LINE_COLOR: str = "black"

    label_text: str = "UntitledLabel"
    line_color: str = DEFAULT_LINE_COLOR
    
    def __init__(self, name: str, base_pair: int) -> None:
        super().__init__(name=name, base_pair=base_pair)
        self.set_label_text(self.get_name())
        
    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        
        degrees = (self.get_base_pair() / p_total_base_pairs) * 360
        radians = np.deg2rad(degrees)
        
        # Calculate the x and y of the inner coordinates of the line
        start_xy = ((p_center[0] + p_radius) * np.sin(radians),
                    (p_center[1] + p_radius) * np.cos(radians))

        print(start_xy)

        # To make the line for the label stick out farther we times by a scale factor to make the radius larger
        end_xy = (((p_center[0] + p_radius + p_line_width) * self.DEFAULT_LINE_LENGTH_SF) * np.sin(radians),
                  ((p_center[1] + p_radius + p_line_width) * self.DEFAULT_LINE_LENGTH_SF) * np.cos(radians))

        print(end_xy)

        print(f"degrees{degrees}")
        align= "left" if (degrees <= 180) else "right"
        print(f"align={align}")

        ax.plot([start_xy[0], end_xy[0]], [start_xy[1], end_xy[1]], color=self.get_line_color(), alpha=self.DEFAULT_LINE_ALPHA)
        ax.text(x=end_xy[0], y=end_xy[1], s=self.get_label_text(), fontsize=self.DEFAULT_FONT_SIZE, horizontalalignment=align, verticalalignment="center")

    def set_label_text(self, label_text:str) -> None:
        self.label_text = label_text

    def get_label_text(self) -> str:
        return self.label_text

    # TODO - Should color be of colour type or just str?
    def get_line_color(self) -> str:
        return self.line_color

    def set_line_color(self, line_color: str) -> None:
        self.color = line_color

# Currently just an alias for a SinglePairLabel
class RestrictionSite(SinglePairLabel):
    
    def __init__(self, name:str , base_pair: int) -> None:
        super().__init__(name, base_pair)
        self.set_label_text(f"{self.get_name()} ({self.get_base_pair()})")


class RectangleFeature(MultiPairFeature):

    SUPPORTED_LABEL_STYLES = ["on-circle", "off-circle", "inside-circle"]
    label_style = ["on-circle"]

    # Center and radius for plotting the curved rectangle against plasmid circle
    #center
    radius: float
    
    # TODO - Look at whether I should just be doing deg2rad instead of the moronic clockwise analogies

    line_width_scale_factor: float = 1

    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        
        # Creating a rectangle
        # =============================================
        
        # 1 - Work out angles from 12 o clock
        # theta1 is the angle from the start of the plasmid (0th bp) to the start of the plasmid
        twelve_to_start_angle = (self.get_start_pairs() / p_total_base_pairs) * 360
        print(f"twelve_to_start_angle={twelve_to_start_angle}")
        # theta2 is the angle from the start of the feature to the end of the feature
        start_to_end_angle = twelve_to_start_angle + ((self.length() / p_total_base_pairs) * 360)
        print(f"start_to_end_angle={start_to_end_angle}")

        # 2 - Convert angles starting from 12oclock clockwise into ones starting from 3oclock counterclockwise
        # Span from end to start to reverse the normal counterclockwise behaviour of matplotlib when
        theta_1=to_counter_clockwise(start_to_end_angle)
        print(f"theta_1={theta_1}")
        theta_2=to_counter_clockwise(twelve_to_start_angle)
        print(f"theta_2={theta_2}")

        # 3 - TODO - Adjust placement of the rectangle
        # For the rectangle to sit on top of the plasmid ring, the radius must be beyond 
        # the radius of the plasmid + half the plasmid circle width
        r_radius: float = p_radius
        
        # rwidth is the width of the rectangle, calculated by adding a scale factor to the width of the plasmid line
        r_width: float = p_radius

        # TODO - Support on / off the circle placement by moving radius

        # 4 - Place the wedge
        rectangle = Wedge(center=p_center, r=r_radius, theta1=theta_1, theta2=theta_2, width=(p_line_width * self.line_width_scale_factor), label="blah testing")
        ax.add_patch(rectangle)

        # Labelling 
        # ==================================================
        
        # Add all relevant label styles, removing duplicates
        for style in  list(dict.fromkeys(self.label_style)):
            if style not in self.SUPPORTED_LABEL_STYLES:
                raise ValueError(f"{style} is not in the list of supported styles {self.SUPPORTED_LABEL_STYLES}")
            
            if style == "on-circle":
                # Text must curve
                curved_text = CurvedText()
            
            elif style == "off-circle":
                # Add a label
                label_base_pair_location: int = round((self.get_start_pairs() + self.get_end_pairs()) / 2)
                label_text = f"{self.get_name()} ({self.get_start_pairs()} - {self.get_end_pairs()})"
                label = SinglePairLabel(label_text, label_base_pair_location)
                # Set line color to the same as the feature colour
                # TODO - MAKE THIS CHANGE
                label.set_line_color("blue")
                label.render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)
                
        
    
    def set_line_width_scale_factor(self, sf: float) -> None:
        self.line_width_scale_factor = sf
        

class ArrowFeature(RectangleFeature):

    # TODO - Add direction
    direction: int = 1

    def render(self, ax: Axes, p_total_base_pairs: int, p_center, p_radius: float, p_line_width: float) -> None:
        super().render(ax, p_total_base_pairs, p_center, p_radius, p_line_width)
        
        # The angle from which the arrow will center at
        angle = (self.get_end_pairs() / p_total_base_pairs) * 360
        print(f"name={self.name}, angle={angle}")
        
        radians = np.deg2rad(angle)
        
        # the angle from the start of the plasmid (0th bp) to the start of the plasmid
        twelve_to_start_angle = (self.get_start_pairs() / p_total_base_pairs) * 360
        
        #  the angle from the start of the feature to the end of the feature
        start_to_end_angle = twelve_to_start_angle + ((self.length() / p_total_base_pairs) * 360)
        
        angle_radians_theta_1 = np.radians(twelve_to_start_angle)
        angle_radians_theta_2 = np.radians(start_to_end_angle)
        
        #========Line
        ## TRY CHANGE ANGLE TO 90
        # arc = Arc(p_center,p_radius,p_radius,angle=angle,
        #     theta1=0,theta2=angle_radians_theta_2,capstyle='round',linestyle='-',lw=10,color="black")
        # ax.add_patch(arc)


        #========Create the arrow head
        
        endX=p_center[0]+(p_radius-(p_line_width/2))*np.sin(radians) #Do trig to determine end position
        endY=p_center[1]+(p_radius-(p_line_width/2))*np.cos(radians)

        print(f"endX={endX}")
        print(f"endY={endY}")

        print(f"rotated at {np.radians(angle+angle_radians_theta_2)}")

        ax.add_patch(                    #Create triangle as arrow head
            RegularPolygon(
                (endX, endY),            # (x,y)
                3,                       # number of vertices
                p_line_width/2,                # radius
                np.radians(angle),     # orientation
                color="black"
            )
        )
        #ax.set_xlim([p_center[0]-p_radius,p_center[1]+p_radius]) and ax.set_ylim([p_center[0]-p_radius,p_center[1]+p_radius]) 
        # Make sure you keep the axes scaled or else arrow will distort


class CircleFeature(MultiPairFeature):

    def __init__(self, start_pair: int, end_pair: int) -> None:
        super().__init__(start_pair, end_pair)
