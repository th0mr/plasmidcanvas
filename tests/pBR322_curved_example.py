"""
A copy of example 3 from the usage guide
This file can be ran individually or as part of the test
"""

from plasmidcanvas.plasmid import Plasmid
from plasmidcanvas.feature import ArrowFeature, RectangleFeature, RestrictionSite, MultiPairFeature

def get_curved_plasmid() -> Plasmid:

    # Define a plasmid of X base pairs long, with a name
    curved_plasmid = Plasmid("pBR322_3", 4361)
    curved_plasmid.set_marker_style("auto")
    curved_plasmid.set_feature_label_font_size(7)
    curved_plasmid.set_plasmid_line_width_sf(1.25)

    # Adding tcr
    tcr = ArrowFeature("tcr", 86,1276)
    curved_plasmid.add_feature(tcr)

    # Add rop protein for pBR322
    rop = ArrowFeature("rop", 1915,2106)
    rop.set_line_width_scale_factor(1.5)
    rop.set_color("purple")
    curved_plasmid.add_feature(rop)

    # Add a rectangle, base of mobility for pBR322
    bom = RectangleFeature("bom", 2208,2348)
    curved_plasmid.add_feature(bom)

    # Add ori
    ori = ArrowFeature("ori", 2534, 3122, -1)
    ori.set_color("orange")
    curved_plasmid.add_feature(ori)

    # # Add ampr
    ampr = ArrowFeature("ampr", 3293, 4153, -1)
    ampr.set_color("red")
    curved_plasmid.add_feature(ampr)

    for feature in curved_plasmid.get_features():
        if issubclass(feature.__class__, MultiPairFeature):
            feature.set_label_styles(["on-circle"]) 

    # # Add ampr promoter as an arrow
    ampr_promoter = ArrowFeature("ampr promoter", 4154, 4258, -1)
    ampr_promoter.set_color("darkred")

    ampr_promoter.set_line_width_scale_factor(0.75)
    curved_plasmid.add_feature(ampr_promoter)

    # Add the sites to the plasmid
    curved_plasmid.add_feature(RestrictionSite("BamHI", 375))
    curved_plasmid.add_feature(RestrictionSite("BfuAI - BspMI", 1054))
    curved_plasmid.add_feature(RestrictionSite("Bpu10I", 1581))
    curved_plasmid.add_feature(RestrictionSite("AflIII - PciI", 2473))
    curved_plasmid.add_feature(RestrictionSite("AhdI", 3366))

    return curved_plasmid

# Plot the plasmid
#get_curved_plasmid().save_to_file("curved_plasmid.png")