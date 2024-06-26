"""
A copy of example 2 from the usage guide

This file can be ran individually or as part of the test
"""

from plasmidcanvas.plasmid import Plasmid
from plasmidcanvas.feature import ArrowFeature, RectangleFeature, RestrictionSite

def get_overlapping_plasmid():
    ol_plasmid = Plasmid("pBR322_2", 4361)

    # Adding an arrow
    # for pBR322 this is TcR
    tcr = ArrowFeature("TcR", 86,1276)
    # # # Customise the thinkness of the line relative to the thickness of the plasmid circle
    # # tcr.set_line_width_scale_factor(1.0)
    ol_plasmid.add_feature(tcr)

    # # Add rop protein for pBR322
    rop = ArrowFeature("rop", 1915,2106)
    ol_plasmid.add_feature(rop)

    # # Add a rectangle, base of mobility for pBR322
    bom = RectangleFeature("bom", 2208,2348)
    ol_plasmid.add_feature(bom)

    # # Add ori
    ori = ArrowFeature("ori", 2534, 3122, -1)
    ori.set_color("orange")
    ol_plasmid.add_feature(ori)

    # # Add ampr - technically this arrow should have a portion segmented for its signal sequence
    ampr = ArrowFeature("ampr", 3293, 4153, -1)
    ampr.set_color("red")
    ol_plasmid.add_feature(ampr)

    # # Add ampr promoter as an arrow
    ampr_promoter = ArrowFeature("ampr promoter", 4154, 4258, -1)
    ampr_promoter.set_color("darkred")
    ol_plasmid.add_feature(ampr_promoter)

    overlapping = ArrowFeature("overlapping feature", 3500, 4300)
    overlapping.set_color("darkblue")
    ol_plasmid.add_feature(overlapping)

    overlapping = ArrowFeature("overlapping feature2", 3366, 3440)
    overlapping.set_color("darkgreen")
    ol_plasmid.add_feature(overlapping)

    overlapping = ArrowFeature("overlapping feature3", 3400, 3800)
    overlapping.set_color("darkgreen")
    ol_plasmid.add_feature(overlapping)

    overlapping = ArrowFeature("overlapping feature4", 2900, 3100)
    overlapping.set_color("darkgreen")
    ol_plasmid.add_feature(overlapping)

    overlapping = ArrowFeature("overlapping feature5", 3600, 3700)
    overlapping.set_color("darkgreen")
    ol_plasmid.add_feature(overlapping)

    overlapping = RectangleFeature("overlapping feature6", 2600, 3200)
    overlapping.set_color("darkgreen")
    ol_plasmid.add_feature(overlapping)

    return ol_plasmid

#get_overlapping_plasmid().save_to_file("overlapping_plasmid.png")