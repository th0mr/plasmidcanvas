"""
A copy of example 1 from the usage guide showing how to build pBR322 in plasmidcanvas
This file can be ran individually or as part of the test
"""

from plasmidcanvas.plasmid import Plasmid
from plasmidcanvas.feature import ArrowFeature, RectangleFeature, RestrictionSite

def get_basic_plasmid() -> Plasmid:
    basic_plasmid = Plasmid("pBR322_1", 4361)

    # Adding features
    tcr = ArrowFeature("TcR", 86, 1276)
    basic_plasmid.add_feature(tcr)

    bom = RectangleFeature("bom", 2208,2348)
    basic_plasmid.add_feature(bom)

    ori = ArrowFeature("ori", 2534, 3122, direction=-1)
    basic_plasmid.add_feature(ori)

    ampr = ArrowFeature("ampr", 3293, 4153, direction=-1)
    basic_plasmid.add_feature(ampr)

    ampr_promoter = ArrowFeature("ampr promoter", 4154, 4258, direction=-1)
    basic_plasmid.add_feature(ampr_promoter)

    # Add a couple of restriction sites to the plasmid
    restriction_site_1 = RestrictionSite("BamHI", 375)
    restriction_site_2 = RestrictionSite("BfuAI - BspMI", 1054)
    restriction_site_3 = RestrictionSite("Bpu10I", 1581)
    restriction_site_4 = RestrictionSite("AflIII - PciI", 2473)
    restriction_site_5 = RestrictionSite("AhdI", 3366)

    # Add the sites to the plasmid
    basic_plasmid.add_feature(restriction_site_1)
    basic_plasmid.add_feature(restriction_site_2)
    basic_plasmid.add_feature(restriction_site_3)
    basic_plasmid.add_feature(restriction_site_4)
    basic_plasmid.add_feature(restriction_site_5)

    return basic_plasmid

#get_basic_plasmid().save_to_file("basic_plasmid.png")