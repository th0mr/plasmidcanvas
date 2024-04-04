# An example showing how to build pBR322 in plasmidcanvas

from plasmidcanvas.plasmid import Plasmid
from plasmidcanvas.feature import ArrowFeature, RectangleFeature, RestrictionSite

plasmid = Plasmid("pBR322", 4361)

# Adding features
tcr = ArrowFeature("TcR", 86, 1276)
plasmid.add_feature(tcr)

bom = RectangleFeature("bom", 2208,2348)
plasmid.add_feature(bom)

ori = ArrowFeature("ori", 2534, 3122, direction=-1)
plasmid.add_feature(ori)

ampr = ArrowFeature("ampr", 3293, 4153, direction=-1)
plasmid.add_feature(ampr)

ampr_promoter = ArrowFeature("ampr promoter", 4154, 4258, direction=-1)
plasmid.add_feature(ampr_promoter)

# Add a couple of restriction sites to the plasmid
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

plasmid.save_to_file("pBR322_basic.png")