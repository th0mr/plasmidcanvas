
from plasmidcanvas.plasmid import Plasmid
from plasmidcanvas.feature import ArrowFeature, RectangleFeature, RestrictionSite

# Define a plasmid of X base pairs long, with a name
plasmid = Plasmid("pBR322", 4361)
plasmid.set_marker_style("auto")
plasmid.set_feature_label_font_size(7)
plasmid.set_plasmid_line_width_sf(1.25)

# Adding tcr
tcr = ArrowFeature("tcr", 86,1276)
plasmid.add_feature(tcr)

# Add rop protein for pBR322
rop = ArrowFeature("rop", 1915,2106)
rop.set_line_width_scale_factor(1.5)
rop.set_color("purple")
plasmid.add_feature(rop)

# Add a rectangle, base of mobility for pBR322
bom = RectangleFeature("bom", 2208,2348)
plasmid.add_feature(bom)

# Add ori
ori = ArrowFeature("ori", 2534, 3122, -1)
ori.set_color("orange")
plasmid.add_feature(ori)

# # Add ampr
ampr = ArrowFeature("ampr", 3293, 4153, -1)
ampr.set_color("red")
plasmid.add_feature(ampr)

for feature in plasmid.get_features():
    feature.set_label_styles(["on-circle"])

# # Add ampr promoter as an arrow
ampr_promoter = ArrowFeature("ampr promoter", 4154, 4258, -1)
ampr_promoter.set_color("darkred")

ampr_promoter.set_line_width_scale_factor(0.75)
plasmid.add_feature(ampr_promoter)

# Add the sites to the plasmid
plasmid.add_feature(RestrictionSite("BamHI", 375))
plasmid.add_feature(RestrictionSite("BfuAI - BspMI", 1054))
plasmid.add_feature(RestrictionSite("Bpu10I", 1581))
plasmid.add_feature(RestrictionSite("AflIII - PciI", 2473))
plasmid.add_feature(RestrictionSite("AhdI", 3366))

# Plot the plasmid
plasmid.save_to_file("myplasmid")