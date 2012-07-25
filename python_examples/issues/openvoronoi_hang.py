
# 2012-07-25: this reportedly hangs.

import openvoronoi

points = ((-0.7000000000000002, -0.5249999999990376),
        (0.7, -0.5249999999990376),
)

print "OpenVoronoi version: ",openvoronoi.version()

dia = openvoronoi.VoronoiDiagram(0.972222, 2)
for p in points:
    ovp = openvoronoi.Point(*p)
    print "Adding ",ovp
    dia.addVertexSite(ovp)

print dia
print "VD OK?: ",dia.check()

# aewallin, I get:
"""
Adding  (-0.7, -0.525)
Adding  (0.7, -0.525)
VoronoiDiagram 
 num_vertices    = 18
 num_edges       = 31
 num_point_sites = 2
 num_line_sites  = 0
 num_split_vertices  = 0

VD OK?:  True

"""
