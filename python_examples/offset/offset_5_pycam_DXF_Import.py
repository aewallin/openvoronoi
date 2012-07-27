# -*- coding: utf-8 -*-

import math

import openvoronoi
import ovdvtk # functions for drawing VDs with VTK
import offset2vtk # functions for drawing offset-loops with VTK

import pycam.Utils.log
import pycam.Geometry.Line
import pycam.Geometry.Polygon
import pycam.Geometry.Model
#from pycam.Geometry.PointUtils import pdist_sq
from pycam.Geometry.utils import epsilon

_log = pycam.Utils.log.get_logger()

# point_set is a list of 2D points: [ [x0,y0], [x1,y1], ... , [xN,yN] ]
# NOTE: all points in point_set must be unique! i.e. no duplicate points allowed
#
# line_set is a list of line-segments, given as index-pairs into the point_set list
# [ [start0,end0] , [start1,end1], ... , [startM,endM] ]
# NOTE: intersecting line-segments are not allowed.
#
# this defines line-segments   point_set[start0] - points_set[end0]  and so on...
#
# NOTE: currently openvoronoi only supports vertices of degree 2 or lower.
# i.e. a "star" geometry where three or more line-segments connect to a central vertex is forbidden
# SUPPORTED:     point_set = [p0,p1,p2,p3]    line_set = [[0,1], [1,2], [2,3], [3,0]]
# NOT SUPPORTED: point_set = [p0,p1,p2,p3]    line_set = [[0,1], [0,2], [0,3]]   (three line-segments connect to p0!)
def _add_connected_point_set_to_diagram(point_set, line_set, dia):
    # add all points to the diagram
    vpoints = []
    for p in point_set:
        ovp = openvoronoi.Point(*p[:2])
        vpoints.append(dia.addVertexSite(ovp))
    _log.info("all vertices added to openvoronoi!")
    # now add all line-segments
    for segment in line_set:
        start_idx = vpoints[segment[0]]
        end_idx   = vpoints[segment[1]]
        dia.addLineSite(start_idx, end_idx)
    _log.info("all lines added to openvoronoi!")


def _polygons_to_line_set(polygons):
    # all points (unique!)
    point_set = []
    # all line-segments (indexes into the point_set array)
    line_set = []
    previous_point_index = 0
    point_count = 0
    for polygon_index, polygon in enumerate(polygons):
        _log.info("polygon #%d has %d vertices" % (polygon_index, len(polygon)))
        first_point = True
        poly_pts = polygon.get_points()
        # if the polygon is closed, repeat the first point at the end
        if polygon.is_closed and poly_pts:
            poly_pts.append(poly_pts[0])
        for p in poly_pts:
            point_count += 1
            if p not in point_set:
                # this point is a new point we have not seen before
                point_set.append(p)
            current_point_index = point_set.index(p)
            # on the first iteration we have no line-segment
            if not first_point:
                #_log.info(" line from %s to %s" % (previous_point_index, current_point_index))
                line_set.append((previous_point_index, current_point_index))
            else:
                first_point = False
            previous_point_index = current_point_index
    _log.info("point_count: %d" % len(point_set))
    _log.info("point_set size: %d", len(point_set))
    _log.info("number of line-segments: %d" % len(line_set))
    #_log.info("Point set: %s" % str(point_set))
    #_log.info("Line set: %s" % str(line_set))
    return point_set, line_set


# filter diagramn so that either only interior or only exterior or diagram remains valid
# pocket polygons should be input in CW order
# island polygons should be input in CCW order
def _filter_diagram(dia, interior=True):
    polygon_filter = openvoronoi.PolygonInterior( interior )
    dia.filter_graph(polygon_filter)
    return

def _polygon_model_radius(polygons):
    maxx = max([poly.maxx for poly in polygons])
    maxy = max([poly.maxy for poly in polygons])
    minx = min([poly.minx for poly in polygons])
    miny = min([poly.miny for poly in polygons])
    radius = math.sqrt((maxx - minx) ** 2 + (maxy - miny) ** 2) / 1.8
    return radius
    
def pocket_model(polygons, first_offset, second_offset=None, interior=True, N_offsets=1):
    _log.info("number of polygons: %d" % len(polygons))
    _log.info("offset distance: %f" % first_offset)
    if second_offset is None:
        second_offset=first_offset
        
    radius = _polygon_model_radius(polygons) 
    
    _log.info("Radius: %f" % radius)
    bin_size = int(math.ceil(math.sqrt(sum([len(poly.get_points()) for poly in polygons]))))
    _log.info("bin_size: %f" % bin_size)
    dia = openvoronoi.VoronoiDiagram(radius, bin_size)
    point_set, line_set = _polygons_to_line_set(polygons)
    _add_connected_point_set_to_diagram(point_set, line_set, dia)
    _filter_diagram(dia,interior) # choose either interior or exterior to be pocketed
    _log.info("diagram complete")
    _log.info("diagram check: %s" % str(dia.check()))
    offset_dia = openvoronoi.Offset(dia.getGraph())
    _log.info("offset diagram created")
    
    offset_loops = []
    if N_offsets == 1: # one single offset
        _log.info("generating one offset")
        offset_loops = offset_dia.offset(first_offset)
    
    elif N_offsets == -1: # "complete" offsetting of pocket interior
        _log.info("generating complete interior offsets" )        
        offset_loops = []

        current_offset = first_offset 
        offset_loops.extend( offset_dia.offset(current_offset) )
        N = 0
        _log.info(" Offset %d at distance: %f" % (N_offsets, current_offset) )        

        while True:
            N = N + 1
            current_offset = current_offset + second_offset
            current_loops = offset_dia.offset(current_offset)
            if len(current_loops) > 0:
                _log.info(" Offset %d at distance: %f" % (N, current_offset) )        
                offset_loops.extend( current_loops )
            else:
                break # if offset returns no loops then we are done

    
    else: # a specified (>= 2) number of offsets
        _log.info("generating %d offsets" % N_offsets)
        offset_loops = []
        current_offset = first_offset
        while N_offsets > 0:
            current_loops = offset_dia.offset(current_offset)
            if len(current_loops) > 0:
                offset_loops.extend( current_loops )
            N_offsets = N_offsets - 1
            current_offset = current_offset + second_offset
         
    _log.info("got %d loops from openvoronoi" % len(offset_loops))
    return offset_loops


if __name__ == "__main__":
    w=1024
    h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) # a VTK window for drawing 
    ovdvtk.drawOCLtext(myscreen, rev_text=openvoronoi.version() )   # the OpenVoronoi text, revision, and date
    
    # rotate camera for 2D view
    far = 1
    camPos = far
    zmult = 3
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    import sys
    filename = "star2_scaled.dxf"
    offset_distance = 0.1
    sys.argv = ["dummy",filename, offset_distance]
    
    import pycam.Importers.DXFImporter as importer
    print "reading from DXF-file: ",sys.argv[1]
    model = importer.import_model(sys.argv[1])
    model.revise_directions()

    # draw the input geometry
    point_set, line_set = _polygons_to_line_set( model.get_polygons() )
    offset2vtk.drawLinesegs(myscreen, point_set, line_set)
    
    first_offset = 0.03
    second_offset =0.015
    interior = True
    N_offsets = -1
    offset_loops = pocket_model(model.get_polygons(), first_offset,second_offset,interior, N_offsets )
    #offset2vtk.drawOffsets( myscreen, offset_loops)
    
    offset2vtk.drawOffsets2( myscreen, offset_loops)
    
    print "PYTHON All DONE."
    myscreen.render()   
    myscreen.iren.Start()
