"""@camvtk docstring
This module provides classes for visualizing CAD/CAM algorithms using VTK.
This module is part of OpenVoronoi.

Copyright 2010-2012 Anders Wallin (anders.e.e.wallin "at" gmail.com)
Published under the GNU General Public License, see http://www.gnu.org/licenses/
"""

import vtk
import time
import datetime
# import ocl
# import pyocl
import openvoronoi as ovd
import math

white = (1,1,1)
black = (0,0,0)
grey = ( float(127)/255,float(127)/255,float(127)/255)

red= (1,0,0)
pink = ( float(255)/255,float(192)/255,float(203)/255)
orange = ( float(255)/255,float(165)/255,float(0)/255)
yellow= (1,1,0)
yellow2= ( float(180)/255,float(255)/255,float(0)/255)
purple=( float(255)/255,float(0)/255,float(176)/255)

green= (0,1,0)
lgreen = ( float(150)/255,float(255)/255,float(150)/255)
dgreen= ( float(21)/255,float(119)/255,float(41)/255)
grass = ( float(182)/255,float(248)/255,float(71)/255)

blue= (0,0,1)
lblue= ( float(125)/255,float(191)/255,float(255)/255 )
blue2= ( float(30)/255,float(144)/255,float(255)/255 )
blue3= ( float(0)/255,float(255)/255,float(163)/255 )
cyan=  (0,1,1)
mag2 =( float(123)/255 , float(35)/255 , float(251)/255  )
magenta = ( float(153)/255 , float(42)/255 , float(165)/255  )


def drawLine(myscreen, pt1, pt2, lineColor):
    myscreen.addActor( Line(p1=(pt1.x,pt1.y,0),p2=(pt2.x,pt2.y,0),color=lineColor) ) 

def drawVertex(myscreen, p, r, vcolor):
    myscreen.addActor( Sphere( center=(p.x,p.y, 0), radius=r, color=vcolor ) )

def drawArc(myscreen, pt1, pt2, r, cen, cw, arcColor, da=0.1):
    # draw arc as many line-segments
    start = pt1-cen
    end = pt2-cen
    theta1 = math.atan2(start.x,start.y)
    theta2 = math.atan2(end.x,end.y)
    alfa=[] # the list of angles
    #da=0.1
    CIRCLE_FUZZ = 1e-9
    # idea from emc2 / cutsim g-code interp G2/G3
    if (cw == False ): 
        while ( (theta2 - theta1) > -CIRCLE_FUZZ): 
            theta2 -= 2*math.pi
    else:
        while( (theta2 - theta1) < CIRCLE_FUZZ): 
            theta2 += 2*math.pi
    
    dtheta = theta2-theta1
    arclength = r*dtheta
    dlength = min(da, arclength/10)
    steps = int( float(arclength) / float(dlength))
    #print "arc subdivision steps: ",steps
    rsteps = float(1)/float(steps)
    dc = math.cos(-dtheta*rsteps) # delta-cos  
    ds = math.sin(-dtheta*rsteps) # delta-sin
    
    previous = pt1
    tr = [start.x, start.y]
    for i in range(steps):
        #f = (i+1) * rsteps #; // varies from 1/rsteps..1 (?)
        #theta = theta1 + i* dtheta
        tr = rotate(tr[0], tr[1], dc, ds) #; // rotate center-start vector by a small amount
        x = cen.x + tr[0] 
        y = cen.y + tr[1] 
        current = ovd.Point(x,y)
        myscreen.addActor( Line(p1=(previous.x,previous.y,0),p2=(current.x,current.y,0),color=arcColor) )
        previous = current 

# rotate by cos/sin. from emc2 gcodemodule.cc
def rotate(x, y,  c,  s):
    tx = x * c - y * s;
    y = x * s + y * c;
    x = tx;
    return [x,y]
    
def drawOffsets(myscreen, ofs):
    # draw loops
    nloop = 0
    lineColor = lgreen
    arcColor = green #grass
    for lop in ofs:
        n = 0
        N = len(lop)
        first_point=[]
        previous=[]
        for p in lop:
            # p[0] is the Point
            # p[1] is -1 for lines, and r for arcs
            if n==0: # don't draw anything on the first iteration
                previous=p[0]
                #first_point = p[0]
            else:
                cw=p[3]
                cen=p[2]
                r=p[1]
                p=p[0]
                if r==-1:
                    drawLine(myscreen, previous, p, lineColor)
                else:
                    drawArc(myscreen, previous, p, r,cen,cw, arcColor)
                #myscreen.addActor( ovdvtk.Line(p1=(previous.x,previous.y,0),p2=(p.x,p.y,0),color=loopColor) )
                previous=p
            n=n+1
        print "rendered loop ",nloop, " with ", len(lop), " points"
        nloop = nloop+1
        

class VD:
    def __init__(self, myscreen, vd, scale=1, textscale=0.06, vertexradius=0.04):
        self.myscreen = myscreen
        self.vd=vd
        self.scale=scale
        
        self.gen_pts=[ovd.Point(0,0)]
        self.generators = PointCloud(pointlist=self.gen_pts)
        self.verts=[]
        self.far=[]
        self.edges =[]
        self.gens =[]
        self.actors=[]
        self.pointsiteColor = yellow
        self.generatorColor = yellow
        self.vertexColor = blue
        self.seedColor = pink
        self.edgeColor = cyan
        self.vertexRadius = vertexradius
        
        self.vdtext  = Text()
        self.vdtext.SetPos( (50, myscreen.height-70) )
        myscreen.addActor(self.vdtext)
        self.actors.append( self.vdtext )
        
        self.vdtext2  = Text()
        self.vdtext2.SetPos( (myscreen.width-500, 50) ) 
        self.vdtext2.SetText( "--" )
        myscreen.addActor(self.vdtext2)
        self.actors.append( self.vdtext2 )
        
        self.gittext  = Text()
        self.gittext.SetPos( (50, 50) )
        self.gittext_text = "github.com/aewallin"
        self.gittext.SetText( self.gittext_text )
        myscreen.addActor(self.gittext)
        
        self.N_pointgen = 0
        self.N_linegen = 0
        self.vdtext_text = "CPU Time:"
        self.setVDText()
        self.drawClearanceDisk = 0
        self.textScale = textscale
        self.drawVertexIndex=1
        self.drawVertices=1
        self.drawGenerators=1
        self.offsetEdges = 0
        self.drawNullEdges = 1
        
    def getActors(self):
        acts=[]
        for a in self.verts:
            acts.append(a)
        for a in self.far:
            acts.append(a)
        for a in self.edges:
            acts.append(a)
        for a in self.gens:
            acts.append(a)
        for a in self.actors:
            acts.append(a)
            
        return acts
        #self.verts=[]
        #self.far=[]
        #self.edges =[]
        #self.gens =[]
        
    def setVertexRadius(self, r):
        self.vertexRadius=r
    
    def drawFarCircle(self, circleColor=magenta):
        r=self.vd.getFarRadius()
        self.myscreen.addActor( Circle( center=(0,0,0), radius=r, color=circleColor ) )
    
    def setVDText(self):
        self.N_pointgen = self.vd.numPointSites()
        self.N_linegen = self.vd.numLineSites()
        self.N_arcgen = self.vd.numArcSites()
        
        #self.vdtext_text = " "
        self.vdtext_text = "Voronoi-Diagram with :\n"
        self.vdtext_text += str(self.N_pointgen) + " point-sites.\n"
        self.vdtext_text += str(self.N_linegen) + " line-sites.\n"
        self.vdtext_text += str(self.N_arcgen) + " arc-sites.\n"
        #self.vdtext_text += "YELLOW = New point-generator/site\n"
        #self.vdtext_text += "PINK = Seed vertex\n"
        #self.vdtext_text += "RED = Delete vertices/edges\n"
        #self.vdtext_text += "GREEN = Modified VD edges\n"
        self.vdtext.SetText( self.vdtext_text )
    
    def setVDText2(self,times):
        self.N_pointgen = self.vd.numPointSites()
        self.N_linegen = self.vd.numLineSites()
        pts = self.N_pointgen
        #print times
        lns = self.N_linegen
        if lns==0 or lns ==1:
            lns=2 # avoid dividing by log(1)=0
        if pts==0 or pts==1:
            pts=2 # avoid dividing by log(1)
        self.vdtext2_text = "VD in {0:.3f} s CPU time:\n".format(sum(times))
        self.vdtext2_text += "N={0} point-sites in {1:.3f} s ".format(pts, times[0])
        self.vdtext2_text += "= {0:.2f} us*N*log2(N) \n".format( 1e6*float( times[0] )/(float(pts)*float(math.log10(pts)/math.log10(2.0))) )
        self.vdtext2_text += "M={0} line-sites in {1:.3f} s ".format(lns, times[1])
        self.vdtext2_text += "= {0:.2f} us*M*log2(M)".format( 1e6*float( times[1] )/(float(lns)*float(math.log10(lns)/math.log10(2.0))) )
        self.vdtext2.SetText( self.vdtext2_text )
        
        
    def setGenerators(self):
        for p in self.gens:
            self.myscreen.removeActor(p)
        self.gens = []
        for pt_data in self.vd.getGenerators():
            pt = pt_data[0] # [ point, dist, status, index ]
            idx = pt_data[3] # index
            
            p = self.scale*pt
            if self.drawGenerators:
                actor = Sphere( center=(p.x,p.y, 0), radius=self.vertexRadius, color=self.generatorColor )
                self.gens.append(actor)
                self.myscreen.addActor( actor )
            
            if self.drawVertexIndex:
                id_text = str(idx)
                factor = FollowerText( text=id_text,center=(p.x,p.y,0), scale = self.textScale, color=self.pointsiteColor)
                factor = FollowerText( text=id_text,center=(p.x,p.y,0), scale = self.textScale, color=self.pointsiteColor)
                self.gens.append(factor)
                self.myscreen.addActor( factor )

        self.myscreen.render() 
    
    def setFar(self, vd):
        for pt in vd.getFarVoronoiVertices():
            p=self.scale*pt[0]
            self.myscreen.addActor( Sphere( center=(p.x,p.y, 0), radius=self.vertexRadius, color=pink ) )
            cir_actor = Circle( center=(p.x,p.y,0), radius=(pt[1])*self.scale, color=self.vertexColor )
            #self.verts.append(cir_actor)
            self.myscreen.addActor(cir_actor)
            
        self.myscreen.render() 
    
    def setVertices(self):
        for p in self.verts:
            self.myscreen.removeActor(p) # remove old actors
        self.verts = []
        for pt in self.vd.getVoronoiVertices(): # create new actors
            # [ position, dist, status, index ]
            p = self.scale*pt[0]
            vcolor = cyan
            status = pt[2]
            idx = pt[3]
            
            if status == ovd.VertexStatus.IN:
                vcolor = red
            elif status == ovd.VertexStatus.NEW:
                vcolor = green
            elif status == ovd.VertexStatus.OUT:
                vcolor = blue
            
            if self.drawVertices:
                actor = Sphere( center=(p.x,p.y, 0), radius=self.vertexRadius, color=vcolor )
                self.verts.append(actor)
                self.myscreen.addActor( actor )
                        
            if self.drawVertexIndex:
                id_text = str(idx)
                factor = FollowerText( text=id_text,center=(p.x,p.y,0), scale = self.textScale, color=vcolor)
                self.verts.append(factor)
                self.myscreen.addActor( factor )
            
            if self.drawClearanceDisk:
                cir_actor = Circle( center=(p.x,p.y,0), radius=pt[1]*self.scale, color=self.vertexColor )
                self.verts.append(cir_actor)
                self.myscreen.addActor(cir_actor)
        self.myscreen.render() 

    def drawErrorVertices(self):
        #print "drawVertexIdx ", index
        for pt in self.vd.getVoronoiVertices():
            p = pt[0] # self.scale*
            vcolor = cyan
            status = pt[2]
            idx = pt[3]
            # 4 type
            err = pt[5]
            if err > 1e-6:
                #id_text = str(idx)
                print "drawErrorVertex ", idx, " pos= ", p, " err=",err
                #factor = FollowerText( text=id_text,center=(p.x,p.y,0), scale = self.textScale, color=vcolor)
                #self.verts.append(factor)
                actor = Sphere( center=(p.x,p.y, 0), radius=self.vertexRadius, color= pink )
                self.verts.append(actor)
                self.myscreen.addActor( actor )
                
                #self.myscreen.addActor( factor )
        
    
    def drawVertexIdx(self, index):
        #print "drawVertexIdx ", index
        for pt in self.vd.getVoronoiVertices():
            p = self.scale*pt[0]
            vcolor = cyan
            status = pt[2]
            idx = pt[3]
            if idx == index:
                id_text = str(idx)
                #print "drawVertexIdx text= ", id_text, " pos= ", p
                factor = FollowerText( text=id_text,center=(p.x,p.y,0), scale = self.textScale, color=vcolor)
                #self.verts.append(factor)
                self.myscreen.addActor( factor )
                
    def drawIncidentVertexIds(self):
        for pt in self.vd.getVoronoiVertices():
            p = self.scale*pt[0]
            vcolor = cyan
            status = pt[2]
            idx = pt[3]
            vtype = pt[4]
            if status == ovd.VertexStatus.IN:
                vcolor = red
            elif status == ovd.VertexStatus.OUT:
                vcolor = blue
            elif status == ovd.VertexStatus.NEW:
                vcolor = green
                    
            if status == ovd.VertexStatus.IN or status == ovd.VertexStatus.OUT or status == ovd.VertexStatus.NEW:
                if ( vtype != ovd.VertexType.SEPPOINT and vtype != ovd.VertexType.ENDPOINT and vtype != ovd.VertexType.POINTSITE):
                    id_text = str(idx)
                    #print "drawVertexIdx text= ", id_text, " pos= ", p
                    factor = FollowerText( text=id_text,center=(p.x,p.y,0), scale = self.textScale, color=vcolor)
                    #self.verts.append(factor)
                    self.myscreen.addActor( factor )
        
    def setEdgesPolydata(self):
        for e in self.edges:
            self.myscreen.removeActor(e)
        self.edges = []
        self.edgePoints = vtk.vtkPoints()
        self.lineCells=vtk.vtkCellArray()
        self.colorLUT = vtk.vtkLookupTable()
        idx = 0
        last_idx = 0
        vd_edges=[]
        if self.offsetEdges == 0:
            vd_edges = self.vd.getVoronoiEdges()
        else:
            vd_edges = self.vd.getVoronoiEdgesOffset()
            
        for e in vd_edges:
            epts  = e[0]  # points
            etype = e[1]  # type

            first = 1
            segs=[]
            for p in epts:
                self.edgePoints.InsertNextPoint( p.x, p.y, 0)
                if first==0:
                    seg = [last_idx,idx]
                    segs.append(seg)
                first = 0
                last_idx = idx
                idx = idx + 1
                
            # create line and cells
            for seg in segs:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, seg[0])
                line.GetPointIds().SetId(1, seg[1])
                #print " indexes: ", seg[0]," to ",seg[1]
                self.lineCells.InsertNextCell(line)
        Colors = vtk.vtkUnsignedCharArray()


        self.colorLUT = vtk.vtkUnsignedCharArray()
        self.colorLUT.SetNumberOfComponents(3)
        self.colorLUT.SetName("Colors")
        # go through edges once more and set colors
        m=0
        for e in vd_edges:
            # e[1]  edgeType
            src_status = e[2] # src status
            trg_status = e[3] # trg status
            ecolor = self.edgeTypeColor( e[1], e[2], e[3] ) 
            for dummy in e[0]: # go through all the points
                self.colorLUT.InsertNextTuple3( 255*ecolor[0], 255*ecolor[1], 255*ecolor[2])
                m=m+1

        linePolyData = vtk.vtkPolyData()
        linePolyData.SetPoints(self.edgePoints)
        linePolyData.SetLines(self.lineCells)
        linePolyData.GetPointData().SetScalars(self.colorLUT)
        linePolyData.Modified() 
        linePolyData.Update()
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(linePolyData)
        self.edge_actor = vtk.vtkActor()
        self.edge_actor.SetMapper(mapper)
        self.myscreen.addActor( self.edge_actor )
        self.edges.append(self.edge_actor)
        
    def edgeStatusColor(self, src_status, trg_status, default_color):
        
        if ( (src_status == ovd.VertexStatus.IN) and (trg_status == ovd.VertexStatus.IN)):
            return red
        elif ( ((src_status == ovd.VertexStatus.OUT) and 
             (trg_status == ovd.VertexStatus.IN) ) or 
             ((src_status == ovd.VertexStatus.IN) and 
             (trg_status == ovd.VertexStatus.OUT) ) ):
            return purple
        elif ( ((src_status == ovd.VertexStatus.OUT) and 
             (trg_status == ovd.VertexStatus.NEW) ) or 
             ((src_status == ovd.VertexStatus.NEW) and 
             (trg_status == ovd.VertexStatus.OUT) ) ):
            return default_color
        elif ( ((src_status == ovd.VertexStatus.IN) and 
             (trg_status == ovd.VertexStatus.NEW) ) or 
             ((src_status == ovd.VertexStatus.NEW) and 
             (trg_status == ovd.VertexStatus.IN) ) ):
            return red    
        elif ( (src_status == ovd.VertexStatus.NEW) and (trg_status == ovd.VertexStatus.NEW)):
            return green
        else:
            return default_color
            
    
    def edgeTypeColor(self, edgeType, src_status, trg_status):
        #print " edgeStatusColor edgeType= ",edgeType
        if (edgeType == ovd.EdgeType.LINELINE):
            return self.edgeStatusColor(src_status,trg_status,lblue)
        if (edgeType == ovd.EdgeType.PARA_LINELINE):
            return self.edgeStatusColor(src_status,trg_status,blue)
        if (edgeType == ovd.EdgeType.LINE):
            return  self.edgeStatusColor(src_status,trg_status,cyan)
        if (edgeType == ovd.EdgeType.OUTEDGE):
            return pink
        if (edgeType == ovd.EdgeType.SEPARATOR):
            return self.edgeStatusColor(src_status,trg_status, mag2)
        if (edgeType == ovd.EdgeType.LINESITE):
            return yellow
        if (edgeType == ovd.EdgeType.ARCSITE):
            return yellow2
        if (edgeType == ovd.EdgeType.PARABOLA):
            return self.edgeStatusColor(src_status,trg_status, blue2)
        if (edgeType == ovd.EdgeType.HYPERBOLA):
            return self.edgeStatusColor(src_status,trg_status, blue3)
        if (edgeType == ovd.EdgeType.NULLEDGE):
            return white
        else:
            print "UNKOWN edge type = ", edgeType
            return white
            
    def setEdges(self):
        for e in self.edges:
            self.myscreen.removeActor(e)
        self.edges = []
        for e in self.vd.getVoronoiEdges():
            epts  = e[0]  # points
            etype = e[1]  # type
            src_status = e[2] # src status
            trg_status = e[3] # trg status
            ecolor = pink # self.edgeColor
            #print "drawing etype=", etype
            actor = 0
            
            if (etype == ovd.VoronoiEdgeType.LINELINE):
                ecolor = self.edgeStatusColor(src_status,trg_status,lblue)
                for n in range( len(epts)-1 ):
                    p1 = self.scale*epts[n]  
                    p2 = self.scale*epts[n+1] 
                    #print "line ",n," : ",p1," to ",p2
                    actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                    self.myscreen.addActor(actor)
                    self.edges.append(actor)
            elif (etype == ovd.VoronoiEdgeType.LINE):
                ecolor = self.edgeStatusColor(src_status,trg_status,cyan)
                for n in range( len(epts)-1 ):
                    p1 = self.scale*epts[n]  
                    p2 = self.scale*epts[n+1] 
                    actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                    self.myscreen.addActor(actor)
                    self.edges.append(actor)
            elif (etype == ovd.VoronoiEdgeType.OUTEDGE):
                ecolor = pink
                p1 = self.scale*epts[0]  
                p2 = self.scale*epts[1] 
                actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                self.myscreen.addActor(actor)
                self.edges.append(actor)
            elif (etype == ovd.VoronoiEdgeType.SEPARATOR):
                ecolor = self.edgeStatusColor(src_status,trg_status, orange)
                p1 = self.scale*epts[0]  
                p2 = self.scale*epts[1] 
                actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                self.myscreen.addActor(actor)
                self.edges.append(actor)
            elif (etype == ovd.VoronoiEdgeType.LINESITE):
                ecolor = yellow
                p1 = self.scale*epts[0]  
                p2 = self.scale*epts[1] 
                actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                self.myscreen.addActor(actor)
                self.edges.append(actor)
            elif (etype == ovd.VoronoiEdgeType.PARABOLA):
                ecolor = self.edgeStatusColor(src_status,trg_status, blue2)
                eactor = PolyLine(pointList=epts,color=ecolor)
                self.myscreen.addActor(eactor)
                self.edges.append(eactor)
            elif (etype == ovd.VoronoiEdgeType.ARCSITE):
                ecolor = yellow #self.edgeStatusColor(src_status,trg_status, blue2)
                eactor = PolyLine(pointList=epts,color=ecolor)
                self.myscreen.addActor(eactor)
                self.edges.append(eactor)

        self.myscreen.render() 
        
    def setAll(self):
        self.setVDText()
        self.setGenerators()
        #self.setFar(vd)
        self.setVertices()
        self.setEdgesPolydata()
        #self.setEdges()
        

def drawOCLtext(myscreen, rev_text=" "):
    t = Text()
    t.SetPos( (myscreen.width-250, myscreen.height-70) )
    date_text = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    #rev_text = ovd.version()
    t.SetText( "OpenVoronoi\n" + rev_text + "\n" + date_text )
    myscreen.addActor(t)

def drawBB( myscreen, vol ):
    """ draw a bounding-box """
    lines = []
    lines.append( Line(p1=(vol.bb.minx, vol.bb.miny, vol.bb.minz) , p2=(vol.bb.maxx, vol.bb.miny, vol.bb.minz))  )
    lines.append( Line(p1=(vol.bb.minx, vol.bb.maxy, vol.bb.minz) , p2=(vol.bb.maxx, vol.bb.maxy, vol.bb.minz))  )
    lines.append( Line(p1=(vol.bb.minx, vol.bb.maxy, vol.bb.maxz) , p2=(vol.bb.maxx, vol.bb.maxy, vol.bb.maxz))  )
    lines.append( Line(p1=(vol.bb.minx, vol.bb.miny, vol.bb.maxz) , p2=(vol.bb.maxx, vol.bb.miny, vol.bb.maxz))  )
    
    lines.append( Line(p1=(vol.bb.minx, vol.bb.miny, vol.bb.minz) , p2=(vol.bb.minx, vol.bb.miny, vol.bb.maxz))  )
    lines.append( Line(p1=(vol.bb.maxx, vol.bb.miny, vol.bb.minz) , p2=(vol.bb.maxx, vol.bb.miny, vol.bb.maxz))  )
    
    lines.append( Line(p1=(vol.bb.minx, vol.bb.maxy, vol.bb.minz) , p2=(vol.bb.minx, vol.bb.maxy, vol.bb.maxz))  )
    lines.append( Line(p1=(vol.bb.maxx, vol.bb.maxy, vol.bb.minz) , p2=(vol.bb.maxx, vol.bb.maxy, vol.bb.maxz))  )
    
    lines.append( Line(p1=(vol.bb.minx, vol.bb.miny, vol.bb.minz) , p2=(vol.bb.minx, vol.bb.maxy, vol.bb.minz))  )
    lines.append( Line(p1=(vol.bb.maxx, vol.bb.miny, vol.bb.minz) , p2=(vol.bb.maxx, vol.bb.maxy, vol.bb.minz))  )
    
    lines.append( Line(p1=(vol.bb.minx, vol.bb.miny, vol.bb.maxz) , p2=(vol.bb.minx, vol.bb.maxy, vol.bb.maxz))  )
    lines.append( Line(p1=(vol.bb.maxx, vol.bb.miny, vol.bb.maxz) , p2=(vol.bb.maxx, vol.bb.maxy, vol.bb.maxz))  )    
    
    for l in lines:
        myscreen.addActor(l)

def drawTree(myscreen,t,color=red,opacity=0.2, offset=(0,0,0)):
    """ draw an octree """
    nodes = t.get_nodes()
    #nmax=len(nodes)
    #i=0
    for n in nodes:
        cen = n.point() # center of cube
        scale = n.get_scale() # scale of cube
        cube = camvtk.Cube(center=(cen.x+offset[0], cen.y+offset[1], cen.z+offset[2]), length= scale, color=color)
        cube.SetOpacity(opacity)
        #cube.SetPhong()
        cube.SetGouraud()
        #cube.SetWireframe()
        myscreen.addActor( cube )
        #if (nmax>100):
        #    print "i=", i
        #    print "div=", (float(nmax)/10)
        #    if ( (i % (float(nmax)/10))==0):
        #        print ".",
        #i=i+1
    #print "done."



def drawTree2(myscreen,t,color=red,opacity=0.2):
    """ draw an octree as an STLSurface """
    tlist = pyocl.octree2trilist(t)
    surf = STLSurf(triangleList=tlist)
    surf.SetColor(color)
    surf.SetOpacity(opacity)
    myscreen.addActor(surf)
    
def drawArrows(myscreen,center=(0,0,0)):
    # X Y Z arrows
    arrowcenter=center
    xar = Arrow(color=red,   center=arrowcenter, rotXYZ=(0,0,0))
    yar = Arrow(color=green, center=arrowcenter, rotXYZ=(0,0,90))
    zar = Arrow(color=blue,  center=arrowcenter, rotXYZ=(0,-90,0))
    myscreen.addActor(xar)
    myscreen.addActor(yar)
    myscreen.addActor(zar)

def drawCylCutter(myscreen, c, p):
    cyl = Cylinder(center=(p.x,p.y,p.z), radius=c.radius,
                            height=c.length,
                            rotXYZ=(90,0,0), color=grey)
    cyl.SetWireframe()
    myscreen.addActor(cyl) 

def drawBallCutter(myscreen, c, p):
    cyl = Cylinder(center=(p.x,p.y,p.z+c.getRadius() ), radius=c.getRadius(),
                            height=c.getLength(),
                            rotXYZ=(90,0,0), color=red)
    #cyl.SetWireframe()
    sph = Sphere(center=(p.x,p.y,p.z+c.getRadius()), radius=c.getRadius(), color=red)
    myscreen.addActor(cyl) 
    myscreen.addActor(sph)
    acts=[]
    acts.append(cyl)
    acts.append(sph)
    return acts 


class VTKScreen():
    """
    a vtk render window for displaying geometry
    """
    def __init__(self, width=1280, height=720): # 1080p:  1920x1080
        """ create a screen """
        self.width=width
        self.height=height

        self.ren = vtk.vtkRenderer()
        self.renWin = vtk.vtkRenderWindow()
        self.renWin.AddRenderer(self.ren)
        self.renWin.SetSize(self.width,self.height)
        
        self.iren = vtk.vtkRenderWindowInteractor()
        self.iren.SetRenderWindow(self.renWin)
        interactorstyle = self.iren.GetInteractorStyle()
        interactorstyle.SetCurrentStyleToTrackballCamera()     
           
        self.camera = vtk.vtkCamera()
        self.camera.SetClippingRange(0.01, 1000)
        self.camera.SetFocalPoint(0, 0, 0)
        self.camera.SetPosition(0, 35, 5)
        self.camera.SetViewAngle(30)
        self.camera.SetViewUp(0, 0, 1)
        self.ren.SetActiveCamera(self.camera)
        self.iren.Initialize()
        
        
    def setAmbient(self, r, g, b):
        """ set ambient color """
        self.ren.SetAmbient(r, g, b)
                    
    def addActor(self, actor):
        """ add an actor """
        self.ren.AddActor(actor)
    
    def removeActor(self, actor):
        """ remove an actor"""
        #actor.Delete()
        self.ren.RemoveActor(actor)
        

    def render(self):
        """ render scene"""
        self.renWin.Render()
        
    def GetLights(self):
        return self.ren.GetLights()
    def CreateLight(self):
        self.ren.CreateLight()
    def MakeLight(self):
        return self.ren.MakeLight()
    def AddLight(self,l):
        self.ren.AddLight(l)
    def RemoveAllLights(self):
        self.ren.RemoveAllLights()
    def SetLightCollection(self,lights):
        self.ren.SetLightCollection(lights)
    def Close(self):
        self.iren.TerminateApp()
        

class CamvtkActor(vtk.vtkActor):
    """ base class for actors"""
    def __init__(self):
        """ do nothing"""
        pass
    
    def Delete(self):
        self.Delete()
    
    def SetColor(self, color):
        """ set color of actor"""
        self.GetProperty().SetColor(color)
    
    def SetOpacity(self, op=0.5):
        """ set opacity of actor, 0 is see-thru (invisible)"""
        self.GetProperty().SetOpacity(op)   
    
    def SetWireframe(self):
        """ set surface to wireframe"""
        self.GetProperty().SetRepresentationToWireframe()
        
    def SetSurface(self):
        """ set surface rendering on"""
        self.GetProperty().SetRepresentationToSurface() 
        
    def SetPoints(self):
        """ render only points"""
        self.GetProperty().SetRepresentationToPoints()
        
    def SetFlat(self):     
        """ set flat shading"""
        self.GetProperty().SetInterpolationToFlat()
    
    def SetGouraud(self):
        """ set gouraud shading"""
        self.GetProperty().SetInterpolationToGouraud()
    
    def SetPhong(self):
        """ set phong shading"""
        self.GetProperty().SetInterpolationToPhong()

    # possible TODOs
    # specular
    # diffuse
    # ambient



class FollowerText(vtk.vtkFollower):
    """ 3D text """
    def __init__(self,text="test",color=cyan,center=(0,0,0),scale=1):
        self.textSource = vtk.vtkVectorText()
        self.textSource.SetText( text )
        self.scale = scale

        self.transform = vtk.vtkTransform()
        
        self.transform.Translate(center[0], center[1], center[2] )
        self.transform.Scale(self.scale, self.scale, self.scale)
        self.transformFilter=vtk.vtkTransformPolyDataFilter()
        self.transformFilter.SetTransform(self.transform)
        self.transformFilter.SetInputConnection(self.textSource.GetOutputPort())
        self.transformFilter.Update()
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection( self.transformFilter.GetOutputPort() )
        self.SetMapper(self.mapper)
        self.SetColor(color)
        
    def SetScale(self, scale):
        self.scale = scale
        self.transform.Scale(self.scale, self.scale, self.scale)
        self.transformFilter.Update()
        
    def SetText(self,text):
        self.textSource.SetText( text )
        
    def SetColor(self, color):
        """ set color of actor"""
        self.GetProperty().SetColor(color)
        
        
class Cone(CamvtkActor):
    """ a cone"""
    def __init__(self,  center=(-2,0,0), radius = 1, angle=45, height=0.4, color=(1,1,0) , resolution=60):
        """ cone"""
        self.src = vtk.vtkConeSource()
        self.src.SetResolution(resolution)
        self.src.SetRadius( radius ) 
        #self.src.SetAngle( angle )
        self.src.SetHeight( height )
        #self.src.SetCenter(center)
        
        transform = vtk.vtkTransform()
        transform.Translate(center[0], center[1], center[2] - self.src.GetHeight()/2)
        #transform.RotateX(rotXYZ[0])
        transform.RotateY( -90 )
        #transform.RotateZ(rotXYZ[2])
        transformFilter=vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(self.src.GetOutputPort())
        transformFilter.Update()
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(transformFilter.GetOutput())
        
        
        #self.mapper = vtk.vtkPolyDataMapper()
        #self.mapper.SetInput(self.src.GetOutput())
        self.SetMapper(self.mapper)
        self.SetColor(color)

class Sphere(CamvtkActor):
    """ a sphere"""
    def __init__(self, radius=1, resolution=20, center=(0,2,0),
                color=(1,0,0)):
        """ create sphere"""
        self.src = vtk.vtkSphereSource()
        self.src.SetRadius(radius)
        self.src.SetCenter(center)
        self.src.SetThetaResolution(resolution)
        self.src.SetPhiResolution(resolution)

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src.GetOutput())
        self.SetMapper(self.mapper)
        self.SetColor(color)

class Cube(CamvtkActor):
    """ a cube"""
    def __init__(self,center=(2,2,0) , length=1, color=(0,1,0) ):
        """ create cube"""
        self.src = vtk.vtkCubeSource()
        self.src.SetCenter(center)
        self.src.SetXLength(length)
        self.src.SetYLength(length)
        self.src.SetZLength(length)
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src.GetOutput())
        self.SetMapper(self.mapper)
        self.SetColor(color)

class Cylinder(CamvtkActor):
    """ cylinder """
    def __init__(self,center=(0,-2,0) , radius=0.5, height=2, color=(0,1,1),
                    rotXYZ=(0,0,0), resolution=50 ):
        """ cylinder """
        self.src = vtk.vtkCylinderSource()
        self.src.SetCenter(0,0,0)
        self.src.SetHeight(height)
        self.src.SetRadius(radius)
        self.src.SetResolution(resolution)
        # SetResolution
        # SetCapping(int)
        # CappingOn() CappingOff()
        
        # this transform rotates the cylinder so it is vertical
        # and then translates the lower tip to the center point
        transform = vtk.vtkTransform()
        transform.Translate(center[0], center[1], center[2]+height/2)
        transform.RotateX(rotXYZ[0])
        transformFilter=vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(self.src.GetOutputPort())
        transformFilter.Update()

        
        self.mapper = vtk.vtkPolyDataMapper()
        #self.mapper.SetInput(self.src.GetOutput())
        self.mapper.SetInput( transformFilter.GetOutput() )
        self.SetMapper(self.mapper)
        self.SetColor(color)


class Line(CamvtkActor):
    """ line """
    def __init__(self,p1=(0,0,0) , p2=(1,1,1), color=(0,1,1) ):   
        """ line """
        self.src = vtk.vtkLineSource()
        self.src.SetPoint1(p1)
        self.src.SetPoint2(p2)
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src.GetOutput())
        self.SetMapper(self.mapper)
        self.SetColor(color)

class PolyLine(CamvtkActor):
    def __init__(self, pointList=[], color=(1,1,1) ):
        self.src=[]
        points = vtk.vtkPoints()
        polyline = vtk.vtkCellArray()
        
        idx = 0
        first = 1
        last_idx = 0
        segs=[]
        for p in pointList:
            points.InsertNextPoint(p.x, p.y, 0)
            #print "p = ",p
            if first==0:
                seg = [last_idx,idx]
                segs.append(seg)
            first = 0
            last_idx = idx
            idx = idx + 1

        for seg in segs:
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, seg[0])
            line.GetPointIds().SetId(1, seg[1])
            #print " indexes: ", seg[0]," to ",seg[1]
            polyline.InsertNextCell(line)
            

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetLines(polyline)
        polydata.Modified()
        polydata.Update()
        self.src=polydata
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src)
        self.SetMapper(self.mapper)
        self.SetColor(color)
        polydata.Modified()
        polydata.Update()

        # SetScaleFactor(double)
        # GetOrigin


class Tube(CamvtkActor):
    """ line with tube filter"""
    def __init__(self,p1=(0,0,0) , p2=(1,1,1), radius=0.1, color=(0,1,1) ):   
        self.src = vtk.vtkLineSource()
        self.src.SetPoint1(p1)
        self.src.SetPoint2(p2)
        
        self.tubefilter = vtk.vtkTubeFilter()
        self.tubefilter.SetInput( self.src.GetOutput() )
        self.tubefilter.SetRadius( radius )
        self.tubefilter.SetNumberOfSides( 30 )
        self.tubefilter.Update()
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection( self.tubefilter.GetOutputPort() )
        self.SetMapper(self.mapper)
        self.SetColor(color)


class Circle(CamvtkActor):
    """ circle"""
    def __init__(self,center=(0,0,0) , radius=1, color=(0,1,1), resolution=50 ):   
        """ create circle """
        lines =vtk.vtkCellArray()
        id = 0
        points = vtk.vtkPoints()
        for n in xrange(0,resolution):
            line = vtk.vtkLine()
            angle1 = (float(n)/(float(resolution)))*2*math.pi
            angle2 = (float(n+1)/(float(resolution)))*2*math.pi
            p1 = (center[0]+radius*math.cos(angle1), center[1]+radius*math.sin(angle1), center[2])
            p2 = (center[0]+radius*math.cos(angle2), center[1]+radius*math.sin(angle2), center[2])
            points.InsertNextPoint(p1)
            points.InsertNextPoint(p2)
            line.GetPointIds().SetId(0,id)
            id=id+1
            line.GetPointIds().SetId(1,id)
            id=id+1
            lines.InsertNextCell(line)
            

        self.pdata = vtk.vtkPolyData()
        self.pdata.SetPoints(points)
        self.pdata.SetLines(lines)
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.pdata)
        self.SetMapper(self.mapper)
        self.SetColor(color)
        
class Tube(CamvtkActor):
    """ a Tube is a line with thickness"""
    def __init__(self, p1=(0,0,0) , p2=(1,1,1), radius=0.2, color=(0,1,1) ):   
        """ tube"""
        points = vtk.vtkPoints()
        points.InsertNextPoint(p1)
        points.InsertNextPoint(p2)
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0,0)
        line.GetPointIds().SetId(1,1)
        lines =vtk.vtkCellArray()
        lines.InsertNextCell(line)
        self.pdata = vtk.vtkPolyData()
        self.pdata.SetPoints(points)
        self.pdata.SetLines(lines)

        tubefilter=vtk.vtkTubeFilter()
        tubefilter.SetInput(self.pdata)
        tubefilter.SetRadius(radius)
        tubefilter.SetNumberOfSides(50)
        tubefilter.Update()
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(tubefilter.GetOutput())
        self.SetMapper(self.mapper)
        self.SetColor(color)


class Point(CamvtkActor):
    """ point"""
    def __init__(self, center=(0,0,0), color=(1,2,3) ):   
        """ create point """
        self.src = vtk.vtkPointSource()
        self.src.SetCenter(center)
        self.src.SetRadius(0)
        self.src.SetNumberOfPoints(1)

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src.GetOutput())
        self.SetMapper(self.mapper)
        self.SetColor(color)

class Arrow(CamvtkActor):
    """ arrow """
    def __init__(self, center=(0,0,0), color=(0,0,1), rotXYZ=(0,0,0) ):
        """ arrow """
        self.src = vtk.vtkArrowSource()
        #self.src.SetCenter(center)
        
        transform = vtk.vtkTransform()
        transform.Translate(center[0], center[1], center[2])
        transform.RotateX(rotXYZ[0])
        transform.RotateY(rotXYZ[1])
        transform.RotateZ(rotXYZ[2])
        transformFilter=vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(self.src.GetOutputPort())
        transformFilter.Update()
        

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput( transformFilter.GetOutput() )
        self.SetMapper(self.mapper)
        self.SetColor(color)


class Text(vtk.vtkTextActor):
    """ 2D text, HUD-type"""
    def __init__(self, text="text",size=18,color=(1,1,1),pos=(100,100)):
        """create text"""
        self.SetText(text)
        self.properties=self.GetTextProperty()
        self.properties.SetFontFamilyToArial()
        self.properties.SetFontSize(size)
        
        self.SetColor(color)
        self.SetPos(pos)
    
    def SetColor(self,color):
        """ set color of text """
        self.properties.SetColor(color)
    
    def SetPos(self, pos):
        """ set position on screen """
        self.SetDisplayPosition(pos[0], pos[1])

    def SetText(self, text):
        """ set text to be displayed """
        self.SetInput(text)
        
    def SetSize(self, size):
        self.properties.SetFontSize(size)
        
class Text3D(vtk.vtkFollower):
    """ 3D text rendered in the scene"""
    def __init__(self, color=(1,1,1), center=(0,0,0), text="hello", scale=1, camera=[]):
        """ create text """
        self.src = vtk.vtkVectorText()
        self.SetText(text)
        #self.SetCamera(camera)
        transform = vtk.vtkTransform()
        
        transform.Translate(center[0], center[1], center[2])
        transform.Scale(scale, scale, scale)
        #transform.RotateY(90)
        #transform2 = vtk.vtkTransform()
        #transform.Concatenate(transform2)
        #transformFilter=vtk.vtkTransformPolyDataFilter()
        #transformFilter.SetTransform(transform)
        #transformFilter.SetInputConnection(self.src.GetOutputPort())
        #transformFilter.Update()
        
        #follower = vtk.vtkFollower()
        #follower.SetMapper
        
        self.SetUserTransform(transform)
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.src.GetOutputPort())
        self.SetMapper(self.mapper)
        self.SetColor(color)
        
    def SetText(self, text):
        """ set text to be displayed"""
        self.src.SetText(text)
        
    def SetColor(self,color):
        """ set color of text"""
        self.GetProperty().SetColor(color)        

class Axes(vtk.vtkActor):
    """ axes (x,y,z) """
    def __init__(self, center=(0,0,0), color=(0,0,1) ):
        """ create axes """
        self.src = vtk.vtkAxes()
        #self.src.SetCenter(center)

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src.GetOutput())
        self.SetMapper(self.mapper)

        self.SetColor(color)
        self.SetOrigin(center)
        # SetScaleFactor(double)
        # GetOrigin
     

    def SetColor(self, color):
        self.GetProperty().SetColor(color)
    
    def SetOrigin(self, center=(0,0,0)):
        self.src.SetOrigin(center[0], center[1], center[2])

class Toroid(CamvtkActor):
    def __init__(self, r1=1, r2=0.25, center=(0,0,0), rotXYZ=(0,0,0), color=(1,0,0)):
        self.parfun = vtk.vtkParametricSuperToroid()
        self.parfun.SetRingRadius(r1)
        self.parfun.SetCrossSectionRadius(r2)
        self.parfun.SetN1(1)
        self.parfun.SetN2(1)
         
        self.src = vtk.vtkParametricFunctionSource()
        self.src.SetParametricFunction(self.parfun)
        
        transform = vtk.vtkTransform()
        transform.Translate(center[0], center[1], center[2])
        transform.RotateX(rotXYZ[0])
        transform.RotateY(rotXYZ[1])
        transform.RotateZ(rotXYZ[2])
        transformFilter=vtk.vtkTransformPolyDataFilter()
        transformFilter.SetTransform(transform)
        transformFilter.SetInputConnection(self.src.GetOutputPort())
        transformFilter.Update()
        
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(transformFilter.GetOutput())
        self.SetMapper(self.mapper)
        self.SetColor(color)    

"""
class TrilistReader(vtk.vtkPolyDataAlgorithm):
    def __init__(self, triangleList):
        vtk.vtkPolyDataAlgorithm.__init__(self)
        self.FileName = None
        self.SetNumberOfInputPorts(0)
        self.SetNumberOfOutputPorts(1)
        
    def FillOutputPortInfornmation(self, port, info):
        if port == 0:
            info.Set( vtk.vtkDataObject.DATA_TYPE_NAME(), "vtkPolyData")
            return 1
        return 0
        
    def RequestData(self, request, inputVector, outputVector):
        outInfo = outputVector.GetInformationObject(0)
        output = outInfo.Get( vtk.vtkDataObject.DATA_OBJECT() )
        polydata = vtk.vtkPolyData()
        points = vtk.vtkPoints()
        points.InsertNextPoint(0,0,0)
        polydata.SetPoints(points)
        
        output.ShallowCopy(polydata)
        return 1
"""

class STLSurf(CamvtkActor):
    def __init__(self, filename=None, triangleList=[], color=(1,1,1) ):
        self.src=[]
        if filename is None:
            points = vtk.vtkPoints()
            triangles = vtk.vtkCellArray()
            n=0
            for t in triangleList:
                triangle = vtk.vtkTriangle()
                for p in t.getPoints():
                    points.InsertNextPoint(p.x, p.y, p.z)
                triangle.GetPointIds().SetId(0,n)
                n=n+1
                triangle.GetPointIds().SetId(1,n)
                n=n+1
                triangle.GetPointIds().SetId(2,n)
                n=n+1
                triangles.InsertNextCell(triangle)
            polydata= vtk.vtkPolyData()
            polydata.SetPoints(points)
            polydata.SetPolys(triangles)
            polydata.Modified()
            polydata.Update()
            self.src=polydata
            self.mapper = vtk.vtkPolyDataMapper()
            self.mapper.SetInput(self.src)
            self.SetMapper(self.mapper)
            
        else: # a filename was specified
            self.src = vtk.vtkSTLReader()
            self.src.SetFileName(filename)
            self.src.Update()
            self.mapper = vtk.vtkPolyDataMapper()
            self.mapper.SetInput(self.src.GetOutput())
            self.SetMapper(self.mapper)

        self.SetColor(color)
        # SetScaleFactor(double)
        # GetOrigin

class PointCloud(CamvtkActor):
    def __init__(self, pointlist=[]):
        points = vtk.vtkPoints()
        cellArr = vtk.vtkCellArray()
        #Colors = vtk.vtkUnsignedCharArray()
        #Colors.SetNumberOfComponents(3)
        #Colors.SetName("Colors")
        self.zheight = 0
        
        
        n=0
        for p in pointlist:
            vert = vtk.vtkVertex()
            points.InsertNextPoint(p.x, p.y, self.zheight)
            vert.GetPointIds().SetId(0,n)
            cellArr.InsertNextCell( vert )
            #col = clColor(p.cc())
            #Colors.InsertNextTuple3( float(255)*col[0], float(255)*col[1], float(255)*col[2] )
            n=n+1
        
        
        polydata= vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetVerts( cellArr )
        #polydata.GetPointData().SetScalars(Colors)

        polydata.Modified()
        polydata.Update()
        self.src=polydata
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src)
        self.SetMapper(self.mapper)
        #self.SetColor(color)
        

class CLPointCloud(CamvtkActor):
    def __init__(self, pointlist=[]):
        points = vtk.vtkPoints()
        cellArr = vtk.vtkCellArray()
        Colors = vtk.vtkUnsignedCharArray()
        Colors.SetNumberOfComponents(3)
        Colors.SetName("Colors")
        
        n=0
        for p in pointlist:
            vert = vtk.vtkVertex()
            points.InsertNextPoint(p.x, p.y, p.z)
            vert.GetPointIds().SetId(0,n)
            cellArr.InsertNextCell( vert )
            col = clColor(p.cc())
            Colors.InsertNextTuple3( float(255)*col[0], float(255)*col[1], float(255)*col[2] )
            n=n+1
            
        polydata= vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.SetVerts( cellArr )
        polydata.GetPointData().SetScalars(Colors)

        polydata.Modified()
        polydata.Update()
        self.src=polydata
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src)
        self.SetMapper(self.mapper)
        #self.SetColor(color)



class Plane(CamvtkActor):
    def __init__(self, center=(0,0,0), color=(0,0,1) ):
        self.src = vtk.vtkPlaneSource()
        #self.src.SetCenter(center)
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.src.GetOutput())
        self.SetMapper(self.mapper)

        self.SetColor(color)
        self.SetOrigin(center)
        # SetScaleFactor(double)
        # GetOrigin
     


# TODO:
# vtkArcSource
# vtkDiskSource
# vtkFrustumSource
# vtkOutlineSource
# vtkParametricFunctionSource
# PlatonicSolid
# ProgrammableSource (?)
# PSphereSource
# RegularPolygon

#----------------------------------------------------------------

#---- misc helper functions
def vtkPolyData2OCLSTL(vtkPolyData,oclSTL):
    """ read vtkPolyData and add each triangle to an ocl.STLSurf """
    for cellId in range(0,vtkPolyData.GetNumberOfCells()):
        cell = vtkPolyData.GetCell(cellId)
        points = cell.GetPoints()
        plist = []
        for pointId in range(0,points.GetNumberOfPoints()):
            vertex = points.GetPoint(pointId)
            p = ocl.Point(vertex[0],vertex[1],vertex[2])
            plist.append(p)
        t = ocl.Triangle(plist[0],plist[1],plist[2])
        oclSTL.addTriangle(t)
