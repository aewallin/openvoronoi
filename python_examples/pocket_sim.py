import vtk
import math
import ovdvtk
import openvoronoi as ovd

class PocketSim:
    def __init__(self, myscreen,  textscale=0.06):
        self.myscreen = myscreen
        self.actors=[]
        self.xsize = 512
        self.ysize = 512
        self.vol = vtk.vtkImageData()
        self.vol.SetDimensions(512,512 ,1)
        #print float(1)/float(self.xsize)
        self.vol.SetSpacing(float(1)/float(self.xsize),float(1)/float(self.ysize),float(1)/float(self.ysize))
        #self.vol.SetSpacing(1,1,1)
        self.vol.SetOrigin(-0.5,-0.5,-0.01)
        self.vol.SetScalarTypeToUnsignedChar()
        self.vol.SetNumberOfScalarComponents(3)
        self.vol.AllocateScalars()
        
        self.scalars = vtk.vtkCharArray()
        # initialize the image
        for n in range(512):
            for m in range(512):
                x = 512/2 -n
                y = 512/2 -m
                d = math.sqrt(x*x+y*y)
                if ( d < 100 ):
                    #print "red"
                    col = ovdvtk.red
                elif (d< 200 ):
                    col = ovdvtk.green
                else:
                    col = ovdvtk.blue
                for c in range(3):
                    self.scalars.InsertTuple1( n*(512*3) + m*3 +c, 255*col[c] )
        
        self.vol.GetPointData().SetScalars( self.scalars )
        self.vol.Update()
        self.ia = vtk.vtkImageActor()
        self.ia.SetInput(self.vol)
        self.ia.InterpolateOff()
        self.vol.Update()
        self.actors.append( self.ia )
        self.myscreen.addActor( self.ia )
        
        """
        self.verts=[]
        self.far=[]
        self.edges =[]
        self.gens =[]
        
        self.pointsiteColor = yellow
        self.generatorColor = yellow
        self.vertexColor = blue
        self.seedColor = pink
        self.edgeColor = cyan
        self.vertexRadius = vertexradius
        """
        self.topleft_text  = ovdvtk.Text()
        self.topleft_text.SetPos( (50, myscreen.height-70) )
        myscreen.addActor(self.topleft_text)
        self.actors.append( self.topleft_text )
        
        self.vdtext2  = ovdvtk.Text()
        self.vdtext2.SetPos( (myscreen.width-500, 50) ) 
        self.vdtext2.SetText( "PocketSim" )
        myscreen.addActor(self.vdtext2)
        self.actors.append( self.vdtext2 )
        
        self.gittext  = ovdvtk.Text()
        self.gittext.SetPos( (50, 50) )
        self.gittext_text = "github.com/aewallin"
        self.gittext.SetText( self.gittext_text )
        myscreen.addActor(self.gittext)
        
        #self.N_pointgen = 0
        #self.N_linegen = 0
        #self.vdtext_text = "CPU Time:"
        #self.setVDText()
        
        #self.drawClearanceDisk = 0
        #self.textScale = textscale
        #self.drawVertexIndex=1
        #self.drawVertices=1
        #self.drawGenerators=1
        #self.offsetEdges = 0
        #self.drawNullEdges = 1
        
        self.pen_up = True
        self.setCutter(0.03)
        self.p = ovd.Point(0,0)
        self.ca = ovdvtk.Circle(center=(self.p.x,self.p.y,0) , radius=self.r, color=ovdvtk.green, resolution=50 )
        self.myscreen.addActor(self.ca)
        
        self.down_color = ovdvtk.dgreen
        self.up_color = ovdvtk.green
        self.p = ovd.Point(0,0)
        
    def setCutter(self,r):
        self.r = r
    
    def penUp(self):
        self.pen_up = True
        self.ca.SetColor( self.up_color )
        self.Update()
        
    def penDown(self):
        self.pen_up = False
        self.ca.SetColor( self.down_color )
        self.Cut()
        self.Update()
        
    def Update(self):
        self.myscreen.render()

    def drawFarCircle(self, circleColor=ovdvtk.magenta):
        self.myscreen.addActor( ovdvtk.Circle( center=(0,0,0), radius=1, color=circleColor ) )
        
    def Cut(self):
        j=0
        for n in range(512):
            for m in range(512):
                x = 512/2 -n
                y = 512/2 -m
                d = math.sqrt(x*x+y*y)
                if ( d < 55 ):
                    #print "black"
                    col = ovdvtk.yellow
                    print n," , ",m ," set black", col
                    j=j+1
                    for c in range(3):
                        self.scalars.SetValue( n*(512*3) + m*3 +c, chr(255*col[c]) )
                        #self.vol.GetPointData().SetValue( n*(512*3) + m*3 +c, chr(col[c]) )
        self.vol.GetPointData().SetScalars( self.scalars )
        self.vol.Update()
        self.Update()
        print j,"cut"
        
    """
    def setVDText(self):
        self.N_pointgen = self.vd.numPointSites()
        self.N_linegen = self.vd.numLineSites()

        #self.vdtext_text = " "
        self.vdtext_text = "Voronoi-Diagram with :\n"
        self.vdtext_text += str(self.N_pointgen) + " point-sites.\n"
        self.vdtext_text += str(self.N_linegen) + " line-sites.\n"
        #self.vdtext_text += "YELLOW = New point-generator/site\n"
        #self.vdtext_text += "PINK = Seed vertex\n"
        #self.vdtext_text += "RED = Delete vertices/edges\n"
        #self.vdtext_text += "GREEN = Modified VD edges\n"
        self.vdtext.SetText( self.vdtext_text )
    """
    def setTopLeftText(self,text):
        self.topleft_text.SetText( text )
        
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
            if status == ovd.VoronoiVertexStatus.IN:
                vcolor = red
            elif status == ovd.VoronoiVertexStatus.OUT:
                vcolor = blue
            elif status == ovd.VoronoiVertexStatus.NEW:
                vcolor = green
                    
            if status == ovd.VoronoiVertexStatus.IN or status == ovd.VoronoiVertexStatus.OUT or status == ovd.VoronoiVertexStatus.NEW:
                if ( vtype != ovd.VoronoiVertexType.SEPPOINT and vtype != ovd.VoronoiVertexType.ENDPOINT and vtype != ovd.VoronoiVertexType.POINTSITE):
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
        
        if ( (src_status == ovd.VoronoiVertexStatus.IN) and (trg_status == ovd.VoronoiVertexStatus.IN)):
            return red
        elif ( ((src_status == ovd.VoronoiVertexStatus.OUT) and 
             (trg_status == ovd.VoronoiVertexStatus.IN) ) or 
             ((src_status == ovd.VoronoiVertexStatus.IN) and 
             (trg_status == ovd.VoronoiVertexStatus.OUT) ) ):
            return purple
        elif ( ((src_status == ovd.VoronoiVertexStatus.OUT) and 
             (trg_status == ovd.VoronoiVertexStatus.NEW) ) or 
             ((src_status == ovd.VoronoiVertexStatus.NEW) and 
             (trg_status == ovd.VoronoiVertexStatus.OUT) ) ):
            return default_color
        elif ( ((src_status == ovd.VoronoiVertexStatus.IN) and 
             (trg_status == ovd.VoronoiVertexStatus.NEW) ) or 
             ((src_status == ovd.VoronoiVertexStatus.NEW) and 
             (trg_status == ovd.VoronoiVertexStatus.IN) ) ):
            return red    
        elif ( (src_status == ovd.VoronoiVertexStatus.NEW) and (trg_status == ovd.VoronoiVertexStatus.NEW)):
            return green
        else:
            return default_color
            
    
    def edgeTypeColor(self, edgeType, src_status, trg_status):
        #print " edgeStatusColor edgeType= ",edgeType
        if (edgeType == ovd.VoronoiEdgeType.LINELINE):
            return self.edgeStatusColor(src_status,trg_status,lblue)
        if (edgeType == ovd.VoronoiEdgeType.PARA_LINELINE):
            return self.edgeStatusColor(src_status,trg_status,blue)
        if (edgeType == ovd.VoronoiEdgeType.LINE):
            return  self.edgeStatusColor(src_status,trg_status,cyan)
        if (edgeType == ovd.VoronoiEdgeType.OUTEDGE):
            return pink
        if (edgeType == ovd.VoronoiEdgeType.SEPARATOR):
            return self.edgeStatusColor(src_status,trg_status, mag2)
        if (edgeType == ovd.VoronoiEdgeType.LINESITE):
            return yellow
        if (edgeType == ovd.VoronoiEdgeType.PARABOLA):
            return self.edgeStatusColor(src_status,trg_status, blue2)
        else:
            #print "UNKOWN edge type = ", edgeType
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
                
            if (etype == ovd.VoronoiEdgeType.LINE):
                ecolor = self.edgeStatusColor(src_status,trg_status,cyan)
                for n in range( len(epts)-1 ):
                    p1 = self.scale*epts[n]  
                    p2 = self.scale*epts[n+1] 
                    actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                    self.myscreen.addActor(actor)
                    self.edges.append(actor)

            if (etype == ovd.VoronoiEdgeType.OUTEDGE):
                ecolor = pink
                p1 = self.scale*epts[0]  
                p2 = self.scale*epts[1] 
                actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                self.myscreen.addActor(actor)
                self.edges.append(actor)
            if (etype == ovd.VoronoiEdgeType.SEPARATOR):
                ecolor = self.edgeStatusColor(src_status,trg_status, orange)
                p1 = self.scale*epts[0]  
                p2 = self.scale*epts[1] 
                actor = Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ecolor)
                self.myscreen.addActor(actor)
                self.edges.append(actor)
            if (etype == ovd.VoronoiEdgeType.LINESITE):
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

        self.myscreen.render() 
        
    def setAll(self):
        self.setVDText()
        self.setGenerators()
        #self.setFar(vd)
        self.setVertices()
        self.setEdgesPolydata()
        #self.setEdges()
 
