import openvoronoi as ovd
import ovdvtk
import time
import vtk
import datetime
import math
import random
import os

def drawVertex(myscreen, p, vertexColor, rad=1):
    myscreen.addActor( ovdvtk.Sphere( center=(p.x,p.y,p.z), radius=rad, color=vertexColor ) )

def drawEdge(myscreen, e, edgeColor=ovdvtk.yellow):
    p1 = e[0]
    p2 = e[1]
    myscreen.addActor( ovdvtk.Line( p1=( p1.x,p1.y,p1.z), p2=(p2.x,p2.y,p2.z), color=edgeColor ) )

def drawFarCircle(myscreen, r, circleColor):
    myscreen.addActor( ovdvtk.Circle( center=(0,0,0), radius=r, color=circleColor ) )

def drawDiagram( myscreen, vd ):
    drawFarCircle(myscreen, vd.getFarRadius(), ovdvtk.pink)
    
    for v in vd.getGenerators():
        drawVertex(myscreen, v, ovdvtk.green, 2)
    for v in vd.getVoronoiVertices():
        drawVertex(myscreen, v, ovdvtk.red, 1)
    for v in vd.getFarVoronoiVertices():
        drawVertex(myscreen, v, ovdvtk.pink, 10)
    vde = vd.getVoronoiEdges()
    
    print " got ",len(vde)," Voronoi edges"
    for e in vde:
        drawEdge(myscreen,e, ovdvtk.cyan)

class VD:
    def __init__(self, myscreen, vd, scale=1):
        self.myscreen = myscreen
        self.gen_pts=[ovd.Point(0,0)]
        self.generators = ovdvtk.PointCloud(pointlist=self.gen_pts)
        self.verts=[]
        self.far=[]
        self.edges =[]
        self.generatorColor = ovdvtk.green
        self.vertexColor = ovdvtk.red
        self.edgeColor = ovdvtk.cyan
        self.vdtext  = ovdvtk.Text()
        self.vertexRadius = scale/100
        self.vdtext.SetPos( (50, myscreen.height-150) )
        self.Ngen = 0
        self.vdtext_text = ""
        self.setVDText(vd)
        self.scale=scale
        
        myscreen.addActor(self.vdtext)
        
    def setVDText(self, vd):
        self.Ngen = len( vd.getGenerators() )-3
        self.vdtext_text = "VD with " + str(self.Ngen) + " generators.\n"
        self.vdtext_text += "YELLOW = New point-generator/site\n"
        self.vdtext_text += "PINK = Seed vertex\n"
        self.vdtext_text += "RED = Delete vertices/edges\n"
        self.vdtext_text += "GREEN = Modified VD edges\n"
        self.vdtext.SetText( self.vdtext_text )
        
        
    def setGenerators(self, vd):
        if len(self.gen_pts)>0:
            myscreen.removeActor( self.generators ) 
        #self.generators=[]
        self.gen_pts = []
        for p in vd.getGenerators():
            self.gen_pts.append(self.scale*p)
        self.generators= ovdvtk.PointCloud(pointlist=self.gen_pts) 
        self.generators.SetPoints()
        myscreen.addActor(self.generators)
        self.setVDText(vd)
        myscreen.render() 
    
    def setFar(self, vd):
        for p in vd.getFarVoronoiVertices():
            p=self.scale*p
            myscreen.addActor( ovdvtk.Sphere( center=(p.x,p.y, 0), radius=self.vertexRadius, color=ovdvtk.pink ) )
        myscreen.render() 
    
    def setVertices(self, vd):
        for p in self.verts:
            myscreen.removeActor(p)
        self.verts = []
        for pt in vd.getVoronoiVertices():
            p=self.scale*pt
            actor = ovdvtk.Sphere( center=(p.x,p.y, 0), radius=self.vertexRadius, color=self.vertexColor )
            self.verts.append(actor)
            myscreen.addActor( actor )
            #draw clearance-disk
            """
            cir_actor = ovdvtk.Circle( center=(p.x,p.y,p.z), radius=math.sqrt(pt[1])*self.scale, color=self.vertexColor )
            self.verts.append(cir_actor)
            myscreen.addActor(cir_actor)
            """
        myscreen.render() 
        
    def setEdgesPolydata(self, vd):
        self.edges = []
        self.edges = vd.getEdgesGenerators()
        self.epts = vtk.vtkPoints()
        nid = 0
        lines=vtk.vtkCellArray()
        for e in self.edges:
            p1 = self.scale*e[0]
            p2 = self.scale*e[1] 
            self.epts.InsertNextPoint( p1.x, p1.y, p1.z)
            self.epts.InsertNextPoint( p2.x, p2.y, p2.z)
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0,nid)
            line.GetPointIds().SetId(1,nid+1)
            nid = nid+2
            lines.InsertNextCell(line)
        
        linePolyData = vtk.vtkPolyData()
        linePolyData.SetPoints(self.epts)
        linePolyData.SetLines(lines)
        
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(linePolyData)
        
        self.edge_actor = vtk.vtkActor()
        self.edge_actor.SetMapper(mapper)
        self.edge_actor.GetProperty().SetColor( ovdvtk.cyan )
        myscreen.addActor( self.edge_actor )
        myscreen.render() 

    def setEdges(self, vd):
        for e in self.edges:
            myscreen.removeActor(e)
        self.edges = []
        for e in vd.getEdgesGenerators():
            p1 = self.scale*e[0]  
            p2 = self.scale*e[1] 
            actor = ovdvtk.Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=self.edgeColor )
            myscreen.addActor(actor)
            self.edges.append(actor)
        myscreen.render() 
        
    def setAll(self, vd):
        self.setGenerators(vd)
        #self.setFar(vd)
        self.setVertices(vd)
        self.setEdges(vd)


def writeFrame( w2if, lwr, n ):
    w2if.Modified() 
    current_dir = os.getcwd()
    filename = current_dir + "/frames/vd500_zoomout"+ ('%05d' % n)+".png"
    lwr.SetFileName( filename )
    #lwr.Write()

def regularGridGenerators(far, Nmax):
    # REGULAR GRID
    rows = int(math.sqrt(Nmax))
    print "rows= ",rows
    gpos=[-0.7*far ,  1.4*far/float(rows-1) ]  # start, stride
    plist = []
    for n in range(rows):
        for m in range(rows):
            x=gpos[0]+gpos[1]*n
            y=gpos[0]+gpos[1]*m
            # rotation
            alfa = 0
            xt=x
            yt=y
            x = xt*math.cos(alfa)-yt*math.sin(alfa)
            y = xt*math.sin(alfa)+yt*math.cos(alfa)
            plist.append( ovd.Point(x,y) )
    random.shuffle(plist)
    return plist

def randomGenerators(far, Nmax):
    pradius = (1.0/math.sqrt(2))*far
    plist=[]
    for n in range(Nmax):
        x=-pradius+2*pradius*random.random()
        y=-pradius+2*pradius*random.random()
        plist.append( ovd.Point(x,y) )
    return plist
    
def circleGenerators(far, Nmax):
    # POINTS ON A CIRCLE
    #"""
    #cpos=[50,50]
    #npts = 100
    dalfa= float(2*math.pi)/float(Nmax-1)
    #dgamma= 10*2*math.pi/npts
    #alfa=0
    #ofs=10
    plist=[]
    radius=0.81234*float(far)
    for n in range(Nmax):
        x=float(radius)*math.cos(float(n)*float(dalfa))
        y=float(radius)*math.sin(float(n)*float(dalfa))
        plist.append( ovd.Point(x,y) )
    #random.shuffle(plist)
    return plist
    
    
if __name__ == "__main__":  
    #print ocl.revision()
    myscreen = ovdvtk.VTKScreen()
    ovdvtk.drawOCLtext(myscreen)
    
    w2if = vtk.vtkWindowToImageFilter()
    w2if.SetInput(myscreen.renWin)
    lwr = vtk.vtkPNGWriter()
    lwr.SetInput( w2if.GetOutput() )
    #w2if.Modified()
    #lwr.SetFileName("tux1.png")
    
    scale=1
    myscreen.render()
    random.seed(42)
    far = 1
    
    camPos = far
    zmult = 5
    myscreen.camera.SetPosition(camPos/float(1000), camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    vd = ovd.VoronoiDiagram(far,1200)
    
    vod = VD(myscreen,vd,scale)
    #vod.setAll(vd)
    drawFarCircle(myscreen, scale*vd.getFarRadius(), ovdvtk.orange)
    
    Nmax = 2000
    
    #plist = randomGenerators(far, Nmax)
    plist = regularGridGenerators(far, Nmax)
    #plist = circleGenerators(far, Nmax)
    
    #print plist[169]
    #exit()
    n=1
    t_before = time.time() 
    #delay = 1.5 # 0.533
    delay = 0.1 # 0.533
    #ren = [1,2,3,4,5,59,60,61,62]
    #ren = [16,17]
    ren = range(0,Nmax,100)
    ren = range(0,Nmax,200)
    #ren = [0,1,Nmax-2, Nmax-1, Nmax]
    nf=0
    for p in plist: #[0:20]:
        print "**********"
        print "PYTHON: adding generator: ",n," at ",p
        
        if n in ren:
            vod.setAll(vd)
            myscreen.render()
            time.sleep(delay)
            writeFrame( w2if, lwr, nf )
            nf=nf+1
        
        #GENERATOR
        #"""
        vertexRadius = float(far)/float(100)
        gp=scale*p
        gen_actor = ovdvtk.Sphere( center=(gp.x,gp.y, 0), radius=vertexRadius, color=ovdvtk.yellow )
        if n in ren:
            myscreen.addActor(gen_actor)
            myscreen.render()
            time.sleep(delay)
            writeFrame( w2if, lwr, nf )
            nf=nf+1
        #"""
            
        #SEED
        #"""
        sv = scale*vd.getSeedVertex(p)
        print " seed vertex is ",sv
        seed_actor = ovdvtk.Sphere( center=(sv.x,sv.y, 0), radius=vertexRadius, color=ovdvtk.pink )
        if n in ren:
            myscreen.addActor(seed_actor)
            myscreen.render()
            time.sleep(delay)
            writeFrame( w2if, lwr, nf )
            nf=nf+1
        #"""

        #DELETE-SET
        #"""
        delset = vd.getDeleteSet(p)
        #print " seed vertex is ",sv
        p_actors = []
        if n in ren:
            for pd in delset:
                pos = scale*pd[0]
                type = pd[1]
                p_actor = ovdvtk.Sphere( center=(pos.x,pos.y, 0), radius=vertexRadius, color=ovdvtk.red )
                p_actors.append(p_actor)
            for a in p_actors:
                myscreen.addActor(a)
            myscreen.render()
            time.sleep(delay)
            writeFrame( w2if, lwr, nf )
            nf=nf+1
        #"""
        
        #DELETE-EDGES
        #"""
        delEdges = vd.getDeleteEdges(p)
        modEdges = vd.getModEdges(p)
        #print " seed vertex is ",sv
        edge_actors = []
        if n in ren:
            for e in delEdges:
                p1 = scale*e[0]
                p2 = scale*e[1]
                e_actor = ovdvtk.Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ovdvtk.red ) 
                edge_actors.append(e_actor)
            for e in modEdges:
                p1 = scale*e[0]
                p2 = scale*e[1]
                e_actor = ovdvtk.Line( p1=( p1.x,p1.y, 0), p2=(p2.x,p2.y, 0), color=ovdvtk.green ) 
                edge_actors.append(e_actor)
            for a in edge_actors:
                myscreen.addActor(a)
            myscreen.render()
            time.sleep(delay)
            writeFrame( w2if, lwr, nf )
            nf=nf+1
        #"""
        

        vd.addVertexSite( p )
        
        #w2if.Modified() 
        #lwr.SetFileName("frames/vd500_"+ ('%05d' % n)+".png")
        #lwr.Write()
        
        if n in ren:
            vod.setAll(vd)
            myscreen.render()
            time.sleep(delay)
            writeFrame( w2if, lwr, nf )
            nf=nf+1
        
        if n in ren:
            myscreen.removeActor(gen_actor)
            myscreen.removeActor(seed_actor)
            for a in p_actors:
                myscreen.removeActor(a)
            for a in edge_actors:
                myscreen.removeActor(a)
        
        n=n+1
        #print "**********"
        
    t_after = time.time()
    calctime = t_after-t_before
    print " VD done in ", calctime," s, ", calctime/Nmax," s per generator"
        
    print "PYTHON All DONE."


    
    myscreen.render()    
    myscreen.iren.Start()
