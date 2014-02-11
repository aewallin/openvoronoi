import ovdvtk
import math
import vtk

pts = []
# calculate Hilbert curve, using recursion
# based on: http://www.andyshelley.co.uk/axishilbert/index.php
#    Copyright 2010   Andy Shelley <andy@andyshelley.co.uk>
#    GPL v2 or later
# self.Hilbert    (0, self.Blocks, 90, self.Side, self.Feed)
def Hilbert ( ca, level, angle, size, p): 
    if level == 0:
        return

    ca = (ca + 360 - angle) % 360
    Hilbert( ca, level - 1, -angle, size,p )
    step = size
    if (ca == 0):
        p[0]=p[0]+step
    elif (ca == 90):
        p[1]=p[1]-step
    elif (ca == 180):
        p[0]=p[0]-step
    else:
        p[1]=p[1]+step
    pts.append( [ p[0], p[1] ] )

    ca = (ca + angle) % 360
    Hilbert ( ca, level - 1, angle, size, p)
    if (ca == 0):
        p[0]=p[0]+step
    elif (ca == 90):
        p[1]=p[1]-step
    elif (ca == 180):
        p[0]=p[0]-step
    else:
        p[1]=p[1]+step
    pts.append( [ p[0], p[1] ] )
    Hilbert ( ca, level - 1, angle, size, p)
    
    ca = (ca + angle) % 360
    if (ca == 0):
        p[0]=p[0]+step 
    elif (ca == 90):
        p[1]=p[1]-step
    elif (ca == 180):
        p[0]=p[0]-step
    else:
        p[1]=p[1]+step
    pts.append( [ p[0], p[1] ] )
    Hilbert ( ca, level - 1, -angle, size, p)
    ca = (ca + 360 - angle) % 360

def drawCurve(myscreen,curve,loopColor):
    n = 0
    N = len(curve)
    first_point=[]
    previous=[]
    for p in curve:
        if n==0: # don't draw anything on the first iteration
            previous=p 
            first_point = p
        else:
            myscreen.addActor( ovdvtk.Line(p1=(previous[0],previous[1],0),p2=(p[0],p[1],0),color=loopColor) )
            previous=p
        n=n+1
    print "rendered curve  with ", len(curve), " points"

def drawCurve2(myscreen, curve, curvecolor):
    oPoints = vtk.vtkPoints()
    lineCells=vtk.vtkCellArray()
    idx = 0
    last_idx = 0
    segs=[]
    first = 1
    print " curve with ", len(curve)," points"
    for p in curve:
        oPoints.InsertNextPoint( p[0], p[1], 0)
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
        lineCells.InsertNextCell(line)
    
    linePolyData = vtk.vtkPolyData()
    linePolyData.SetPoints(oPoints)
    linePolyData.SetLines(lineCells)
    linePolyData.Modified() 
    linePolyData.Update()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(linePolyData)
    edge_actor = vtk.vtkActor()
    edge_actor.SetMapper(mapper)
    edge_actor.GetProperty().SetColor( curvecolor )
    myscreen.addActor( edge_actor )

if __name__ == "__main__":
    
    w=1024
    h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    #ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    # draw a unit-circle
    ca = ovdvtk.Circle(center=(0,0,0) , radius=1, color=(0,1,1), resolution=50 )
    myscreen.addActor(ca)   
    
    scale=1
    far = 1
    camPos = far
    zmult = 3
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0, 0)
    
    size = 1.1
    
    # bounding-box
    #myscreen.addActor( ovdvtk.Line(p1=(-size/2,-size/2,0),p2=(-size/2,size/2,0),color=ovdvtk.orange) )
    #myscreen.addActor( ovdvtk.Line(p1=(-size/2,size/2,0),p2=(size/2,size/2,0),color=ovdvtk.orange) )
    #myscreen.addActor( ovdvtk.Line(p1=(size/2,size/2,0),p2=(size/2,-size/2,0),color=ovdvtk.orange) )
    #myscreen.addActor( ovdvtk.Line(p1=(size/2,-size/2,0),p2=(-size/2,-size/2,0),color=ovdvtk.orange) )
    
    startpt = [-size/2,-size/2]
    level = 1
    step = size / math.pow(2, level)

    Hilbert( 0, level, 90, step, startpt)
    #print pts
    drawCurve2(myscreen,pts,ovdvtk.yellow)
    
    print "PYTHON All DONE."

    myscreen.render()        
    myscreen.iren.Start()
