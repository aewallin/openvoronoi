import vtk
import math
import ovdvtk
import pocket_sim
import openvoronoi as ovd

w=1920
h=1080
    
myscreen = ovdvtk.VTKScreen(width=w, height=h) 
camPos = 3
zmult=1
myscreen.camera.SetFocalPoint(0, 0, 0)
myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)

ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
myscreen.render()   

ps = pocket_sim.PocketSim(myscreen)
ps.drawFarCircle()
ps.penUp()
ps.penDown()

ps.setTopLeftText("mytext")
ps.Cut()

"""
vol = vtk.vtkImageData()
vol.SetDimensions(512,512,1)
vol.SetSpacing(float(1)/float(512),float(1)/float(512),float(1)/float(512))
vol.SetOrigin(-0.5,-0.5,-0.01)
vol.AllocateScalars()
vol.SetNumberOfScalarComponents(3)
vol.SetScalarTypeToUnsignedChar()

scalars = vtk.vtkCharArray()
red = [255,0,0]
blue= [0,0,255]
green= [0,255,0]
black= [0,0,0]
for n in range(512):
    for m in range(512):
        x = 512/2 -n
        y = 512/2 -m
        d = math.sqrt(x*x+y*y)
        if ( d < 100 ):
            col = red
        elif (d< 200 ):
            col = green
        else:
            col = blue
        for c in range(3):
            scalars.InsertTuple1( n*(512*3) + m*3 +c, col[c] )
 
vol.GetPointData().SetScalars(scalars)
vol.Update()
"""


"""
for n in range(512):
    for m in range(512):
        x = 512/2 -n
        y = 512/2 -m
        d = math.sqrt(x*x+y*y)
        if ( d < 55 ):
            col = black
            for c in range(3):
                scalars.SetValue( n*(512*3) + m*3 +c, chr(col[c]) )

vol.Update()
"""
"""
ia = vtk.vtkImageActor()
ia.SetInput(vol)
ia.InterpolateOff()

myscreen.addActor(ia)
  

"""
"""
ca = ovdvtk.Circle(center=(0,0,0) , radius=1, color=ovdvtk.magenta, resolution=50 )
myscreen.addActor(ca)

ca = ovdvtk.Circle(center=(0,0,0) , radius=0.123, color=ovdvtk.magenta, resolution=50 )
myscreen.addActor(ca)
"""

myscreen.render()   

print "done."
myscreen.iren.Start()
