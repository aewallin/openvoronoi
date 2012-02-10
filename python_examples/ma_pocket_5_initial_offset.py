import openvoronoi as ovd
import ovdvtk

import time
import vtk
import datetime
import math
import random
import os
import sys
import pickle
import gzip

def drawCircle(myscreen, c, r, circlecolor):
    ca = ovdvtk.Circle(center=(c.x,c.y,0) , radius=r, color=circlecolor, resolution=50 )
    myscreen.addActor(ca)
    
def drawPoint( myscreen, c, pcolor , rad = 0.002):
    ca = ovdvtk.Sphere(center=(c.x,c.y,0) , radius=rad, color=pcolor)
    myscreen.addActor(ca)


# rotate by cos/sin. from emc2 gcodemodule.cc
def rotate(x, y,  c,  s):
    tx = x * c - y * s;
    y = x * s + y * c;
    x = tx;
    return [x,y]
    
def drawArc(myscreen, pt1, pt2, r, cen,cw,arcColor):
    # draw arc as many line-segments
    start = pt1-cen
    end = pt2-cen
    theta1 = math.atan2(start.x,start.y)
    theta2 = math.atan2(end.x,end.y)
    alfa=[] # the list of angles
    da=0.1
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
    dlength = min(0.01, arclength/10)
    steps = int( float(arclength) / float(dlength))
    rsteps = float(1)/float(steps)
    dc = math.cos(-dtheta*rsteps) # delta-cos  
    ds = math.sin(-dtheta*rsteps) # delta-sin
    
    previous = pt1
    tr = [start.x, start.y]
    for i in range(steps):
        tr = rotate(tr[0], tr[1], dc, ds) #; // rotate center-start vector by a small amount
        x = cen.x + tr[0] 
        y = cen.y + tr[1] 
        current = ovd.Point(x,y)
        myscreen.addActor( ovdvtk.Line(p1=(previous.x,previous.y,0),p2=(current.x,current.y,0),color=arcColor) )
        previous = current 

def rapid_to_next(myscreen, prv_tang, nxt_tang, c1, r1, c2, r2, prv, nxt):
    # rapid from prev, to nxt
    # while staying inside c1(r1) and c2(r)

    rad_default = 0.03
    rad = min( rad_default, 0.9*r1 , 0.9*r2)
    
    prv_tang.normalize()
    nxt_tang.normalize()

    prv_normal = -1*prv_tang.xy_perp()
    nxt_normal = nxt_tang.xy_perp()
    
    cen1 = prv + rad*prv_normal # + rad1*prv_tang
    cen2 = nxt - rad* nxt_normal #rapid_tang # + rad1*prv_tang

    rapid_tang = cen2-cen1
    rapid_tang.normalize()

    trg1 = cen1 + rad* rapid_tang.xy_perp() #prv_tang
    src2 = cen2 + rad* rapid_tang.xy_perp() 

    drawArc(myscreen, prv, trg1, rad, cen1, True, ovdvtk.blue) # lead-out arc
    ovdvtk.drawLine(myscreen, trg1, src2, ovdvtk.cyan) # rapid
    drawArc(myscreen, src2, nxt, rad, cen2, True, ovdvtk.lblue) # lead-in arc

def rapid_to_new_branch(myscreen, prv_tang, nxt_tang, c1, r1, c2, r2, prv, nxt):
    # rapid from prev, to nxt
    # while staying inside c1(r1) and c2(r)
    rad_default = 0.03
    rad1 = min( rad_default, 0.9*r1 ) # wrong? we get the new-branch r1 here, while we would want the old-branch r1
    rad2 = min( rad_default, 0.9*r2 )    
    prv_tang.normalize()
    nxt_tang.normalize()

    prv_normal = -1*prv_tang.xy_perp()
    nxt_normal = nxt_tang.xy_perp()
    
    cen1 = prv + rad1*prv_normal # + rad1*prv_tang
    
    cen2 = nxt - rad2* nxt_normal #rapid_tang # + rad1*prv_tang

    rapid_tang = cen2-cen1
    rapid_tang.normalize()
    
    trg1 = cen1 + rad1* prv_tang 
    src2 = cen2 - rad2* nxt_tang 

    drawArc(myscreen, prv, trg1, rad1, cen1, True, ovdvtk.orange) # lead-out arc
    ovdvtk.drawLine(myscreen, trg1, src2, ovdvtk.magenta) # rapid
    drawArc(myscreen, src2, nxt, rad2, cen2, True, ovdvtk.mag2) # lead-in arc

def final_lead_out(myscreen, prv_tang, nxt_tang, c1, r1, c2, r2, prv, nxt):
    # rapid from prev, to nxt
    # while staying inside c1(r1) and c2(r)
    rad_default = 0.03
    rad1 = min( rad_default, 0.9*r1 ) # wrong? we get the new-branch r1 here, while we would want the old-branch r1
    #rad2 = min( rad_default, r2 )    
    prv_tang.normalize()
    #nxt_tang.normalize()

    prv_normal = -1*prv_tang.xy_perp()
    #nxt_normal = nxt_tang.xy_perp()
    
    cen1 = prv + rad1*prv_normal # + rad1*prv_tang
    
    #cen2 = nxt - rad2* nxt_normal #rapid_tang # + rad1*prv_tang

    #rapid_tang = cen2-cen1
    #rapid_tang.normalize()
    
    trg1 = cen1 + rad1* prv_tang 
    #src2 = cen2 - rad2* nxt_tang 

    drawArc(myscreen, prv, trg1, rad1, cen1, True, ovdvtk.red) # lead-out arc
    #ovdvtk.drawLine(myscreen, trg1, src2, ovdvtk.magenta) # rapid
    #drawArc(myscreen, src2, nxt, rad2, cen2, True, ovdvtk.mag2) # lead-in arc

def spiral_clear(myscreen, out_tangent, in_tangent, c1, r1, c2, r2, out1, in1):
    print "spiral clear!"
    # end spiral at in1
    # archimedean spiral
    # r = a + b theta
    
    in1_dir = in1-c1
    in1_theta = math.atan2(in1_dir.y,in1_dir.x)
    #in1_theta = in1_theta
    print "c1 =", c1
    print "in1 = ",in1
    print " end theta = ",in1_theta
    drawPoint( myscreen, c1, ovdvtk.red )
    drawPoint( myscreen, in1, ovdvtk.blue, 0.006 )
    # width = 2*pi*b
    # => b = width/(2*pi)
    b=0.01/(2*math.pi)
    # r = a + b in1_theta = r_max
    # =>
    # a = r_max-b*in1_theta
    a = r1 - b*in1_theta
    
    # figure out the start-angle
    theta_min = in1_theta
    theta_max = in1_theta
    dtheta = 0.1
    min_r = 0.001
    while True:
        r = a+b*theta_min
        if r < min_r:
            break
        else:
            theta_min = theta_min - dtheta
    print "start_theta = ", theta_min
    
    Npts = (theta_max - theta_min)/dtheta
    Npts = int(Npts)
    print "spiral has ",Npts," points"
    p = ovd.Point(c1)
    for n in range(Npts+1):
        theta = theta_min + n*dtheta
        r = a + b*theta
        theta = theta - 2* abs(in1_theta - math.pi/2 )
        trg = c1 + r*ovd.Point( -math.cos(theta), math.sin(theta) )
        ovdvtk.drawLine(myscreen, p,trg,ovdvtk.pink)
        p = trg
    
if __name__ == "__main__":  
    #w=2500
    #h=1500

    w=1920
    h=1080
    #w=1024
    #h=1024
    myscreen = ovdvtk.VTKScreen(width=w, height=h) 
    ovdvtk.drawOCLtext(myscreen, rev_text=ovd.version() )
    
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
    zmult = 1.8
    # camPos/float(1000)
    myscreen.camera.SetPosition(0, -camPos/float(1000), zmult*camPos) 
    myscreen.camera.SetClippingRange(-(zmult+1)*camPos,(zmult+1)*camPos)
    myscreen.camera.SetFocalPoint(0.0, 0.22, 0)
    
    vd = ovd.VoronoiDiagram(far,120)
    print ovd.version()
    
    # for vtk visualization
    vod = ovdvtk.VD(myscreen,vd,float(scale), textscale=0.01, vertexradius=0.003)
    vod.drawFarCircle()
    
    vod.textScale = 0.02
    vod.vertexRadius = 0.0031
    vod.drawVertices=0
    vod.drawVertexIndex=0
    vod.drawGenerators=0
    vod.offsetEdges = 0
    vd.setEdgeOffset(0.05)
    
    linesegs = 1 # switch to turn on/off line-segments
    
    segs = []
    #ovd.Point(1,1)
    eps=0.9
    p1=ovd.Point(-0.1,-0.2)
    p2=ovd.Point(0.2,0.1)
    p3=ovd.Point(0.4,0.2)
    p4=ovd.Point(0.6,0.6)
    p5=ovd.Point(-0.6,0.3)

    pts = [p1,p2,p3,p4,p5]
    
    #t_after = time.time()
    #print ".done in {0:.3f} s.".format( t_after-t_before )
    times=[]
    id_list = []
    m=0
    t_before = time.time()
    for p in pts:
        
        id_list.append( vd.addVertexSite( p ) )
        #print m," added vertex", seg_id[0]
        m=m+1
   
    t_after = time.time()
    times.append( t_after-t_before )
    
    print "all point sites inserted. "
    vd.check()
    
    t_before = time.time()
    vd.addLineSite( id_list[0], id_list[1])
    vd.addLineSite( id_list[1], id_list[2])
    vd.addLineSite( id_list[2], id_list[3])
    vd.addLineSite( id_list[3], id_list[4])
    vd.addLineSite( id_list[4], id_list[0])
    t_after = time.time()
    times.append( t_after-t_before )
    vd.check()
    
    pi = ovd.PolygonInterior(True)
    vd.filter_graph(pi)
    ma = ovd.MedialAxis()
    vd.filter_graph(ma)
    
    vod.setVDText2(times)
    
    vod.setAll()
    
    mapocket = ovd.MedialAxisPocket(vd.getGraph())
    mapocket.setWidth(0.01)
    
    maxmic = mapocket.maxMic()
    
    #print maxmic
    previous_center = maxmic[0]
    previous_radius = maxmic[1]
    cl = ovd.Point(0,0)
    
    # the initial largest MIC. to be cleared with a spiral-path
    drawCircle( myscreen, maxmic[0], maxmic[1] , ovdvtk.red )
    
    # the rest of the MICs are then cleared
    nframe=0
    first = True
    previous_out1 = ovd.Point()
    out_tangent = ovd.Point()
    in_tangent = ovd.Point()
    while True:
        mic = mapocket.nxtMic()
        if 0: #nframe == 5:
            break
        if len(mic) >= 2:
            cen2 = mic[0]
            r2 = mic[1]
            #drawCircle( myscreen, mic[0], mic[1] , ovdvtk.green )
            previous_center = mic[6]
            previous_radius = mic[7]
            new_branch = mic[8] # true/false indicates if we are starting on new branch
            prev_branch_center = mic[9]
            prev_branch_radius = mic[10] # old branch MIC radius

            # these are the bi-tangent points
            if not mic[2].is_right( previous_center, cen2 ):
                in1 = mic[2]
                in2 = mic[4]
                out2 = mic[5]
                out1 = mic[3]
            else:
                in1 = mic[3]
                in2 = mic[5]
                out2 = mic[4]
                out1 = mic[2]
                
            drawPoint( myscreen, in1, ovdvtk.red )
            drawPoint( myscreen, out1, ovdvtk.pink )
            drawPoint( myscreen, in2, ovdvtk.green )
            drawPoint( myscreen, out2, ovdvtk.grass )
            in_tangent = in2-in1
            # rapid traverse to in1
            if not first:
                if new_branch:
                    # new branch re-position move
                    rapid_to_new_branch(myscreen, out_tangent, in_tangent, prev_branch_center, prev_branch_radius , cen2, r2, previous_out1, in1)
                else:
                    # normal arc-rapid-arc to next MIC
                    rapid_to_next(myscreen, out_tangent, in_tangent, previous_center, previous_radius, cen2, r2, previous_out1, in1)
            else:
                # spiral-clear the start-MIC. The spiral should end at in1
                spiral_clear(myscreen, out_tangent, in_tangent, previous_center, previous_radius, cen2, r2, previous_out1, in1)
                print "No rapid-move on first-iteration."
                first = False

                
            # in bi-tangent
            ovdvtk.drawLine(myscreen, in1, in2, ovdvtk.green)
            # draw arc
            drawArc(myscreen, in2, out2, r2, cen2, True, ovdvtk.green)
            # out bi-tangent
            ovdvtk.drawLine(myscreen, out2, out1, ovdvtk.green)
            
            previous_out1 = out1 # this is used as the start-point for the rapid on the next iteration
            #previous_center = cen2
            #previous_radius = r2
            out_tangent = out1-out2
            
            #uncomment to write screenshots to disk
            #w2if.Modified()
            #lwr.SetFileName("frames/%06d.png" % ( nframe ) )
            #lwr.Write()
        else:
            # end of operation. do a final lead-out arc.
            final_lead_out(myscreen, out_tangent, in_tangent, previous_center, previous_radius, cen2, r2, previous_out1, in1)
            print "Final lead-out arc"
            break
        nframe = nframe+1
        
    print "mic-pocket done."
    
    

        
    print "PYTHON All DONE."

    myscreen.render()   

     
    myscreen.iren.Start()
