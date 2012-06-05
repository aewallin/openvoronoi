
clearance_height= 20
feed_height = 10
feed = 200
plunge_feed = 100
metric = True
scale = 1

def line_to(x,y,z):
    print "G1 X% 8.6f Y% 8.6f Z% 8.6f F%.0f" % (x, y, z, feed)

def xy_line_to(x,y):
    print "G1 X% 8.4f Y% 8.4f " % (scale*x, scale*y)

# (endpoint, radius, center, cw?)
def xy_arc_to( x,y, r, cx,cy, cw ):
    if (cw):
        print "G2 X% 8.5f Y% 8.5f R% 8.5f" % (scale*x, scale*y, r)
    else:
        print "G3 X% 8.5f Y% 8.5f R% 8.5f" % (scale*x, scale*y, r)
    # FIXME: optional IJK format arcs
    
def xy_rapid_to(x,y):
    print "G0 X% 8.4f Y% 8.4f " % (scale*x, scale*y)

def pen_up():
    print "G0Z% 8.4f " % (clearance_height)

def pen_down(z=0):
    print "G0Z% 8.4f" % (feed_height)
    plunge(z)

def plunge(z):
    print "G1 Z% 8.4f F% 8.0f" % (z, plunge_feed)

def preamble():
    if (metric):
        print "G21 F% 8.0f" % (feed) # G20 F6 for inch
    else:
        print "G20 F% 8.0f" % (feed) # G20 F6 for inch
        
    print "G64 P0.001"
    pen_up()
    print "G0 X0 Y0"

def postamble():
    pen_up()
    print "M2"

if __name__ == "__main__":
    print "Nothing to see here."
