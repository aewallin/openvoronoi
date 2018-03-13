import openvoronoi as ovd
import math


def regular_n_polygon(n=3, verbose=True):
    far = 1.0
    vd = ovd.VoronoiDiagram(far, 120)

    r = 0.9
    # input points (vertices/sites)
    if verbose:
        print "regular polygon n=%d" % (n)
    id_list = []
    for i in range(n):
        a = math.pi * i / n * 2.0
        x = r * math.cos(a)
        y = r * math.sin(a)
        id_list.append(vd.addVertexSite(ovd.Point(x, y)))
        if verbose:
            print "%d angle=%1.3f point=(%1.3f, %1.3f) id=%d" % (i, a, x, y, id_list[-1])

    for i in range(n):
        if verbose:
            print "add line", i, ":", id_list[i - 1], id_list[i]
        # on first iteration negative index, but that's ok
        vd.addLineSite(id_list[i - 1], id_list[i])
    check = vd.check()
    if verbose:
        print "VD check: ", check
    return check


if __name__ == "__main__":
    regular_n_polygon(3)  # OK
    regular_n_polygon(4)  # warnings/errors when adding line 3
    # regular_n_polygon(5) # OK
    # regular_n_polygon(6) # warnings/errors when adding line 4
    # regular_n_polygon(7) # OK
    # regular_n_polygon(8) # warnings/errors when adding line 5
    # regular_n_polygon(9) # OK
    # regular_n_polygon(10) # warnings/errors when adding line 6
    # regular_n_polygon(11) # OK
    # regular_n_polygon(12) # OK (? no warning for line 7)
    # regular_n_polygon(13) # OK
    # regular_n_polygon(14) # OK (? no warning for line 8)
    # regular_n_polygon(15) # OK
    # regular_n_polygon(16) # warnings/errors when adding line 9
    # regular_n_polygon(17) # OK
    # regular_n_polygon(18) # OK (? no warning for line 10?)
    # regular_n_polygon(19) # OK
    # regular_n_polygon(20) # warnings/errors when adding line 11
