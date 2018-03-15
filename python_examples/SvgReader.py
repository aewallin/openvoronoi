from svgpathtools import svg2paths2
import math

def sub(p1, p2):
    return (p1[0] - p2[0], p1[0] - p2[0])
def add(p1, p2):
    return (p1[0] + p2[0], p1[0] + p2[0])
def mul(p, f):
    return (p[0] *f, p[1]*f)
def mag(p):
    return math.sqrt(p[0]*p[0]+p[1]*p[1])
def normalize(v):
    return v.mul(1.0/mag(v))
def tup(p):
    return (p.real, p.imag)
def dist(p1, p2):
    return mag(sub(p1, p2))

class SvgReader:

    def __init__(self, filename, error_threshold = .01, min_length = 0.1):
        self.paths, self.attributes, svg_attributes = svg2paths2(filename)
        self.polys = []
        self.miny = self.minx = float("+inf")
        self.maxy = self.maxx = float("-inf")
        self.error_threshold = error_threshold
        self.min_length = min_length

        #print("svg_attributes[idx]:{}".format(svg_attributes))

    def path_to_poly(self, path):
        """Return a list of points that forms a closed polygon

        TODO:
        If the path is not a line the steps parameter will divide the
        segment into \p steps steps
        """
        pts = []
        for segment in path:
            pts.extend(self.adaptiveAdd(segment, 0, 1))
        if not path.isclosed():
            print("Path is not closed!:")
            p = path[0].point(0)
            self.addPoint(pts, p.real, p.imag)
        return pts

    def adaptiveAdd(self, segment, t0, t1):
        """Adaptively sample segment


        """
        p = tup(segment.point(t0))
        p2 = tup(segment.point(t1))
        v = sub(p2, p)
        d = mag(v)
        if d ==0:
            print("no dist")
            return [p]
        t1_2 = (t0+t1)/2.0
        half = add(p,mul(mul(v, 1.0/d), d*.5))
        half_p = tup(segment.point(t1_2))
        diff = dist(half, half_p)
        error = diff/d
        if d < self.min_length  or error < self.error_threshold:
            self.registerPoint(*p)
            return [p]
        else:
            #print("Error {}".format(int(100*error)))
            self.registerPoint(*half_p)
            first = self.adaptiveAdd(segment, t0, t1_2)
            res = []
            if (first):
                res = first
            second = self.adaptiveAdd(segment, t1_2, t1)
            if second:
                res.extend(second)
            return res


    def registerPoint(self, x,y):
        self.miny = min(y, self.miny)
        self.minx = min(x, self.minx)
        self.maxy = max(y, self.maxy)
        self.maxx = max(x, self.maxx)

    def addPoint(self, ls, x, y):
        self.registerPoint(x,y)
        ls.append((x,y))

    def parse(self):
        """Parse all paths to linear polygons

        Determine AABB, center and store transform
        """
        for path in self.paths:
            for continuous_subpath in path.continuous_subpaths():
                poly = self.path_to_poly(continuous_subpath)
                if poly:
                    self.polys.append(poly)

        self.center = (self.maxx + self.minx/2, (self.maxy + self.miny)/2)
        w = self.maxx - self.minx
        h = self.maxy - self.miny
        self.radius = math.sqrt(w*w + h*h)


    def offsetPoint(self, p):
        return (p[0] - self.offsetx,  p[1] - self.offsety)

    def centerPolys(self):
        self.offsetx = self.center[0]
        self.offsety = self.center[1]

        adjusted = []
        for poly in self.polys:
            ls = []
            for p in poly:
                ls.append(self.offsetPoint(p))
   #             ls.reverse()
            adjusted.append(ls)
        self.polys = adjusted

    def stats(self):
        res = ["{} paths".format(len(self.polys))]
        res += ["miny:{}".format(self.miny)]
        res += ["minx:{}".format(self.minx)]
        res += ["maxy:{}".format(self.maxy)]
        res += ["maxx:{}".format(self.maxx)]
        res += ["center:{}".format(self.center)]
        res += ["radius:{}".format(self.radius)]
        cnt = 0;
        for idx in range(len(self.polys)):
            p = self.polys[idx]
            s = len(p)
            res += ["\t Path{} with {} vertices".format(idx +1, s)]
            cnt += s
        res += ["Total {} vertices".format(cnt)]
        return "\n".join(res);

    def __getitem__(self, index):
        return self.polys[index]
