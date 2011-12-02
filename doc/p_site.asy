
// OpenVoronoi documentation
// asymptote source file for 
// point-point edge figure

import geometry;
size(200,200);

// two point-sites
pair p0=(0,0);

// dx = +1
// dy = -0.6
// so normal is
// dx = +0.6
// dy = 1



dot("$p$" ,p0 , NE,  black);
path g = unitcircle;
draw(unitcircle,  dashed+red );
label("$O_P^+(t)$" ,unitcircle ,   NE,red);

path radiusarrow = (0,0)--(-sqrt(2)/2,sqrt(2)/2);
draw( radiusarrow,blue,Arrows,PenMargin(5,5));
label( "$t$", radiusarrow , NE, blue);
