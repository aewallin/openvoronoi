
// OpenVoronoi documentation
// asymptote source file for 
// line site figure

import geometry;
size(200,200);


// endpoints for line
pair p1=2*(-0.5,0.5);
pair p2=2*(0.5,-0.1);

// dx = +1
// dy = -0.6
// so normal is
// dx = +0.6
// dy = 1
pair tangent=(1,-0.6);

pair normal=(+0.6,1);
real ofs = 0.4;
pair pos_p1 = p1 + ofs*normal;
pair pos_p2 = p2 + ofs*normal;
pair neg_p1 = p1 - ofs*normal;
pair neg_p2 = p2 - ofs*normal;

real t = 0.1;
pair arrow_p1 = p1+t*tangent;
pair arrow_p2 = pos_p1+t*tangent;
pair arrow_p3 = neg_p1+t*tangent;
path pos_arrow = arrow_p1--arrow_p2;
path neg_arrow = arrow_p1--arrow_p3;
draw( pos_arrow,blue,Arrows,PenMargin(5,5));
draw( neg_arrow,blue,Arrows,PenMargin(5,5));
label( "$t$" , pos_arrow , NW, blue);
label( "$t$" , neg_arrow , NW, blue);

// dot("$p$" ,p1--p2 , NE,  black);

draw(p1--p2); 
label( "$l: ax+by+c=0$" , p1--p2 , NE);

draw(pos_p1--pos_p2,dashed+red); 
draw(neg_p1--neg_p2,dashed+red); 
label( "$O_L^+(t): ax+by+c+t=0$" , pos_p1--pos_p2 , NE,red);
label( "$O_L^-(t): ax+by+c-t=0$" , neg_p1--neg_p2 , SW,red);
