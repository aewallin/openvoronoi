
// OpenVoronoi documentation
// asymptote source file for 
// point-point edge figure

import geometry;

// two point-sites
pair p1=(-0.5,0.5);
pair p2=(0.5,-0.1);

// the apex point
pair pa = 0.5*(p1+p2);

// dx = +1
// dy = -0.6
// so normal is
// dx = +0.6
// dy = 1
size(400,200);

pair n1 = (+0.6, 1);
pair pos_edge = pa+0.8*n1;
pair neg_edge = pa-0.8*n1;

picture pic1;
real size1=180;
size(pic1,size1);
label(pic1,"a)" ,(-0.6,1) ,   black);


draw(pic1,interp(p1,p2,-0.3)--interp(p1,p2,1.3),dashed); 

draw(pic1,  pa--pos_edge ); 
draw(pic1,  pa--neg_edge ); 

perpendicular(pic1, pa,SE,pa--pos_edge,blue);

pair ep = pa + 0.5*n1;
dot(pic1,"$e_{PP}^{+}(t)$" ,ep ,   black);

draw(pic1,  p2--ep, blue, Arrows, PenMargin(5,5) ); 

label(pic1, "$t$", p2--ep ,E , blue);
// label(pic1, "$t_{min}$", p2--pa , S);

// the t_min arrow
real ofs=0.05;
draw(pic1, p2-ofs*n1--pa-ofs*n1,blue,Arrows,PenMargins);
label(pic1, "$t_{min}$", p2-ofs*n1--pa-ofs*n1 , SW, blue);

// the sqrt arrow
pair t1 = (-1,+0.6);
draw(pic1, pa+ofs*t1--ep+ofs*t1,blue,Arrows,PenMargins);
label(pic1, "${\sqrt{t^2-t_{min}^2}}$", pa+ofs*t1--ep+ofs*t1 , NW, blue);


dot(pic1,"$p_1$" ,p1 , S,  red);
dot(pic1,"$p_2$" ,p2 ,NE,red);
dot(pic1,"$p_A$" ,pa , W,  black);


// second picture
picture pic2;
real size2=180;
size(pic2,size2);

label(pic2,"b)" ,(-0.6,1) ,   black);
draw(pic2,interp(p1,p2,-0.3)--interp(p1,p2,1.3),dashed); 

draw(pic2, "$e_{PP}^{+}(t)$" , pa--pos_edge ); 
draw(pic2, "$e_{PP}^{-}(t)$" , pa--neg_edge ); 

perpendicular(pic2, pa,NE,pa--pos_edge,blue);

dot(pic2,"$p_1$" ,p1 , NE,  red);
dot(pic2,"$p_2$" ,p2 ,NE,red);
dot(pic2,"$p_A$" ,pa ,   black);

add( pic1.fit(),(0,0), W);
add( pic2.fit(),(200,0), E);
