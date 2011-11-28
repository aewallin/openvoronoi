
// OpenVoronoi documentation
// asymptote source file for 
// point-point edge figure

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
size(500,200);

pair n1 = (+0.6, 1);
pair pos_edge = pa+0.8*n1;
pair neg_edge = pa-0.8*n1;

picture pic1;
real size1=200;
size(pic1,size1);
label(pic1,"a)" ,(-0.6,1) ,   black);
dot(pic1,"$p_1$" ,p1 , NE,  red);
dot(pic1,"$p_2$" ,p2 ,NE,red);
dot(pic1,"$p_A$" ,pa ,   black);

draw(pic1,interp(p1,p2,-0.3)--interp(p1,p2,1.3),dashed); 

draw(pic1,  pa--pos_edge ); 
draw(pic1,  pa--neg_edge ); 

pair ep = pa + 0.5*n1;
dot(pic1,"$e_{PP}^{+}(t)$" ,ep ,   black);

draw(pic1,  p2--ep ); 

label(pic1, "$t$", p2--ep , E );
label(pic1, "$t_{min}$", p2--pa , S);

// second picture
picture pic2;
real size2=200;
size(pic2,size2);

label(pic2,"b)" ,(-0.6,1) ,   black);
draw(pic2,interp(p1,p2,-0.3)--interp(p1,p2,1.3),dashed); 

draw(pic2, "$e_{PP}^{+}(t)$" , pa--pos_edge ); 
draw(pic2, "$e_{PP}^{-}(t)$" , pa--neg_edge ); 

dot(pic2,"$p_1$" ,p1 , NE,  red);
dot(pic2,"$p_2$" ,p2 ,NE,red);
dot(pic2,"$p_A$" ,pa ,   black);

add( pic1.fit(),(0,0), W);
add( pic2.fit(),(250,0), E);
