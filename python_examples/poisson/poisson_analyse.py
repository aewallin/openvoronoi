#import openvoronoi as ovd
#import ovdvtk
#import time
#import vtk
#import datetime
import math
#import random
#import os
import gzip
import pickle


import numpy as np
import numpy.numarray as na
import pylab as P

def sortedDict(adict):
    keys = adict.keys()
    keys.sort()
    l = []
    for key in keys:
        d=[]
        d.append(key)
        d.append(adict[key])
        l.append(d)
    return l

pn = []
pn.append( [3, 1.124e-2, 0.000021e-2] )
pn.append( [4, 1.068454e-1, 0.000025e-1] )
pn.append( [5, 2.59444e-1, 0.00007e-1] )
pn.append( [6, 2.94723e-1, 0.00009e-1] )
pn.append( [7, 1.98768e-1, 0.00007e-1] )
pn.append( [8, 9.0131e-2, 0.0006e-2] )
pn.append( [9, 2.9652e-2, 0.0002e-2] )
pn.append( [10, 7.4487e-3, 0.0006e-3] )
pn.append( [11, 1.4817e-3, 0.0002e-3] )
pn.append( [12, 2.4e-4, 0.0002e-4] )
pn.append( [13, 3.2324e-5, 0.0003e-5] )
pn.append( [14, 3.6835e-6, 0.0004e-6] )
pn.append( [15, 3.6017e-7, 0.0004e-7] )
pn.append( [16, 3.0574e-8, 0.0004e-8] )
pn.append( [17, 2.2762e-9, 0.0002e-9] )
pn.append( [18, 1.4989e-10, 0.0002e-10] )
pn.append( [19, 8.7983e-12, 0.0013e-12] )
pn.append( [20, 4.6314e-13, 0.0004e-12] )

def loadData(seed):
    filename = "poisson/N100000_S{0}.pickle.gz".format(seed)
    f = gzip.open(filename, 'rb')
    pstring = f.read()
    data = pickle.loads( pstring )
    f.close()
    print " loaded ",len(data), " sum=",sum(data.values())
    #print data
    return data

def mergeData(data1,data2):
    out = dict()
    for d1 in data1:
        if d1 in data2:
            out[d1]=data1[d1]+data2[d1]
        else:
            out[d1]=data1[d1]
    
    for d2 in data2:
        if d2 in data1:
            out[d2]=data1[d2]+data2[d2]
        else:
            out[d2]=data2[d2]
    
    
    return out

def loadDatas(slist):
    datas=[]
    for s in slist:
        datas.append( loadData(s) )
    # merge all datas into one
    #print datas
    out=dict()
    for d in datas:
        out = mergeData(d,out)
    return out
    
if __name__ == "__main__":  
    
    #data=loadData(0)
    max_seed=10
    data=loadDatas(range(max_seed))
    print " merged ",len(data), " sum=",sum(data.values())
    
    print data

    
    
    data2 = sortedDict(data)
    
    n=[]
    f=[]
    error=[]
    fsum=0
    for row in data2:
        #print row[0], " = ", row[1]
        n.append(row[0])
        f.append( row[1] )
        fsum = fsum+row[1]
        error.append( math.sqrt(row[1]) )
        
    print "fsum=",fsum
    fnorm=[]
    enorm_pos=[]
    enorm_neg=[]
    y_neg_limit = 1e-15
    for m in range(len(f)):
        fn = float(f[m])/float(fsum)
        fnorm.append( fn )
        en = float(error[m])/float(fsum)
        enorm_pos.append( en )
        if ( en >= fn ):
            enorm_neg.append( fn-y_neg_limit )
        else:
            enorm_neg.append(en)
        
        #print fn, " +/- ", en

    psum = sum( fnorm )
    
    # calculate mean
    wsum=0
    sqsum=0
    weight=0
    for m in range(len(f)):
        wsum = wsum + f[m]*n[m]
        sqsum = sqsum + f[m]*n[m]
    
    n_mean = float(wsum)/float(fsum)
    print "psum= ",psum
    print "n_mean= ",n_mean
    
    # calculate asymptotic approximation pn0
    dx=0.1
    naxis = [x*dx for x in range(0,int(20/dx))]
    pn0=[]
    for nv in naxis:
        a=1/(4*math.pi*math.pi)
        b=math.pow((8*math.pi*math.pi),nv)
        c=math.gamma( 2*nv+1)
        pn0.append(a*b/c)
        
    print len(naxis)
    print len(pn0)
    #exit()
    P.figure()
    
    #P.text(5,0.0001,"$\sum{p_n}=1.0$ " )
    
    P.plot(n,fnorm,'ro',label="OpenVoronoi 11.10-171")
    #P.semilogy(n,fnorm,'ro',label="OpenVoronoi 11.10-171")
    P.errorbar(n,fnorm,yerr=[enorm_neg,enorm_pos],fmt='ro')
    #P.plot(naxis,pn0,'b',label="pn0")
    pn_x=[]
    pn_y=[]
    pn_pos_err=[]
    pn_neg_err=[]
    for row in pn:
        pn_x.append(row[0])
        pn_y.append(row[1])
        pn_pos_err.append(row[1]+row[2])
        pn_neg_err.append(row[1]-row[2])
    P.plot(pn_x,pn_y,'g-',label="Hilhorst 2007")
    P.plot(pn_x,pn_pos_err,'g--')
    P.plot(pn_x,pn_neg_err,'g--')

    P.xlim(2, 20)
    P.ylim(0, 0.4)
    P.title("Probability $p_n$ for a face with $n$ edges in a Poisson Voronoi diagram.")
    P.xlabel("n")
    P.ylabel("$p_n$")
    P.text(40,1e-7,"2011 December\nanders.e.e.wallin 'at' gmail.com\ngithub.com/aewallin/openvoronoi")
    P.legend()


    P.show()


    
    print "PYTHON All DONE."

