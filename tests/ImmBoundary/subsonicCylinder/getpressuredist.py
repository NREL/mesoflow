import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy import interpolate
from sklearn.gaussian_process import GaussianProcessRegressor

def interpolate_ib_gp(x,y,cell_i,cell_j,plox,ploy,dx,dy,field,vfrac):

    predictptx=np.array([x])
    predictpty=np.array([y])
    xx=np.array([])
    yy=np.array([])
    fieldval=np.array([])
    gp = GaussianProcessRegressor()
    for i in range(-3,4):
        for j in range(-3,4):
            if(i!=0 and j!=0):
                if(vfrac[cell_i+i][cell_j+j]==1.0):
                    ptx=plox+(cell_i+i+0.5)*dx
                    pty=ploy+(cell_j+j+0.5)*dy
                    xx=np.append(xx,ptx)
                    yy=np.append(yy,pty)
                    fieldval=np.append(fieldval,field[cell_i+i][cell_j+j])

    gp.fit(np.transpose(np.vstack((xx,yy))), fieldval)
    return(gp.predict(np.array([[x,y]])))

def interpolate_ib_idw(x,y,cell_i,cell_j,plox,ploy,dx,dy,field,vfrac):

    sum_dist2=0.0;
    sum_field_dist2=0.0;
    for i in range(-1,2):
        for j in range(-1,2):
            if(i!=0 and j!=0):
                if(vfrac[cell_i+i][cell_j+j]==1.0):
                    
                    ptx=plox+(cell_i+i+0.5)*dx
                    pty=ploy+(cell_j+j+0.5)*dy
                    dist2=(x-ptx)**2+(y-pty)**2
                    sum_field_dist2+=field[cell_i+i][cell_j+j]/dist2
                    sum_dist2+=1.0/dist2

    return(sum_field_dist2/sum_dist2)


def _magvort(field, data):
    return (np.abs(data["vely_gradient_x"]-data["velx_gradient_y"]))
    
fn_pattern= argv[1]
fn_list = glob.glob(fn_pattern)
fn_list.sort()
print(fn_list)
pinf=1.0
rhoinf=1.0
vinf=0.236

for fnum, fn in enumerate(fn_list):
    
    ds=yt.load(fn)
    prob_lo=ds.domain_left_edge.d
    prob_hi=ds.domain_right_edge.d
    probsize=prob_hi-prob_lo
    ncells=ds.domain_dimensions
    print(ncells)
    minlev=0
    maxlev=ds.index.max_level
    lengths=prob_hi-prob_lo

    print(maxlev)
    covgrid_lev=maxlev
    res=np.array([ncells[0]* (2**covgrid_lev),ncells[1]* (2**covgrid_lev),ncells[2]* (2**covgrid_lev)])
    dx_frb=probsize[0]/(res[0])
    dy_frb=probsize[1]/(res[1])
    midz=int(res[2]/2)
    fields_load=["volfrac","pressure"]
    #ad = ds.covering_grid(level=maxlev, left_edge=prob_lo, dims=res, fields=fields_load)
    ad = ds.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load)

    volfrac   = np.array(ad["volfrac"])
    pressure  = np.array(ad["pressure"])
    print(np.mean(volfrac))

    print(volfrac.shape)
    xlocs=np.array([])
    ylocs=np.array([])
    cploc=np.array([])
    presloc=np.array([])

    vfrac_slice=volfrac[:,:,midz]
    pres_slice=pressure[:,:,midz]

    nangles=25
    angle=np.linspace(-np.pi,0.0,nangles)
    rad=0.5

    for i in range(nangles):
        
        x=rad*np.cos(angle[i])
        y=rad*np.sin(angle[i])

        #find enclosing cell
        cell_i=int(np.floor((x-prob_lo[0])/dx_frb))
        cell_j=int(np.floor((y-prob_lo[1])/dy_frb))

        pres_here=interpolate_ib_idw(x,y,cell_i,cell_j,prob_lo[0],prob_lo[1],dx_frb,dy_frb,\
                pres_slice,vfrac_slice)
        cp_here=(pres_here-pinf)/(0.5*rhoinf*vinf**2)
        presloc=np.append(presloc,pres_here)
        cploc=np.append(cploc,cp_here)

    #for i in range(res[0]):
    #    for j in range(res[1]):
    #        if(vfrac_slice[i,j]*(vfrac_slice[i,j]-1.0)<0.0):
    #            xlocs=np.append(xlocs,prob_lo[0]+(i+0.5)*dx_frb)
    #            ylocs=np.append(ylocs,prob_lo[1]+(j+0.5)*dy_frb)
    #            presloc=np.append(presloc,pres_slice[i][j])
    #            cploc=np.append(cploc,(pres_slice[i][j]-pinf)/(0.5*rhoinf*vinf**2))
                #print(vfrac_slice[i,j])

    #angle=np.arctan(ylocs/xlocs)
    ziparr=sorted(zip(angle+np.pi,presloc,cploc))
    np.savetxt("presdata_%3.3d"%(fnum),ziparr,delimiter="  ")
    #np.savetxt("presdata_%3.3d"%(fnum),np.transpose(np.vstack((xlocs,ylocs,presloc))),delimiter="  ")
