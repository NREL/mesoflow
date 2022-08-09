import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob

def interpolate_ib(x,y,cell_i,cell_j,plox,ploy,dx,dy,field,vfrac):

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
sphererad=0.5

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
    dz_frb=probsize[2]/(res[2])
    midz=int(res[2]/2)
    midy=int(res[1]/2)
    fields_load=["volfrac","pressure"]
    plo=np.array([-5.0*sphererad,-2.0*sphererad,-2.0*sphererad])
    phi=plo+4.0*sphererad
    res1=np.array([4.0*sphererad/dx_frb,4.0*sphererad/dy_frb,4.0*sphererad/dz_frb]).astype(int)
    print(res1)
    print(plo)
    #ad = ds.covering_grid(level=maxlev, left_edge=prob_lo, dims=res, fields=fields_load)
    ad = ds.covering_grid(level=covgrid_lev, left_edge=plo, dims=res1, fields=fields_load)

    volfrac   = np.array(ad["volfrac"])
    pressure  = np.array(ad["pressure"])
    print(np.mean(volfrac))

    print(volfrac.shape)
    xlocs=np.array([])
    ylocs=np.array([])
    cploc=np.array([])
    presloc=np.array([])

    vfrac_line=volfrac[:,int(res1[1]/2),int(res1[2]/2)]
    pres_line=pressure[:,int(res1[1]/2),int(res1[2]/2)]
    xdist=np.linspace(plo[0],phi[0],res1[0])
    pres_grad=np.gradient(pres_line,xdist[1]-xdist[0])
    print("standoff distance:",-sphererad-xdist[np.argmax(pres_grad)])

    #angle=np.arctan(ylocs/xlocs)
    np.savetxt("presdata_%3.3d"%(fnum),np.transpose(np.vstack((xdist,pres_line,vfrac_line,pres_grad))),delimiter="  ")
    #np.savetxt("presdata_%3.3d"%(fnum),np.transpose(np.vstack((xlocs,ylocs,presloc))),delimiter="  ")
