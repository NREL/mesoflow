import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob

def _magvort(field, data):
    return (np.sqrt( \
             (data["velz_gradient_y"]-data["vely_gradient_z"])**2 \
            +(data["velx_gradient_z"]-data["velz_gradient_x"])**2 \
            +(data["vely_gradient_x"]-data["velx_gradient_y"])**2))
    
fn_pattern= argv[1]
fn_list = glob.glob(fn_pattern)
fn_list.sort()
print(fn_list)
set_minmax=False
if(len(argv)>3):
    minval=float(argv[2])
    maxval=float(argv[3])
    set_minmax=True

for i, fn in enumerate(fn_list):
    ds=yt.load(fn)
    grad_fields1 = ds.add_gradient_fields(("boxlib","velx"))
    grad_fields2 = ds.add_gradient_fields(("boxlib","vely"))
    grad_fields3 = ds.add_gradient_fields(("boxlib","velz"))

    #note: yt thinks velx,vely,velz are dimensionless
    #so, units of vorticity comes out to be 1/cm
    ds.add_field(("boxlib", "magvort"), function=_magvort, units="1/cm")
    slc = yt.SlicePlot(ds, 'z', "magvort")
    slc.set_log("magvort",False)
    if(set_minmax):
        slc.set_zlim("magvort",minval,maxval)
    slc.annotate_grids()
    slc.save("vort"+"_%4.4d"%(i))
