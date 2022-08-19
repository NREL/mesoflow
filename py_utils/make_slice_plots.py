import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob
    
    
fn_pattern= argv[1]
fieldname=argv[2]
fn_list = glob.glob(fn_pattern)
fn_list.sort()

set_minmax=False

print(fn_list)
if(len(argv) > 4):
    minval=float(argv[4])
    maxval=float(argv[5])
    set_minmax=True

for i, fn in enumerate(fn_list):
    ds=yt.load(fn)
    slc = yt.SlicePlot(ds, 'z', fieldname)
    slc.set_log(fieldname,False)
    if(set_minmax):
        slc.set_zlim(fieldname,minval,maxval)
    slc.annotate_grids()
    slc.save(fieldname+"_%4.4d"%(i))
