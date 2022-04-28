import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    
fieldname=argv[2]
ds=yt.load(argv[1])
slc = yt.SlicePlot(ds, 'z', fieldname)
slc.set_log(fieldname,False)
slc.annotate_grids()
slc.annotate_title(argv[3])
slc.save(argv[4])
