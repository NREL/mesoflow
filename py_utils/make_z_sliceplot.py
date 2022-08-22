import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    
ds=yt.load(argv[1])
fieldname=argv[2]
slc = yt.SlicePlot(ds, 'z', fieldname)
slc.set_log(fieldname,False)
slc.annotate_grids()
slc.save(argv[2])

