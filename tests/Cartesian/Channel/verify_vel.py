import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

channeldir = int(argv[2])
clength   = 1.0
cwidth    = 0.125
cdepth    = 0.125
ar     = (clength/cwidth)

channeldir_char=chr(ord('x')+channeldir)

ds=yt.load(argv[1])
fieldname="vel"+channeldir_char
res=100
slicedepth = cdepth/2
slc = ds.slice((channeldir+2)%3,slicedepth)
frb = slc.to_frb((1.0,'cm'),[res,int(res/ar)],height=(cwidth,'cm'))
y = np.linspace(0,cwidth,int(res/ar))
fld = np.array(frb[fieldname])[:,int(res/ar/2)]

exactsoln=y*(cwidth-y)*(4.0/cwidth**2)
#=======================================

#=======================================
#Plot solutions
#=======================================
print("plotting")
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,3))
ax1.plot(y,fld/np.max(fld),'k',label="mflo")
ax1.plot(y,exactsoln,'r*',label="exact")
ax1.legend(loc="best")

im=ax2.imshow(np.array(frb[fieldname]),origin="lower")
fig.colorbar(im, ax=ax2)
fig.suptitle("Channel flow along "+channeldir_char)

plt.savefig("vel_channel_"+channeldir_char+".png")
#=======================================

