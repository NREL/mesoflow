import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

#=======================================
#read exact soln
#=======================================
infile=open("Lv_Ihme_JCP_2014",'r')

infile.readline()

x_exact = np.array([])
rho_exact = np.array([])

for line in infile:
    splt=line.split()
    x_exact=np.append(x_exact,float(splt[0]))
    rho_exact=np.append(rho_exact,float(splt[1]))

infile.close()
#=======================================

#=======================================
#read mflo solution
#=======================================
ds=yt.load(argv[1])

clength   = 1.0
cwidth    = 0.125
cdepth    = 0.125

fieldname="density"
res = 100
slicedepth = cdepth/2
slc = ds.slice((int(argv[2])+2)%3,slicedepth)
frb = slc.to_frb((1,'cm'), res)
x = np.linspace(0,1,res)
rho_left=0.96 #kg/m3
fld = np.array(frb[fieldname])[int(res/2),:]/rho_left
#=======================================

#=======================================
#Plot solutions
#=======================================
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(8,4))
ax1.plot(x_exact,rho_exact,'r-',label="Lv, Ihme, JCP, 2014")
ax1.plot(x,fld,'k*',label="mflo",markersize=2)
ax1.legend(loc="best")

im=ax2.imshow(np.array(frb[fieldname]),origin="lower")
fig.colorbar(im, ax=ax2)

dir_char=chr(ord('x')+int(argv[2]))
fig.suptitle("N2-He Shock tube density solution along "+dir_char+" direction")
plt.savefig("ms_stube_density_"+dir_char+".png")
#=======================================

