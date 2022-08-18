import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

axialdir = int(argv[2])
l1=0.7 #hard coded for now
clength   = 1.0
cwidth    = 0.125
cdepth    = 0.125
ar     = (clength/cwidth)

axialdir_char=chr(ord('x')+axialdir)

ds=yt.load(argv[1])
res=100
slicedepth = cdepth/2
slc = ds.slice((axialdir+2)%3,slicedepth)
frb = slc.to_frb((1.0,'cm'),res)
x = np.linspace(0,clength,res)
fld_S1 = np.array(frb["S1"])[res//2,:]
fld_S2 = np.array(frb["S2"])[res//2,:]

v=1.0;
D=0.5;
k=5;
w=np.sqrt(k/D)
l=clength

F=1/(np.exp(v*l1/D)-1)
A=v*(1+F)/(v*F*(np.exp(w*l1)-np.exp(w*(2*l-l1)))-D*w*(np.exp(w*l1)+np.exp(w*(2*l-l1))))
B=-A*np.exp(2*w*l)
c_l1=A*np.exp(w*l1)+B*np.exp(-w*l1)
c1=(c_l1-1)*F
exactsoln=(c1*(np.exp(v*x/D)-1)+1)*np.heaviside(l1-x,1.0)+(A*np.exp(w*x)+B*np.exp(-w*x))*np.heaviside(x-l1,1.0)
print(A,B,c1)
#=======================================

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(2,2,figsize=(8,4))
ax[0][0].plot(x,exactsoln,'r-',label="Exact solution")
ax[0][0].plot(x,fld_S1,'k*',label="mflo",markersize=2)
ax[0][0].legend(loc="best")

im=ax[0][1].imshow(np.array(frb["S1"]),origin="lower")
fig.colorbar(im, ax=ax[0][1])

#ax[1][0].plot(x,0.5*(x-x**2),'r-',label="Exact solution")
ax[1][0].plot(x,fld_S2,'k*',label="mflo",markersize=2)
ax[1][0].legend(loc="best")

im=ax[1][1].imshow(np.array(frb["S2"]),origin="lower")
fig.colorbar(im, ax=ax[1][1])

dir_char=chr(ord('x')+int(argv[2]))
fig.suptitle("S1 and S2 solution along "+dir_char+" direction")
plt.savefig("species_"+dir_char+".png")
#=======================================

