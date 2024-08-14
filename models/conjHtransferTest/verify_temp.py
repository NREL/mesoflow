import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

axialdir_char="x"

ds=yt.load(argv[1])
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d
probsize=prob_hi-prob_lo
ncells=ds.domain_dimensions
clength   = probsize[0]
cwidth    = probsize[1]
cdepth    = probsize[2]

print("probsize,clength:",probsize,clength)

res=ncells[0]
dx=clength/res
eps=1e-10
lb = yt.LineBuffer(ds, (prob_lo[0]+0.5*dx, (0.5+eps)*(prob_lo[1]+prob_hi[1]), (0.5+eps)*(prob_lo[2]+prob_hi[2])), \
    (prob_hi[0]-0.5*dx, (0.5+eps)*(prob_lo[1]+prob_hi[1]), (0.5+eps)*(prob_lo[2]+prob_hi[2])), res)
x = np.linspace(prob_lo[0]+0.5*dx,prob_hi[0]-0.5*dx,res)
fld_temp = lb["temperature"]

Tleft=1100.0
Tright=300.0
kgas=500.0
ksolid=5000.0
lgas=0.7*clength
lsolid=0.3*clength

Tinterface=((kgas/lgas)*Tleft+(ksolid/lsolid)*Tright)/(kgas/lgas+ksolid/lsolid)
Tempexact=np.zeros(res)
for i in range(res):
    if(x[i]<=lgas):
        Tempexact[i]=(x[i]/lgas)*(Tinterface-Tleft)+Tleft
    else:
        xint=x[i]-lgas
        Tempexact[i]=(xint/lsolid)*(Tright-Tinterface)+Tinterface

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,1,figsize=(8,4))
ax.plot(x,fld_temp,'k*',label="mflo",markersize=3,markevery=1)
ax.plot(x,Tempexact,'r-',label="exact",markersize=3,markevery=1)
ax.legend(loc="best")

dir_char="x"
fig.suptitle("temperature solution along "+dir_char+" direction")
plt.savefig("conjht_"+dir_char+".png")
#=======================================

