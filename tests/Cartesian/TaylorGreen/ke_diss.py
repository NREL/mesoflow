import numpy as np
import yt
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sys import argv
import matplotlib as mpl
import matplotlib.colors as colors
import glob
from matplotlib.colors import LogNorm
from scipy import stats
from mpi4py import MPI

fn_pattern = argv[1]
fn_list=[]
try:
    fn_list = sorted(glob.glob(fn_pattern), key=lambda f: int(f.split("plt")[1]))
except:
    if(fn_list==[]):
        print("using file of plotfiles..")
        infile=open(argv[1],'r')
        for line in infile:
            fn_list.append(line.split()[0])
        infile.close()

# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
lfn_list = fn_list[rank::nprocs]

ds = yt.load(lfn_list[0])
fields_load=["velx","vely","velz"]
mean_ke=np.zeros(len(lfn_list))
times=np.zeros(len(lfn_list))

print("total number of files:",len(lfn_list))

for i, fn in enumerate(lfn_list):
    print("doing file:",i)
    ds = yt.load(fn)
    filenum=int(fn.split("plt")[1])
    ad = ds.all_data();

    u    = np.array(ad["velx"])
    v    = np.array(ad["vely"])
    w    = np.array(ad["velz"])

    ke = 0.5*(u**2+v**2+w**2)

    mean_ke[i]=np.mean(ke);
    times[i] = ds.current_time;

temp=comm.gather(times)
if(rank==0):
    globaltime=np.array([y for x in temp for y in x])

temp=comm.gather(mean_ke)
if(rank==0):
    globalmeanke=np.array([y for x in temp for y in x])

if(rank==0):
    ziparr=np.array(sorted(zip(globaltime,globalmeanke)))
    np.savetxt("tg_ke.dat",ziparr)
    #np.savetxt("tg_diss.dat",np.transpose(np.array([time_new,diss])))
    print(ziparr.shape)
    time_new=0.5*(ziparr[0:-1,0]+ziparr[1:,0])
    diss=np.diff(ziparr[:,1])/(ziparr[1:,0]-ziparr[0:-1,0])
    arr1=np.transpose(np.array((time_new,diss)))
    #np.savetxt("tg_diss.dat",np.transpose(np.array([time_new,diss])))
    np.savetxt("tg_diss.dat",arr1)
