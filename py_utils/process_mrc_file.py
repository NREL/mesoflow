import matplotlib.pyplot as plt
import mrcfile
import numpy as np
from sys import argv
from scipy.ndimage import zoom


#=====================================================================
def write_paraview_file_cartmesh(fname,dx,prob_lo,N,ncdata,ccdata):
    
    zero=0
    one=1
    outfile=open(fname,'w')
    outfile.write("<?xml version=\"1.0\"?>\n")
    outfile.write("<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
    outfile.write("<RectilinearGrid WholeExtent=\"%d\t%d\t%d\t%d\t%d\t%d\">\n"%(one,N[0],one,N[1],one,N[2]))
    outfile.write("<Piece Extent=\"%d\t%d\t%d\t%d\t%d\t%d\">\n"%(one,N[0],one,N[1],one,N[2]))

    outfile.write("<PointData>\n")
    n_ncdata=ncdata.shape[0]
    if(n_ncdata > 0):
        for ndataset in range(n_ncdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Point_data%d\" format=\"ascii\">\n"%(ndataset))

            for i in range(ncdata.shape[1]):
                outfile.write("%e "%(ncdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</PointData>\n")

    outfile.write("<CellData>\n")
    n_ccdata=ccdata.shape[0]
    if(n_ccdata > 0):
        for ndataset in range(n_ccdata):
            outfile.write("<DataArray type=\"Float32\" Name=\"Cell_data%d\" format=\"ascii\">\n"%(ndataset))
            for i in range(ccdata.shape[1]):
                outfile.write("%e "%(ccdata[ndataset][i]))
            outfile.write("\n</DataArray>\n")
    outfile.write("</CellData>\n")

    outfile.write("<Coordinates>\n")

    outfile.write("<DataArray type=\"Float32\" Name=\"X\"  format=\"ascii\">\n")
    for i in range(N[0]):
        outfile.write("%e\t"%(prob_lo[0]+i*dx[0]))
    outfile.write("\n</DataArray>\n")

    outfile.write("<DataArray type=\"Float32\" Name=\"Y\"  format=\"ascii\">\n")
    for i in range(N[1]):
        outfile.write("%e\t"%(prob_lo[1]+i*dx[1]))
    outfile.write("\n</DataArray>\n")
    
    outfile.write("<DataArray type=\"Float32\" Name=\"Z\"  format=\"ascii\">\n")
    for i in range(N[2]):
        outfile.write("%e\t"%(prob_lo[2]+i*dx[2]))
    outfile.write("\n</DataArray>\n")
    
    outfile.write("</Coordinates>\n")
    outfile.write("</Piece>\n")
    outfile.write("</RectilinearGrid>\n")
    outfile.write("</VTKFile>")

    outfile.close()
#=====================================================================

part_data=mrcfile.open(argv[1])


ad=part_data.data

try:
    restrict_ratio=float(argv[2])
except:
    restrict_ratio=1.0
ad_restricted = zoom(ad, (restrict_ratio, restrict_ratio, restrict_ratio))

nsize=np.shape(ad)
maxval=np.max(ad)

dx=np.array([0.1,0.1,0.1])
prob_lo=np.array([0.0,0.0,0.0])

ccdata=np.array([])
ncdata=np.array([ad_restricted.flatten()])
npts_xyz=np.shape(ad_restricted)
print(npts_xyz)
outfile=open(argv[2],"w")
outfile.write("%d\t%d\t%d\n"%(npts_xyz[2],npts_xyz[1],npts_xyz[0]))

for k in range(nsize[0]):
    for j in range(nsize[1]):
        for i in range(nsize[2]):
            outfile.write("%e\n"%(ad_restricted[k,j,i]))

write_paraview_file_cartmesh(argv[3],dx,prob_lo,np.array([npts_xyz[2],npts_xyz[1],npts_xyz[0]]),ncdata,ccdata)
