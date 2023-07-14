import yt
ds=yt.load("plt00050")
print(ds.field_list)

slc = yt.SlicePlot(ds, "z", "CH4")
slc.save()

#plot = yt.LinePlot(ds, [("boxlib","CH4"),("boxlib","C2H6")], (0,0.0625,0.0625),(1, 0.0625,0.0625), 100)
#plot.annotate_legend([("boxlib","CH4"),("boxlib","C2H6")])
#plot.set_x_unit("m")
#plot.save()