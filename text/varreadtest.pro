@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 4, /silent
tvlct, r, g, b, /get

snap=41
snapnom=string(snap,format='("l3dsnap",i2,".ninfo")')

print, "loading lare B-fields"
ds=getdata(snap,/bx, /by, /bz, wkdir="../julie_dat_files/Data")

laregrid=ds.grid

print, "correcting for staggered grid"
larefieldcorrect, ds.bx,ds.by,ds.bz, ds.grid.npts, newBx,newBy,newBz
undefine, ds

bv=dblarr(laregrid.npts[0],laregrid.npts[1],laregrid.npts[2],3)
bv[*,*,*,0]=newbx
bv[*,*,*,1]=newby
bv[*,*,*,2]=newbz
undefine, newbx, newby, newbz









END
