;common somename, tt, FRon
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

nulls=nullfinder(bv)
newton_improve, bv, nulls
;print, 'andrews null finder:'
print, nulls
;print, '--'

;ninfo3d, reform(bv[*,*,*,0]), reform(bv[*,*,*,1]), reform(bv[*,*,*,2]), laregrid.x, laregrid.y, laregrid.z, laregrid.npts, ninfo, ncan, px, py, pz, candidates
;;ninfo3d, ds.bx, ds.by, ds.bz, laregrid.x, laregrid.y, laregrid.z, laregrid.npts, ninfo, ncan, px, py, pz, candidates
;print, "ncandidates:", ncan
;print, "coords:", ninfo.coords

npts=8*8*5
nom="Data/"
nidl="IDL/"

save, nulls, filename=nidl+snapnom, /compress
restore, filename=nidl+snapnom, /verbose

STOP




;nom="../relcodev1.0/Data/"
;nom="../teq0/Data/"

print, 'loading: ', nom
ds=getrdata(1,/ke)
inike=ds.ke(0)

kestore=dblarr(npts)
for i=1,npts DO BEGIN
 print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(i,/ke)  
 maxke=max(ds.ke)
 ;maxketime=ds.time(maxkeloc)
 kestore(i-1)=maxke
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

save, kestore, cstore, /compress, /verbose, filename=nom+'../IDL/colorstore_test.sav'
restore, filename=nom+'../IDL/colorstore_test.sav',/verbose

ref=alog10(1.0e-12)

alke=alog10(kestore-2.0d0)
;astore=fix((alke-min(alke))/max(alke-min(alke))*254)
astore=fix((alke-ref)/max(alke-ref)*254)

t1=ref
;t1=alog10(min(kestore)-2.0d0)
t2=alog10(max(kestore)-2.0d0)

;make a new window for the colorbar
;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], title='Particle peak energy gain (eV)', /vertical, /ylog

;make a new window for the colorbar
cgwindow
cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], title='Particle peak energy gain (eV)', /vertical, /ylog

particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e6, /Mm, myxyrange=[-10,10], myzrange=[-17.5,27.5], /rel
for i=2,npts DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e6, /Mm, /rel
endfor


STOP





iniorbplot, 1, op=0, floc=nom, ocol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e6, /me, /rel, myxyrange=[-40.,40.], myzrange=[-20,20]

for i=2,npts DO BEGIN
print, i, npts, format='(i4,"/",i4)'
 iniorbplot, i, /op, floc=nom, ocol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e6, /me, /rel, myxyrange=[-40.,40.], myzrange=[-20,20]
endfor

        
END
