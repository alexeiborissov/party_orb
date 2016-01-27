;@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 2, /silent

f=41

print, 'getting data..'
blah='../../laredata/seplength4/jcrit25_beta09/'
jnom=blah+string(f, format='("newJ", I2,".sav")')
blah2=blah+'Data/'

;ds=getdata(41,/vx,/vy,/vz,/bx,/by,/bz)
;ds=getdata(f,/grid, wkdir="../../julie_dam_files/Data")
;ds=getdata(f,/grid, wkdir='../../laredata/seplength4/jcrit20/Data/')
ds=getdata(f,/grid, wkdir=blah2)

grid=ds.grid
delx=grid.x[1]-grid.x[0]
dely=grid.y[1]-grid.y[0]
delz=grid.z[1]-grid.z[0]

x=grid.x
nx = n_elements(x)-1
dx = x(1)-x(0)
y=grid.y
ny = n_elements(y)-1
dy = y(1)-y(0)
z=grid.z
nz = n_elements(z)-1
dz = z(1)-z(0)

data = getdata(f,wkdir=blah2,/bx) & bx=data.bx
data = getdata(f,wkdir=blah2,/by) & by=data.by
data = getdata(f,wkdir=blah2,/bz) & bz=data.bz

;Cut B since there is an extra x/y/z position for bx/by/bz
bx=bx(*,1:ny,1:nz)
by=by(1:nx,*,1:nz)
bz=bz(1:nx,1:ny,*)

;Move B to the middle of the box
bx = 0.5*(bx(0:nx-1,*,*)+bx(1:nx,*,*))
by = 0.5*(by(*,0:ny-1,*)+by(*,1:ny,*))
bz = 0.5*(bz(*,*,0:nz-1)+bz(*,*,1:nz))

;Move grid points so they are at the same locations as B in order to run 
;the null finding code
xx = 0.5*(x(1:nx)+x(0:nx-1))
yy = 0.5*(y(1:ny)+y(0:ny-1))
zz = 0.5*(z(1:nz)+z(0:nz-1))

bgrid = dblarr(nx,ny,nz,3)
bgrid[*,*,*,0] = bx
bgrid[*,*,*,1] = by
bgrid[*,*,*,2] = bz

;restore, jnom, /verbose
;modj=sqrt(j[*,*,*,0]*j[*,*,*,0]+j[*,*,*,1]*j[*,*,*,1]+j[*,*,*,2]*j[*,*,*,2])
;window, 0
;contour, modj[*,*,384], grid.x, grid.y, /fill, nlev=20, /iso, xr=[-0.6,0.6], yr=[-0.6,0.6], title='|j(z=0)|', ytitle='y', xtitle='x'
;loadct, 4, /silent
;tvlct, r, g, b, /get

npts=1
nom="Data/"
;nom="../relcodev1.0/Data/"
;nom="../teq0/Data/"

mylscl=1e4

print, 'loading: ', nom
ds=getrdata(1,/ke)
inike=ds.ke(0)
loadct, 2
estore=0
kestore=dblarr(npts)
startpt=dblarr(3,npts)
for i=1,npts DO BEGIN
 tflag=0
 ;print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(i,/ke)
 startpt[*,i-1]=[ds.x(0)/mylscl,ds.y(0)/mylscl,ds.z(0)/mylscl]
   
 maxke=max(ds.ke)
 ;maxketime=ds.time(maxkeloc)
 kestore(i-1)=maxke
; if ds.z[0] eq 0.5e6 THEN oplot, [ds.x[0], ds.x[0]]*1e-6, [ds.y[0],ds.y[0]]*1e-6, thick=2, psym=4, symsize=3
 s=n_elements(ds.t)
 IF (abs(ds.t[s-1]-1.0d0) lt 0.01) THEN BEGIN
  print, i, npts, ds.x[0],ds.y[0], ds.z[0], $
  format='(i4,"/",i4, " normal time exit [R= ", E10.3,",", E10.3, ",", E10.3,"]")' 
  tflag=1
 ENDIF
 IF ((abs(ds.x[s-1]) ge 0.58e4) OR (abs(ds.y[s-1]) ge 0.58e4)) THEN BEGIN
  print, i, npts, ds.x[0],ds.y[0], ds.z[0], $
  format='(i4,"/",i4, " normal spatial exit [R= ", E10.3,",", E10.3, ",", E10.3,"]")' 
  tflag=1
 ENDIF 
 IF (tflag eq 0) THEN BEGIN
  print, i, npts, ds.t[s-1], format='(i4,"/",i4, " unknown exit at t=", f11.4)'
  print, ds.x[0],ds.y[0], ds.z[0], ds.x[s-1],ds.y[s-1], ds.z[s-1], $
  format='(" [Ini pos = ", E10.3,",", E10.3,",", E10.3," || Fin pos = ", E10.3,",", E10.3,",", E10.3," ]")' 
 
  estore=[estore,i]
 ENDIF 
; IF ((tflag eq 0) and (ds.z[0] eq 0.45e4)) THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2, col=50
; IF ((tflag eq 0) and (ds.z[0] eq 0.5e4)) THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2, col=130
; IF ((tflag eq 0) and (ds.z[0] eq 0.55e4)) THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2, col=220
; IF ((tflag eq 0) and (ds.z[0] eq 2e4)) THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2, col=130
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)
epts=n_elements(estore)
;STOP

;save, kestore, cstore, /compress, /verbose, filename=nom+'../IDL/colorstore_test.sav'
;restore, filename=nom+'../IDL/colorstore_test.sav',/verbose

ref=alog10(1.0e-12)
alke=alog10(kestore-2.0d0)
astore=fix((alke-ref)/max(alke-ref)*254)
t1=ref
t2=alog10(max(kestore)-2.0d0)


;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
;title='Particle peak energy gain (eV)', /vertical, /ylog


particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], zscl=1e3, /km, myxyrange=[-10.0,10.0], myzrange=[-10,50.0], /rel, /symb

;stop
;stolen clare code snippet:

ch = 0.01
chmin = 0.001
chmax = 0.1
cepsilon = 1.0d-5
cmxline = 10000
csz_stp = size(startpt)

;trace field line forwards
linef = tracefield_3D(startpt[*,0],bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
ii = where(linef[0,*] gt -8000,nii)
if (ii[0] ne -1) then xplot3d, linef(0,ii)*mylscl/1e3, linef(1,ii)*mylscl/1e3, linef(2,ii)*mylscl/1e3, COLOR=[0,255,0], /OVERPLOT
;trace field line backwards
lineb = tracefield_3D(startpt[*,0],bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
ii = where(lineb[0,*] gt -8000,nii)
if (ii[0] ne -1) then xplot3d, lineb(0,ii)*mylscl/1e3, lineb(1,ii)*mylscl/1e3, lineb(2,ii)*mylscl/1e3, COLOR=[170,0,220], /OVERPLOT

STOP
for i=2,npts-1 DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[0,0,0], zscl=1e3, /km, /rel, /symb,myxyrange=[-10.0,10.0], myzrange=[-10.0,50.0]
 linef = tracefield_3D(startpt[*,i-1],bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
 ii = where(linef[0,*] gt -8000,nii)
 if (ii[0] ne -1) then xplot3d, linef(0,ii)*mylscl/1e3, linef(1,ii)*mylscl/1e3, linef(2,ii)*mylscl/1e3, COLOR=[0,255,0], /OVERPLOT
 lineb = tracefield_3D(startpt[*,i-1],bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
 ii = where(lineb[0,*] gt -8000,nii)
 if (ii[0] ne -1) then xplot3d, lineb(0,ii)*mylscl/1e3, lineb(1,ii)*mylscl/1e3, lineb(2,ii)*mylscl/1e3, COLOR=[170,0,220], /OVERPLOT
endfor

;particletrack, estore[1], floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e3, $
;/km, myxyrange=[-10.0,10.0], myzrange=[-17.5,27.5], /rel, /symb
;for i=2,epts-1 DO BEGIN
; particletrack, estore[i], /op, floc=nom, xyzt, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e3, /km, ;/rel, /symb,myxyrange=[-10.0,10.0], myzrange=[-10.0,50.0]
;endfor

END
