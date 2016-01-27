;@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

;loadct, 2, /silent
LOADCT, 11
tvlct, r, g, b, /get

f=32

print, 'getting data..'
;blah='../../laredata/seplength4/jcrit25_beta09/'
blah='../../laredata/seplength4/jcrit20/'
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

npts=256
nom="Data/"

mylscl=1e3
;mylscl=1e3

print, 'loading: ', nom
ds=getrdata(1,/ke)
inike=ds.ke(0)
estore=0
kestore=dblarr(npts)
startpt=dblarr(3,npts)

for i=1,1 DO BEGIN
 tflag=0
 ;print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(i,/ke)
 startpt[*,i-1]=[ds.x(0)/mylscl,ds.y(0)/mylscl,ds.z(0)/mylscl]
   
 maxke=max(ds.ke)
 kestore(i-1)=maxke
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
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)
epts=n_elements(estore)

ref=alog10(1.0e-12)
alke=alog10(kestore-2.0d0)
astore=fix((alke-ref)/max(alke-ref)*254)
t1=ref
t2=alog10(max(kestore)-2.0d0)

;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
;title='Particle peak energy gain (eV)', /vertical, /ylog

;particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], zscl=1e6, /mm, myxyrange=[-1.0,1.0], myzrange=[-1,5.0], /rel, /symb, /bsymb
;particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e3, /km, myxyrange=[-10.0,10.0], myzrange=[-10,50.0], /rel, /symb, /bsymb

;particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], zscl=1e0, /me, myxyrange=[-1.0,1.0], myzrange=[-1.0,5.0], /rel, /symb, lscl=1e0

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
;bx=bx(*,1:ny,1:nz)
;by=by(1:nx,*,1:nz)
;bz=bz(1:nx,1:ny,*)

;Move B to the middle of the box
;bx = 0.5*(bx(0:nx-1,*,*)+bx(1:nx,*,*))
;by = 0.5*(by(*,0:ny-1,*)+by(*,1:ny,*))
;bz = 0.5*(bz(*,*,0:nz-1)+bz(*,*,1:nz))

bgrid = dblarr(nx,ny,nz,3)
bgrid[*,*,*,0] = bx(0:nx-1,0:ny-1,0:nz-1)
bgrid[*,*,*,1] = by(0:nx-1,0:ny-1,0:nz-1)
bgrid[*,*,*,2] = bz(0:nx-1,0:ny-1,0:nz-1)

print, 'destaggering'
bgrid=destaggerB(bgrid)

;Move grid points so they are at the same locations as B in order to run 
;the null finding code
;xx = 0.5*(x(1:nx)+x(0:nx-1))
;yy = 0.5*(y(1:ny)+y(0:ny-1))
;zz = 0.5*(z(1:nz)+z(0:nz-1))
xx = x(1:nx)
yy = y(1:ny)
zz = z(1:nz)


ch = 0.01
chmin = 0.001
chmax = 0.1
cepsilon = 1.0d-5
cmxline = 10000
csz_stp = size(startpt)

;trace field line forwards
linef = tracefield_3D(startpt[*,0],bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
ii = where(linef[0,*] gt -8000,nii)
;if (ii[0] ne -1) then xplot3d, linef(0,ii)*mylscl/1e6, linef(1,ii)*mylscl/1e6, linef(2,ii)*mylscl/1e6, COLOR=[0,255,0], /OVERPLOT
;trace field line backwards
lineb = tracefield_3D(startpt[*,0],bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
ii = where(lineb[0,*] gt -8000,nii)
;if (ii[0] ne -1) then xplot3d, lineb(0,ii)*mylscl/1e6, lineb(1,ii)*mylscl/1e6, lineb(2,ii)*mylscl/1e6, COLOR=[0,70,0], /OVERPLOT


;for i=2,npts DO BEGIN
; particletrack, i, /op, floc=nom, xyzt, zscl=1e3, /km, myxyrange=[-10.0,10.0], myzrange=[-10,50.0], /rel, /symb, /bsymb, lcol=[r(astore(0)),g(astore(0)),b(astore(0))]
; linef = tracefield_3D(startpt[*,i-1],bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
; ii = where(linef[0,*] gt -8000,nii)
; if (ii[0] ne -1) then xplot3d, linef(0,ii)*mylscl/1e6, linef(1,ii)*mylscl/1e6, linef(2,ii)*mylscl/1e6, COLOR=[0,255,0], /OVERPLOT
; ;if (ii[0] ne -1) then xplot3d, linef(0,ii)*mylscl/1e0, linef(1,ii)*mylscl/1e0, linef(2,ii)*mylscl/1e0, COLOR=[0,255,0], /OVERPLOT
; ;trace field line backwards
; lineb = tracefield_3D(startpt[*,i-1],bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
; ii = where(lineb[0,*] gt -8000,nii)
; if (ii[0] ne -1) then xplot3d, lineb(0,ii)*mylscl/1e6, lineb(1,ii)*mylscl/1e6, lineb(2,ii)*mylscl/1e6, COLOR=[0,70,0], /OVERPLOT
; ;if (ii[0] ne -1) then xplot3d, lineb(0,ii)*mylscl/1e0, lineb(1,ii)*mylscl/1e0, lineb(2,ii)*mylscl/1e0, COLOR=[170,0,220], /OVERPLOT
;endfor
;stop
; can we overplot an isosurface of current?

J=dblarr(nx, ny,nz, 3)

;dbxdx=0.5*(bgrid[2:nx-1,*,*,0]-bgrid[0:nx-3,*,*,0])/delx
dbydx=0.5*(bgrid[2:nx-1,*,*,1]-bgrid[0:nx-3,*,*,1])/delx
dbzdx=0.5*(bgrid[2:nx-1,*,*,2]-bgrid[0:nx-3,*,*,2])/delx
dbxdy=0.5*(bgrid[*,2:ny-1,*,0]-bgrid[*,0:ny-3,*,0])/dely
;dbydy=0.5*(bgrid[*,2:ny-1,*,1]-bgrid[*,0:ny-3,*,1])/dely
dbzdy=0.5*(bgrid[*,2:ny-1,*,2]-bgrid[*,0:ny-3,*,2])/dely
dbxdz=0.5*(bgrid[*,*,2:nz-1,0]-bgrid[*,*,0:nz-3,0])/delz
dbydz=0.5*(bgrid[*,*,2:nz-1,1]-bgrid[*,*,0:nz-3,1])/delz
;dbzdz=0.5*(bgrid[*,*,2:nz-1,2]-bgrid[*,*,0:nz-3,2])/delz

j[1:nx-2,1:ny-2,1:nz-2,0]=dbzdy[1:nx-2,*,1:nz-2]-dbydz[1:nx-2,1:ny-2,*]
j[1:nx-2,1:ny-2,1:nz-2,1]=dbxdz[1:nx-2,1:ny-2,*]-dbzdx[*,1:ny-2,1:nz-2]
j[1:nx-2,1:ny-2,1:nz-2,2]=dbydx[*,1:ny-2,1:nz-2]-dbxdy[1:nx-2,*,1:nz-2]

undefine, dbydz, dbzdx, dbxdy, dbzdy, dbxdz, dbydz


modj=sqrt(j[*,*,*,0]*j[*,*,*,0]+j[*,*,*,1]*j[*,*,*,1]+j[*,*,*,2]*j[*,*,*,2])

undefine, j

 ;modb=sqrt(bgrid[*,*,*,0]*bgrid[*,*,*,0]+bgrid[*,*,*,1]*bgrid[*,*,*,1]+bgrid[*,*,*,2]*bgrid[*,*,*,2])


 lsize=10
 theModel = Obj_New('IDLgrModel')
 theModel -> Scale, (max(xx)-min(xx))/DOUBLE(nx), (max(yy)-min(yy))/DOUBLE(ny), (max(zz)-min(zz))/DOUBLE(nz)	    	; julie's box dimensions
 theModel -> Add, mk_iso(modj,20,scol=[0, 255 ,0], /low)	; note threshold choice as 2nd arg of mk_iso
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=lsize
 if obj_valid(test) THEN test->SetProperty, COLOR = [0, 240, 0]
 xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.0*lsize,-1.0*lsize], symbol=test, /overplot




stop
ds=getrdata(1, /fields)
ppath=[[ds.x],[ds.y],[ds.z]]
npath=size(ppath)
modb=sqrt(ds.b[0,*]*ds.b[0,*]+ds.b[1,*]*ds.b[1,*]+ds.b[2,*]*ds.b[2,*])


;can we plot field line arrows at the same places?   
   
  ; dims=[1.,1.,6./2.]
  ; mysymsize=[mylscl,mylscl,mylscl]*0.05*dims
  ; flc=[0,0,200]
  ; as=[ [ds.x[0],ds.x[0]],[ds.y[0],ds.y[0]],[ds.z[0],ds.z[0]]]/1e3
   
  ; oModela=obj_new('IDLgrModel')
  ; oModela->add,mk_vector([ds.b[0,0],ds.b[1,0],ds.b[2,0]]/modb[0]/dims,color=flc)
  ; oAR = obj_new('IDLgrSymbol', oModela)
  ; IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize*0.005
  ; xplot3d, as[*,0], as[*,1], as[*,2], COLOR=flc, NAME='atest', SYMBOL=oAR, THICK=2, /OVERPLOT   	; plot arrow of |B|

stop
;plot field lines for every point in the particle path
for i=2,npath(1)-1 DO BEGIN
 startpt=[ppath[i-1,0],ppath[i-1,1],ppath[i-1,2]]/mylscl
 ;particletrack, i, /op, floc=nom, xyzt, lcol=[0,0,0], zscl=1e3, /km, /rel, /symb,myxyrange=[-10.0,10.0], myzrange=[-10.0,50.0]
 linef = tracefield_3D(startpt,bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
 ii = where(linef[0,*] gt -8000,nii)
 if (ii[0] ne -1) then xplot3d, linef(0,ii)*mylscl/1e3, linef(1,ii)*mylscl/1e3, linef(2,ii)*mylscl/1e3, COLOR=[0,255,0], /OVERPLOT
 lineb = tracefield_3D(startpt,bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
 ii = where(lineb[0,*] gt -8000,nii)
 if (ii[0] ne -1) then xplot3d, lineb(0,ii)*mylscl/1e3, lineb(1,ii)*mylscl/1e3, lineb(2,ii)*mylscl/1e3, COLOR=[170,0,220], /OVERPLOT
 
 ;as=[ [ds.x[i],ds.x[i]],[ds.y[i],ds.y[i]],[ds.z[i],ds.z[i]]]/1e3
 ;oModela=obj_new('IDLgrModel')
 ;oModela->add,mk_vector([ds.b[0,i],ds.b[1,i],ds.b[2,i]]/modb[i]/dims,color=flc)
 ;oAR = obj_new('IDLgrSymbol', oModela)
 ;IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize*0.005
 ;xplot3d, as[*,0], as[*,1], as[*,2], COLOR=flc, NAME='atest', SYMBOL=oAR, THICK=2, /OVERPLOT   	; plot arrow of |B|
endfor



END
