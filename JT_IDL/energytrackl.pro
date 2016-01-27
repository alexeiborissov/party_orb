;@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1
;@xplot3dJT
loadct, 4, /silent
tvlct, r, g, b, /get

f=32
filename='test.png'
ax=-60
az=30
filewidth=800
transparent=0.5


IF keyword_set(filename) and keyword_set(transparent) THEN BEGIN
 type=strsplit(filename,".",/extract)
 newfilenames=[type[0]+'_a.'+type[1],type[0]+'_b.'+type[1]]
ENDIF

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

bgrid = dblarr(nx,ny,nz,3)
bgrid[*,*,*,0] = bx(0:nx-1,0:ny-1,0:nz-1)
bgrid[*,*,*,1] = by(0:nx-1,0:ny-1,0:nz-1)
bgrid[*,*,*,2] = bz(0:nx-1,0:ny-1,0:nz-1)

print, 'destaggering'
bgrid=destaggerB(bgrid)

;Move grid points so they are at the same locations as B in order to run 
;the null finding code
xx = x(1:nx)
yy = y(1:ny)
zz = z(1:nz)


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
 lsize=1
 theModel = Obj_New('IDLgrModel')
 theModel -> Scale, (max(xx)-min(xx))/DOUBLE(nx), (max(yy)-min(yy))/DOUBLE(ny), (max(zz)-min(zz))/DOUBLE(nz)	    	; julie's box dimensions
 theModel -> Add, mk_iso(modj,20,scol=[255,255,255], /low)	; note threshold choice as 2nd arg of mk_iso
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=lsize
 if obj_valid(test) THEN test->SetProperty, COLOR = [0, 240, 0]
 
 

npts=256
;npts=320
nom="Data/"
;nom="../relcodev1.0/Data/"
;nom="../teq0/Data/"

arbscl=1e6
mylscl=1e6
;mylscl=1e3

print, 'loading: ', nom
ds=getrdata(1,/ke)
inike=ds.ke(0)
;loadct, 2
estore=0
kestore=dblarr(npts)
startpt=dblarr(3,npts)
for i=1,npts DO BEGIN
 tflag=0
 ;print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(i,/ke)
 startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]/arbscl
   
 maxke=max(ds.ke)
 ;maxketime=ds.time(maxkeloc)
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
;print, startpt[*,i-1]
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)
epts=n_elements(estore)

ref=alog10(1.0e-12)
alke=alog10(kestore-2.0d0)
astore=fix((alke-ref)/max(alke-ref)*254)
t1=ref
t2=alog10(max(kestore)-2.0d0)

 ch = 0.01
 chmin = 0.001
 chmax = 0.1
 cepsilon = 1.0d-5
 cmxline = 10000
 csz_stp = size(startpt)

 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]


;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
;title='Particle peak energy gain (eV)', /vertical, /ylog

print, 'making a movie..'
naz=0

FOR iii=0,naz DO BEGIN
 az=30+iii
 filename=type[0]+string(iii,format='(i3.3,".")')+type[1]
 print, iii, naz, format='("frame ", i3.3, "/", i3.3)'
IF (keyword_set(filename) and keyword_set(transparent)) THEN ns=1 ELSE ns=0
FOR ij=0,ns DO BEGIN  ; for transparency (as outputting closes the xWindow)

 ;particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=mylscl, /km, myxyrange=[-1,1]*arbscl/mylscl, myzrange=[-1,5]*arbscl/mylscl, /rel, /symb, /bsymb
 ;particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=mylscl, /km, myxyrange=[-1,1]*arbscl/mylscl, myzrange=[-1,5]*arbscl/mylscl, /rel, /symb, /bsymb
 particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=mylscl, /mm, myxyrange=[-1,1]*arbscl/mylscl, myzrange=[-1,5]*arbscl/mylscl, /rel, ax=ax, az=az, filewidth=filewidth


;trace field line forwards
 linef = tracefield_3D(startpt[*,0],bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
 ii = where(linef[0,*] gt -8000,nii)
 if (ii[0] ne -1) then xplot3d, linef(0,ii)*arbscl/mylscl, linef(1,ii)*arbscl/mylscl, linef(2,ii)*arbscl/mylscl, COLOR=steelblue, /OVERPLOT
 ;trace field line backwards
 lineb = tracefield_3D(startpt[*,0],bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
 ii = where(lineb[0,*] gt -8000,nii)
 if (ii[0] ne -1) then xplot3d, lineb(0,ii)*arbscl/mylscl, lineb(1,ii)*arbscl/mylscl, lineb(2,ii)*arbscl/mylscl, COLOR=steelblue, /OVERPLOT

 for i=2,npts DO BEGIN
  particletrack, i, /op, floc=nom, xyzt, zscl=mylscl, /mm, myxyrange=[-1,1]*arbscl/mylscl, myzrange=[-1,5]*arbscl/mylscl, /rel, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))]
  linef = tracefield_3D(startpt[*,i-1],bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
  ii = where(linef[0,*] gt -8000,nii)
  if (ii[0] ne -1) then xplot3d, linef(0,ii)*arbscl/mylscl, linef(1,ii)*arbscl/mylscl, linef(2,ii)*arbscl/mylscl, COLOR=steelblue, /OVERPLOT
  ;trace field line backwards
  lineb = tracefield_3D(startpt[*,i-1],bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
  ii = where(lineb[0,*] gt -8000,nii)
  if (ii[0] ne -1) then xplot3d, lineb(0,ii)*arbscl/mylscl, lineb(1,ii)*arbscl/mylscl, lineb(2,ii)*arbscl/mylscl, COLOR=steelblue, /OVERPLOT
 endfor

 ;xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.0*lsize,-1.0*lsize], symbol=test, /overplot
  IF ((ij eq 0) and keyword_set(filename) and keyword_set(transparent)) THEN xplot3d, lineb(0,ii)*arbscl/mylscl, lineb(1,ii)*arbscl/mylscl, lineb(2,ii)*arbscl/mylscl, COLOR=steelblue, /OVERPLOT, filename=newfilenames[0]
  IF ((ij eq 1) and keyword_set(filename) and keyword_set(transparent)) THEN xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.0*lsize,-1.0*lsize], symbol=test, /overplot, filename=newfilenames[1]
ENDFOR


IF (keyword_set(filename) AND keyword_set(transparent)) THEN BEGIN
 im1=read_image(newfilenames[0])
 im2=read_image(newfilenames[1])
 window, 0, xsize=filewidth, ysize=filewidth
 tv, im1, true=1
 tv, byte((1-transparent)*im1+transparent*im2), true=1
 write_image,filename,type[1],byte((1-transparent)*im1+transparent*im2)
 wdelete, 0
ENDIF
ENDFOR

END
