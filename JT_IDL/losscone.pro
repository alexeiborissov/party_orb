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
;ax=0
;az=0
filewidth=800
transparent=0.5


;IF keyword_set(filename) and keyword_set(transparent) THEN BEGIN
; type=strsplit(filename,".",/extract)
; newfilenames=[type[0]+'_a.'+type[1],type[0]+'_b.'+type[1]]
;ENDIF

print, 'getting data..'
blah='../../../../2014/bourdin/data/VAR378_B.xdr'
blah2='../../../../2014/bourdin/data/VAR378_J_ABS.xdr'
blah3='../../../../2014/bourdin/data/VAR378_E_PARALLEL.xdr'
;jnom=blah+string(f, format='("newJ", I2,".sav")')
;blah2=blah+'Data/'

;ds=getdata(41,/vx,/vy,/vz,/bx,/by,/bz)
;ds=getdata(f,/grid, wkdir="../../julie_dam_files/Data")
;ds=getdata(f,/grid, wkdir='../../laredata/seplength4/jcrit20/Data/')
;ds=getdata(f,/grid, wkdir=blah2)
restore, /verbose, filename=blah3


;grid=ds.grid
delx=x[1]-x[0]
dely=y[1]-y[0]
;delz=z[1]-z[0]
delz=z-shift(z,-1)
;x=grid.x
nx = n_elements(x)-1
dx = x(1)-x(0)
;y=grid.y
ny = n_elements(y)-1
dy = y(1)-y(0)
;z=grid.z
nz = n_elements(z)-1
dz=z-shift(z,-1)

lscl=1e6

xscale=(max(x)-min(x))/DOUBLE(nx)/lscl
yscale=(max(y)-min(y))/DOUBLE(ny)/lscl
zscale=(max(z)-min(z))/DOUBLE(nz)/lscl
zscale2=(max(z)-z(40))/DOUBLE(nz)/lscl
xt='x(Mm)'
yt='y(Mm)'
zt='z(Mm)'


lee=60	;lowest z-gridpoint included

;c1=[0,0,255]	    	;blue
;c1a=[158, 182, 255]	;light blue
;c2=[255,0,0]	    	;red
;c2a=[255, 139, 100]	;light red

c1a=[255,153,204]
c1=[255,0,127]

;c1a=[255,153,255]	;pink
;c1=[153, 0, 153]	;dark pink
c2a=[153,51,255]	;purple
;c2=[127, 0, 255]	;dark purple
c2=[102, 0, 204]	;dark purple

;study currents above 
theModel = Obj_New('IDLgrModel')
theModel -> Scale, xscale,yscale,zscale2	    	    	
theModel -> Add, mk_iso(E_PARALLEL(*,*,lee:250),0.75,scol=c1, /low, minval=0.05)	
test=obj_new('IDLgrSymbol', theModel)
if obj_valid(test) THEN test->setproperty, size=1
if obj_valid(test) THEN test->SetProperty, COLOR = c1
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test, $
xr=[0,max(x)]/lscl, yr=[0,max(y)]/lscl, zr=[0,max(z)]/lscl, ax=ax, az=az, filewidth=filewidth,$
xtitle=xt, ytitle=yt, ztitle=zt

theModel2 = Obj_New('IDLgrModel')
theModel2 -> Scale, xscale,yscale,zscale2	    	    	
theModel2 -> Add, mk_iso(E_PARALLEL(*,*,lee:250),-0.75,scol=c2, maxval=-0.05)	
test2=obj_new('IDLgrSymbol', theModel2)
if obj_valid(test2) THEN test2->setproperty, size=1
if obj_valid(test2) THEN test2->SetProperty, COLOR = c2
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test2, /overplot

theModel3 = Obj_New('IDLgrModel')
theModel3 -> Scale, xscale,yscale,zscale2	    	    	
theModel3 -> Add, mk_iso(E_PARALLEL(*,*,lee:250),0.05,scol=c1a, /low, minval=0.0005d, transparency=0.2)	
test3=obj_new('IDLgrSymbol', theModel3)
if obj_valid(test3) THEN test3->setproperty, size=1
if obj_valid(test3) THEN test3->SetProperty, COLOR =c1a
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test3, /overplot

theModel4 = Obj_New('IDLgrModel')
theModel4 -> Scale, xscale,yscale,zscale2	    	    	
theModel4 -> Add, mk_iso(E_PARALLEL(*,*,lee:250),-0.05,scol=c2a, maxval=-0.0005d, transparency=0.2)	
test4=obj_new('IDLgrSymbol', theModel4)
if obj_valid(test4) THEN test4->setproperty, size=1
if obj_valid(test4) THEN test4->SetProperty, COLOR = c2a
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test4, /overplot

;stop
result = FILE_TEST('kestore.sav') 
IF (result eq 0) THEN BEGIN
;stop
 npts=20*20*3
 nom="Data/"
 print, 'loading: ', nom
 ds=getrdata(1,/ke)
 inike=ds.ke(0)
 ;loadct, 2
 estore=0
 kestore=dblarr(npts)
 startpt=dblarr(3,npts)
 kestore[0]=inike
 startpt[*,0]=[ds.x(0),ds.y(0),ds.z(0)]
 for i=1,npts DO BEGIN
  tflag=0
  print, i, npts, format='(i4,"/",i4)'
  ds=getrdata(i,/ke)
  startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
  maxke=max(ds.ke)
  kestore[i-1]=maxke
  print, startpt[*,i-1], 'max ke=', maxke
 endfor
 save, npts, kestore, startpt, filename='kestore.sav'
ENDIF ELSE BEGIN
;STOP
restore, filename='kestore.sav', /verbose
ENDELSE
;STOP
loadct, 4
;stretch, 30,270
tvlct, RR, GG, BB, /GET
;RR=REVERSE(RR)
;GG=REVERSE(GG)
;BB=REVERSE(BB)
;tvlct, RR, GG, BB
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

ref=alog10(20.0d0)
alke=alog10(kestore-20.0d0)
astore=fix((alke-ref)/max(alke-ref)*254)
t1=ref
t2=alog10(max(kestore)-20.0d0)

 ch = 0.01*dx 	    	    ;step size
 chmin = 0.001*dx	    	    ;min step size
 chmax = 0.5*dx	    	    ; max step size
 cepsilon = 1.0d-5  	    ; tolerance
 cmxline = 100000    	    ; max nsteps
 csz_stp = size(startpt)

 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]


cgwindow
cgColorbar, Divisions=3, Minor=5, Format='(e7.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
title='Particle peak energy gain (eV)', /xlog

restore, /verbose, filename=blah

   linef = tracefield_3D(startpt[*,24],b,x,y,z,ch,chmin,chmax,cepsilon,cmxline)
   iif = where(linef[0,*] gt -8000,nii)
   if (iif[0] ne -1) then xplot3d, linef(0,iif)/lscl, linef(1,iif)/lscl, linef(2,iif)/lscl, COLOR=steelblue, /OVERPLOT
   ;trace field line backward
   lineb = tracefield_3D(startpt[*,24],b,x,y,z,-ch,chmin,chmax,cepsilon,cmxline)
   iib = where(lineb[0,*] gt -8000,nii)
   if (iib[0] ne -1) then xplot3d, lineb(0,iib)/lscl, lineb(1,iib)/lscl, lineb(2,iib)/lscl, COLOR=steelblue, /OVERPLOT


 particletrack, 1, floc=nom, xyzt, lcol=[RR(astore(25)),GG(astore(25)),BB(astore(25))], /op, zscl=mylscl, /mm, lscl=lscl, /rel, /symb, /bsymb

iibmax=n_elements(iib)
iifmax=n_elements(iif)
lbstore=dblarr(4,iibmax)
lfstore=dblarr(4,iifmax)

FOR j=0, iibmax-1 DO BEGIN 
  temp=linterp(b, x, y, z, [lineb[0,j],lineb[1,j],lineb[2,j]])
  lbstore[0:2,j]=temp
  lbstore[3,j]=sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2])
ENDFOR
FOR j=0, iifmax-1 DO BEGIN 
  temp=linterp(b, x, y, z, [linef[0,j],linef[1,j],linef[2,j]])
  lfstore[0:2,j]=temp
  lfstore[3,j]=sqrt(temp[0]*temp[0]+temp[1]*temp[1]+temp[2]*temp[2])
ENDFOR

bmax=max( [[lbstore],[lfstore]])

bstore=reform([[lbstore[3,*]],[lfstore[3,*]]])


END
