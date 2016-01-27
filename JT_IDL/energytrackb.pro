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

 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]
 
;IF keyword_set(filename) and keyword_set(transparent) THEN BEGIN
; type=strsplit(filename,".",/extract)
; newfilenames=[type[0]+'_a.'+type[1],type[0]+'_b.'+type[1]]
;ENDIF

print, 'getting data..'
;blah='../../laredata/seplength4/jcrit25_beta09/'
;blah='../../laredata/seplength4/jcrit20/'
blah='../../data/VAR378_B.xdr'
blah2='../../data/VAR378_J_ABS.xdr'
;jnom=blah+string(f, format='("newJ", I2,".sav")')
;blah2=blah+'Data/'

;ds=getdata(41,/vx,/vy,/vz,/bx,/by,/bz)
;ds=getdata(f,/grid, wkdir="../../julie_dam_files/Data")
;ds=getdata(f,/grid, wkdir='../../laredata/seplength4/jcrit20/Data/')
;ds=getdata(f,/grid, wkdir=blah2)
restore, /verbose, filename=blah2

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

;make isosurface of current
;theModel = Obj_New('IDLgrModel')
;theModel -> Scale, xscale,yscale,zscale	    	    	
;theModel -> Add, mk_iso(J_ABS,0.005,scol=[255,255,255], /low)	
;test=obj_new('IDLgrSymbol', theModel)
;if obj_valid(test) THEN test->setproperty, size=1
;if obj_valid(test) THEN test->SetProperty, COLOR = [0, 240, 0]
;xplot3d, [0,0], [0,0], [0,0], symbol=test, $
;xr=[0,max(x)]/lscl, yr=[0,max(y)]/lscl, zr=[0,max(z)]/lscl, ax=ax, az=az, filewidth=filewidth,$
;xtitle=xt, ytitle=yt, ztitle=zt

;J_ABS=J_ABS(*,*,40:*)

;study currents above 
theModel = Obj_New('IDLgrModel')
theModel -> Scale, xscale,yscale,zscale2	    	    	
theModel -> Add, mk_iso(J_ABS(*,*,40:*),0.00004,scol=steelblue, /low)	
test=obj_new('IDLgrSymbol', theModel)
if obj_valid(test) THEN test->setproperty, size=1
if obj_valid(test) THEN test->SetProperty, COLOR = [0, 240, 0]
xplot3d, [0,0], [0,0], [z(40),z(40)]/lscl, symbol=test, $
xr=[0,max(x)]/lscl, yr=[0,max(y)]/lscl, zr=[0,max(z)]/lscl, ax=ax, az=az, filewidth=filewidth,$
xtitle=xt, ytitle=yt, ztitle=zt

;STOP

npts=50
nom="Data/"
print, 'loading: ', nom
ds=getrdata(1,/ke)
inike=ds.ke(0)
estore=0
kestore=dblarr(npts)
startpt=dblarr(3,npts)
for i=1,npts DO BEGIN
 tflag=0
 print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(i,/ke)
 startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
 maxke=max(ds.ke)
 kestore[i-1]=maxke
 print, startpt[*,i-1]
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

ref=alog10(1.0e-12)
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

;origin=[[0.0,0.0],[0.0,0.0],[0.0,0.0]]
;xplot3d, origin[*,0], origin[*,1], origin[*,2], xr=[min(x),max(x)], yr=[min(y),max(y)], zr=[min(z),max(z)], /symbol, ax=ax, az=az, filewidth=filewidth

cgwindow
cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
title='Particle peak energy gain (eV)', /vertical, /ylog

;print, 'making a movie..'
;naz=0

;FOR iii=0,naz DO BEGIN
; az=30+iii
; filename=type[0]+string(iii,format='(i3.3,".")')+type[1]
; print, iii, naz, format='("frame ", i3.3, "/", i3.3)'
;IF (keyword_set(filename) and keyword_set(transparent)) THEN ns=1 ELSE ns=0
;FOR ij=0,ns DO BEGIN  ; for transparency (as outputting closes the xWindow)

 ;;;particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=mylscl, /mm, myxyrange=[-1,1]*arbscl/mylscl, myzrange=[-1,5]*arbscl/mylscl, /rel, ax=ax, az=az, filewidth=filewidth
 ;particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], /op, zscl=mylscl, /mm, lscl=lscl, /rel, /symb, /bsymb
 particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], /op, zscl=mylscl, /mm, lscl=lscl, /rel, /symb, /bsymb
;trace field line forwards
 ;newstartpt=dblarr(3,nfieldlines*nfieldlines)
 
restore, /verbose, filename=blah
;;TEST GENERAL FIELD TRACING ABILITY
; nxflines=10
; nyflines=10
; tempx=x(0)+max(x)*findgen(nxflines)/(nxflines+1)
; tempy=y(0)+max(y)*findgen(nyflines)/(nyflines+1)
 ;FOR iii=0,nxflines-1 DO BEGIN
 ; FOR jjj=0,nyflines-1 DO BEGIN
;   tempstartpt=[tempx(iii),tempy(jjj),1e7]
;   linef = tracefield_3D(tempstartpt,b,x,y,z,ch,chmin,chmax,cepsilon,cmxline)
;   ii = where(linef[0,*] gt -8000,nii)
 ;  if (ii[0] ne -1) then xplot3d, linef(0,ii)/lscl, linef(1,ii)/lscl, linef(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT
 ;  ;trace field line backward
 ;  lineb = tracefield_3D(tempstartpt,b,x,y,z,-ch,chmin,chmax,cepsilon,cmxline)
 ;  ii = where(lineb[0,*] gt -8000,nii)
 ;  if (ii[0] ne -1) then xplot3d, lineb(0,ii)/lscl, lineb(1,ii)/lscl, lineb(2,ii)/lscl, COLOR=brightpink, /OVERPLOT
 ; ENDFOR
 ;ENDFOR
 
   linef = tracefield_3D(startpt[*,0],b,x,y,z,ch,chmin,chmax,cepsilon,cmxline)
   ii = where(linef[0,*] gt -8000,nii)
   if (ii[0] ne -1) then xplot3d, linef(0,ii)/lscl, linef(1,ii)/lscl, linef(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT
   ;trace field line backward
   lineb = tracefield_3D(startpt[*,0],b,x,y,z,-ch,chmin,chmax,cepsilon,cmxline)
   ii = where(lineb[0,*] gt -8000,nii)
   if (ii[0] ne -1) then xplot3d, lineb(0,ii)/lscl, lineb(1,ii)/lscl, lineb(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT

 

 for i=2,npts DO BEGIN
  ;particletrack, i, floc=nom, xyzt, lcol=[0,0,0], /op, zscl=mylscl, /mm, lscl=lscl, /rel, /symb, /bsymb
  particletrack, i, floc=nom, xyzt, lcol=[r(astore(i)),g(astore(i)),b(astore(i))], /op, zscl=mylscl, /mm, lscl=lscl, /rel, /symb, /bsymb
  linef = tracefield_3D(startpt[*,i-1],b,x,y,z,ch,chmin,chmax,cepsilon,cmxline)
  ii = where(linef[0,*] gt -8000,nii)
  if (ii[0] ne -1) then xplot3d, linef(0,ii)/lscl, linef(1,ii)/lscl, linef(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT
  ;trace field line backwards
  lineb = tracefield_3D(startpt[*,i-1],b,x,y,z,-ch,chmin,chmax,cepsilon,cmxline)
  ii = where(lineb[0,*] gt -8000,nii)
  if (ii[0] ne -1) then xplot3d, lineb(0,ii)/lscl, lineb(1,ii)/lscl, lineb(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT
 endfor

 ;xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.0*lsize,-1.0*lsize], symbol=test, /overplot
 ; IF ((ij eq 0) and keyword_set(filename) and keyword_set(transparent)) THEN xplot3d, lineb(0,ii)*arbscl/mylscl, lineb(1,ii)*arbscl/mylscl, lineb(2,ii)*arbscl/mylscl, COLOR=steelblue, /OVERPLOT, filename=newfilenames[0]
 ; IF ((ij eq 1) and keyword_set(filename) and keyword_set(transparent)) THEN xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.0*lsize,-1.0*lsize], symbol=test, /overplot, filename=newfilenames[1]
;ENDFOR


;IF (keyword_set(filename) AND keyword_set(transparent)) THEN BEGIN
; im1=read_image(newfilenames[0])
; im2=read_image(newfilenames[1])
; window, 0, xsize=filewidth, ysize=filewidth
; tv, im1, true=1
; tv, byte((1-transparent)*im1+transparent*im2), true=1
; write_image,filename,type[1],byte((1-transparent)*im1+transparent*im2)
; wdelete, 0
;ENDIF
;ENDFOR

END
