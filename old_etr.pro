;common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

filename='test.png'
ax=-55
az=30
filewidth=1000
transparent=0.5

loadct, 4, /silent
tvlct, r, g, b, /get
loadct, 4
;stretch, 30,270
tvlct, RR, GG, BB, /GET

npts=16*16*5

result = FILE_TEST('kestore.sav') 
IF (result eq 0) THEN BEGIN
;stop
; npts=16*16*2
 nom="Data/"
 print, 'loading: ', nom
 ds=getrdata(1,/ke, /gyro)
 inike=ds.ke(0)
 ;loadct, 2
 estore=0
 kestore=dblarr(npts)
 maxgyro=dblarr(npts)
 tstop=dblarr(npts)
 startpt=dblarr(3,npts)
 endpt=dblarr(3,npts)
 ns=n_elements(ds.t)
 kestore[0]=inike
 maxg=max(ds.gyror)
 maxgyro[0]=maxg
 startpt[*,0]=[ds.x[0],ds.y[0],ds.z[0]]
 endpt[*,0]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
 tstop[0]=ds.t[ns-1]
 for i=1,npts DO BEGIN
  tflag=0
  print, i, npts, format='(i4,"/",i4)'
  ds=getrdata(i,/ke,/gyro)
  ns=n_elements(ds.t)
  startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
  endpt[*,i-1]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
  tstop[i-1]=ds.t[ns-1]
  maxke=max(ds.ke)
  kestore[i-1]=maxke
  print, startpt[*,i-1], 'max ke=', maxke
  maxg=max(ds.gyror)
  maxgyro[i-1]=maxg
 endfor
 save, npts, kestore, startpt, endpt, tstop, maxgyro, filename='kestore.sav'
ENDIF ELSE BEGIN
;STOP
restore, filename='kestore.sav', /verbose
ENDELSE
;STOP

ref=alog10(1e-7)
topref=alog10(1e5)

alke=alog10(kestore-20.0d0)
astore=fix((alke-ref)/max(topref-ref)*254)
zerof=where(astore le 0)
astore(zerof)=0
maxf=where(astore ge 256)
astore(maxf)=255
t1=ref
t2=topref

arbscl=1e7
bscl=0.01d0
tscl=100.0d0
q=1.60217653E-19

fac=arbscl*arbscl*bscl/tscl
ofac=1.0d0/fac
ofac10=10.0d0/fac
ofac100=100.0d0/fac
ofac140k=140000.0d0/fac

mytit=' ';'energy gain *50keV [B/10G][L/100km]^2[T/20s]^(-1)'

;STOP
xt='x (Mm [L/10Mm])'
yt='y (Mm [L/10Mm])'
zt='z (Mm [L/10Mm])'

;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1*ofac10,10^t2*ofac10], $ 
;title=mytit, /vertical, /ylog
cgwindow
cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1*ofac140k,10^t2*ofac140k], $ 
title=mytit, /xlog

;window, 0
;plot, astore


axx=-77
azz=30
iniorbplot, 1, op=0, floc=nom, ocol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e6, /km, /rel, myxyrange=[-0.4,0.4], myzrange=[-20,20], ax=axx, az=azz, frac=0.09

for i=2,npts DO BEGIN
print, i, npts, format='(i4,"/",i4)'
 iniorbplot, i, /op, floc=nom, ocol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e6, /km, /rel, frac=0.09
endfor


STOP
xlen=DOUBLE(100-(-100))
ylen=DOUBLE(100-(-100))
zlen=DOUBLE(60-(-60))
frac=1.6
lscl=1e6
arbscl=1e6
iaz=0

 frac2=frac+astore(0)/256.*0.7
 oOrb2 = OBJ_NEW('orb', COLOR=[rr(astore(0)),gg(astore(0)),bb(astore(0))]) 	;final orb - coloured with ke
 ;oOrb2->Scale, frac, frac, zlen/xlen*frac
 oOrb2->Scale, frac2, frac2, zlen/xlen*frac2
 oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 for stop point
 ts=[ [startpt[0,0],startpt[0,0]],[startpt[1,0],startpt[1,0]],[startpt[2,0],startpt[2,0]]]
 te=[ [endpt[0,0],endpt[0,0]],[endpt[1,0],endpt[1,0]],[endpt[2,0],endpt[2,0]]]
 

 xplot3d, [-100,-100]*arbscl/lscl, [-100,-100]*arbscl/lscl, [-60,-60]*arbscl/lscl, $
 xr=[-100,100]*arbscl/lscl,yr=[-100,100]*arbscl/lscl,zr=[-60,60]*arbscl/lscl, $
 ax=ax, az=30+iaz, filewidth=filewidth, xtitle=xt, ytitle=yt, ztitle=zt
 
 fwrap, rad=1e7, frac=50.0, zloc=2e7, flinethick=3, /op

;STOP

 XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
 ;IF (astore(0) gt 128) THEN XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/lscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
 XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/lscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
 ;ii = where(linefs[*,0,0] gt -8000,nii)
 ;if (ii[0] ne -1) then xplot3d, linefs[0,ii,0]*arbscl/lscl, linefs[1,ii,0]*arbscl/lscl, linefs[2,ii,0]*arbscl/lscl, COLOR=grey, /OVERPLOT
 ;ii = where(linebs[*,0,0] gt -8000,nii)
 ;if (ii[0] ne -1) then xplot3d, linebs[0,ii,0]*arbscl/lscl, linebs[1,ii,0]*arbscl/lscl, linebs[2,ii,0]*arbscl/lscl, COLOR=grey, /OVERPLOT

 
 FOR i=2,npts DO BEGIN
  frac2=frac+astore(i-1)/256.
  oOrb2 = OBJ_NEW('orb', COLOR=[rr(astore(i-1)),gg(astore(i-1)),bb(astore(i-1))]) 	;final orb - coloured with ke
  oOrb2->Scale, frac2, frac2, zlen/xlen*frac2
  oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 for stop point
  ;print, i, npts, format='(i4,"/",i4)'
  ts=[ [startpt[0,i-1],startpt[0,i-1]],[startpt[1,i-1],startpt[1,i-1]],[startpt[2,i-1],startpt[2,i-1]]]
  te=[ [endpt[0,i-1],endpt[0,i-1]],[endpt[1,i-1],endpt[1,i-1]],[endpt[2,i-1],endpt[2,i-1]]]
  XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
  XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/lscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
 ENDFOR

;fieldwrapping, zloc=5.2e7, intercept=5e7, midrad=1e6, /op


STOP


STOP
 
particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e6, /mm, myxrange=[-100,100], myyrange=[-100,100], myzrange=[-60,60], /rel
for i=2,npts DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e6, /mm, /rel
endfor

STOP

t1=ref
;t1=alog10(min(kestore)-2.0d0)
t2=alog10(max(kestore)-2.0d0)

;make a new window for the colorbar
;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], title='Particle peak energy gain (eV)', /vertical, /ylog

;make a new window for the colorbar
cgwindow
cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], title='Particle peak energy gain (eV)', /vertical, /ylog
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,$
;range=[2,10^max(alke)], title='Particle peak energy (eV)', /vertical, /ylog
;cgColorbar, Divisions=6, Minor=5, Format='(e8.1)', Ticklen=-0.25, /window,$
;range=[min(kestore)-inike,max(kestore)-inike], title='Gain in eV over initial value (2eV)', /vertical

;STOP


        
END
