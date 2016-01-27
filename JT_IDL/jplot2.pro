@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 2, /silent

f=31

print, 'getting data..'
;ds=getdata(41,/vx,/vy,/vz,/bx,/by,/bz)
;ds=getdata(f,/grid, wkdir="../../julie_dam_files/Data")
ds=getdata(f,/grid, wkdir='../../laredata/seplength4/jcrit20/Data/')

grid=ds.grid
delx=grid.x[1]-grid.x[0]
dely=grid.y[1]-grid.y[0]
delz=grid.z[1]-grid.z[0]


jnom="../../laredata/seplength4/jcrit20/"+string(f, format='("newJ", I2,".sav")')

restore, jnom, /verbose
modj=sqrt(j[*,*,*,0]*j[*,*,*,0]+j[*,*,*,1]*j[*,*,*,1]+j[*,*,*,2]*j[*,*,*,2])

window, 0


contour, modj[*,*,384], grid.x, grid.y, /fill, nlev=20, /iso, xr=[-0.6,0.6], yr=[-0.6,0.6], title='|j(z=0)|', ytitle='y', xtitle='x'

loadct, 4, /silent
tvlct, r, g, b, /get

npts=6*6*1
nom="Data/"
;nom="../relcodev1.0/Data/"
;nom="../teq0/Data/"

print, 'loading: ', nom
ds=getrdata(1,/ke)
inike=ds.ke(0)
loadct, 2
estore=0
kestore=dblarr(npts)
for i=1,npts DO BEGIN
 tflag=0
 ;print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(i,/ke)  
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
 IF ((tflag eq 0) and (ds.z[0] eq 0.45e4)) THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2, col=50
 IF ((tflag eq 0) and (ds.z[0] eq 0.5e4)) THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2, col=130
 IF ((tflag eq 0) and (ds.z[0] eq 0.55e4)) THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2, col=220
 
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

epts=n_elements(estore)
;STOP

;save, kestore, cstore, /compress, /verbose, filename=nom+'../IDL/colorstore_test.sav'
;restore, filename=nom+'../IDL/colorstore_test.sav',/verbose

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
cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
title='Particle peak energy gain (eV)', /vertical, /ylog

;particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], zscl=1.0, /Mm, myxyrange=[-1.00,1.00], myzrange=[-1.75,2.75], /rel

;particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e3, $
particletrack, estore[1], floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e3, $
;/km, myxyrange=[-10.0,10.0], myzrange=[-17.5,27.5], /rel, /symb
/km, myxyrange=[-10.0,10.0], myzrange=[-10,50.0], /rel, /symb


for i=2,epts-1 DO BEGIN
 particletrack, estore[i], /op, floc=nom, xyzt, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e3, /km, /rel, /symb,myxyrange=[-10.0,10.0], myzrange=[-10.0,50.0]
endfor
;for i=2,npts-1 DO BEGIN
; particletrack, i, /op, floc=nom, xyzt, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e3, /km, /rel, /symb,myxyrange=[-10.0,10.0], myzrange=[-10.0,50.0]
;endfor

END
