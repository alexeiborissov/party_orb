
common somename, tt, FRon, len, an

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 1, /silent
tvlct, r, g, b, /get


npts=16*16*5
nom="../Data/"
;nom="../../Areduce5/Data/"
;nom="../../Areduce5/Data/"
;nom="../../Sept2/Data/"

grstore=dblarr(npts-1)

print, 'loading: ', nom

for i=1,npts-1 DO BEGIN
 print, i, npts, format='(i4,"/",i4)'
 ds=getdata(i,wkdir=nom,/gyro)  
 ;maxkegain=max(ds.ek-ds.ek(0),maxkeloc)
 maxgr=max(ds.gyror,maxgrloc)
 ;maxketime=ds.time(maxkeloc)
 grstore(i-1)=maxgr
endfor

;colstore=FIX(kestore/max(kestore)*254)
cstore=fix((grstore-min(grstore))/max(grstore-min(grstore))*254)
;
;particletrack, 1, n_c=32, tscl=0.0, lscl=1e7, op=0, floc=nom, xyzt, linecol=[r(cstore(0)),g(cstore(0)),b(cstore(0))]
;for i=2,npts-1,4 DO BEGIN
; particletrack, i, n_c=32, tscl=0.0, lscl=1e7, op=1, floc=nom, xyzt, linecol=[r(cstore(i-1)),g(cstore(i-1)),b(cstore(i-1))]
;endfor

;make a new window for the colorbar
cgwindow
cgColorbar, Divisions=6, Minor=5, Format='(e8.1)', Ticklen=-0.25, /window,$
range=[min(grstore),max(grstore)], title='peak gyroradius size (m)', /vertical

iniorbplot, 1, n_c=32, tscl=0.0, lscl=0.5e2, op=0, floc=nom, ocol=[r(cstore(0)),g(cstore(0)),b(cstore(0))], zscl=0.4e7
for i=2,npts-1 DO BEGIN
 iniorbplot, i, n_c=32, tscl=0.0, lscl=0.5e2, op=1, floc=nom, ocol=[r(cstore(i-1)),g(cstore(i-1)),b(cstore(i-1))], zscl=0.4e7
endfor

        
END
