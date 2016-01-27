;common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 1, /silent
tvlct, r, g, b, /get


npts=16*16*1
nom="Data/"
;nom="../relcodev1.0/Data/"
;nom="../teq0/Data/"

print, 'loading: ', nom
ds=getndata(1,/ke)
inike=ds.ke(0)

;kestore=dblarr(npts-1)
;for i=1,npts-1 DO BEGIN
; print, i, npts, format='(i4,"/",i4)'
; ds=getndata(i,/ke)  
; maxke=max(ds.ke)
; ;maxketime=ds.time(maxkeloc)
; kestore(i-1)=maxke
;endfor
;cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

;save, kestore, cstore, /compress, /verbose, filename=nom+'../IDL/cstore_n.sav'
restore, filename=nom+'../IDL/cstore_n.sav',/verbose

;STOP

particletrack, 1, floc=nom, xyzt, lcol=[r(cstore(0)),g(cstore(0)),b(cstore(0))], zscl=1e6, /mm, myxyrange=[-100,100], myzrange=[-60,60]
for i=2,npts-1 DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[r(cstore(i-1)),g(cstore(i-1)),b(cstore(i-1))], zscl=1e6, /mm
endfor

STOP

;make a new window for the colorbar
cgwindow
cgColorbar, Divisions=6, Minor=5, Format='(e8.1)', Ticklen=-0.25, /window,$
range=[min(kestore)-inike,max(kestore)-inike], title='Gain in eV over initial value (200eV)', /vertical


;STOP
iniorbplot, 1, op=0, floc=nom, ocol=[r(cstore(0)),g(cstore(0)),b(cstore(0))], zscl=1e6, /me

for i=2,npts-1 DO BEGIN
print, i, npts, format='(i4,"/",i4)'
 iniorbplot, i, /op, floc=nom, ocol=[r(cstore(i-1)),g(cstore(i-1)),b(cstore(i-1))], zscl=1e6, /me
endfor

        
END
