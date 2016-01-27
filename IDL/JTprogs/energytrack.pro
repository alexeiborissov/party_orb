;common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 3, /silent
tvlct, r, g, b, /get


npts=16*16*5
nom="Data/"
nom="../teq0/Data/"

print, 'loading: ', nom
ds=getdata(1,/ke)
inike=ds.ek(0)

;kestore=dblarr(npts-1)
;for i=1,npts-1 DO BEGIN
; print, i, npts, format='(i4,"/",i4)'
; ds=getdata(i,wkdir=nom,/ke)  
; maxke=max(ds.ek)
; ;maxketime=ds.time(maxkeloc)
; kestore(i-1)=maxke
;endfor
;cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

;save, kestore, cstore, /compress, /verbose, filename=nom+'../IDL/cstore2.sav'
restore, filename=nom+'../IDL/cstore2.sav',/verbose

particletrack, 1, tscl=0.0, lscl=1e6, floc=nom, xyzt, linecol=[r(cstore(0)),g(cstore(0)),b(cstore(0))], zscl=1e6, /bsymb, mytitle='particle tracks: a=0.1*z0'
for i=2,npts-1,4 DO BEGIN
 particletrack, i, tscl=0.0, lscl=1e6, /op, floc=nom, xyzt, linecol=[r(cstore(i-1)),g(cstore(i-1)),b(cstore(i-1))], zscl=1e6, /bsymb
endfor

;make a new window for the colorbar
cgwindow
cgColorbar, Divisions=6, Minor=5, Format='(e8.1)', Ticklen=-0.25, /window,$
range=[min(kestore)-inike,max(kestore)-inike], title='Gain in eV over initial value (200eV)', /vertical


;STOP
iniorbplot, 1, tscl=0.0, lscl=1e0, op=0, floc=nom, ocol=[r(cstore(0)),g(cstore(0)),b(cstore(0))], zscl=1e6, mytitle='initial positions: a=0.1*z0'

for i=2,npts-1 DO BEGIN
print, i, npts, format='(i4,"/",i4)'
 iniorbplot, i, tscl=0.0, lscl=1e0, /op, floc=nom, ocol=[r(cstore(i-1)),g(cstore(i-1)),b(cstore(i-1))], zscl=1e6
endfor

        
END
