
common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1


npts=16*16*5
nom="../../Areduce/Data/"
;nom="../Data/"
;nom="../../SeptNoEfield/Data/"
;nom="../../Sept2/Data/"
;particletrack, 1, n_c=32, tscl=0.0, lscl=1e7, op=0, floc=nom, xyzt
particletrack, 1, n_c=32, tscl=0.0, lscl=1e6, op=0, floc=nom, xyzt, zscl=1e6, /bouncesymb, mytitle='particle tracks: a=0.1*z0'

for i=2,npts-4,4 DO BEGIN
;for i=14,npts-1,9 DO BEGIN
 particletrack, i, n_c=32, tscl=0.0, lscl=1e6, op=1, floc=nom, xyzt , zscl=1e6, /bouncesymb
;print, reform(xyzt[1023,3])*100.;
endfor

;FRon=0
;;npts=8*8
;nom="../../EqONLY/Data/"
;particletrack, 1, n_c=32, tscl=0.0, op=0, floc=nom, xyzt;
;
;for i=2,npts-1 DO BEGIN
;;;for i=14,npts-1,9 DO BEGIN
; particletrack, i, n_c=32, tscl=0.0, op=1, floc=nom, xyzt
;endfor

        
END
