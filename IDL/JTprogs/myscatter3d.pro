;myscatter3d - 

@ODEINT
common somename, tt, FRon, len

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1
len=1e7

npts=16*16*5
nom="../Data/"
particletrack, 1, n_c=32, tscl=0.0, op=0, floc=nom, xyzt

 FOR it =0,1023,2 DO BEGIN  
   ystart=reform(xyzt[it,0:2])
   tt=reform(xyzt[it,3]);*100./95.0
   yp=ystart
   ODEINT, ystart, eps, h1, string1, yp
   xplot3d, yp(0,*), yp(1,*), yp(2,*), COLOR=[0,0,255], /OVERPLOT
  ENDFOR
  
;FOR j=3,9,2 DO BEGIN
;particletrack, j, n_c=32, tscl=0.0, op=1, floc=nom, xyzt
;
; FOR it =0,1023,2 DO BEGIN  
;   ystart=reform(xyzt[it,0:2])
;   tt=reform(xyzt[it,3]);*100./95.0
;   yp=ystart
;   ODEINT, ystart, eps, h1, string1, yp
;   xplot3d, yp(0,*), yp(1,*), yp(2,*), COLOR=[0,0,255], /OVERPLOT
;  ENDFOR  
;ENDFOR

;FRon=0
;
;npts=8*8
;nom="../../nrcode_Eqonly/Data/"
;particletrack, 1, n_c=32, tscl=0.0, op=0, floc=nom, xyzt
;
; FOR it =0,1023,2 DO BEGIN  
;   ystart=reform(xyzt[it,0:2])
;   tt=reform(xyzt[it,3]);*100./95.0
;   ;tt=0.0
;   yp=ystart
;   ODEINT, ystart, eps, h1, string1, yp
;   xplot3d, yp(0,*), yp(1,*), yp(2,*), COLOR=[0,0,255], /OVERPLOT
;  ENDFOR


;for i=2,npts-1 DO BEGIN
;;for i=14,npts-1,9 DO BEGIN
; particletrack, i, n_c=23, tscl=100, op=1, floc=nom 
;endfor

        
END
