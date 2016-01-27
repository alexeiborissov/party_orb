@ODEINT
@mk_vector
heap_gc
common somename, tt, FRon, an

eps=0.0000001
h1=0.01
string1='derivs'
len=1e6
FRon=1
tt=0.0

particletrack, 10, n_c=32, tscl=0.0, lscl=1e8, op=0, floc=nom, xyzt, linecol=[255,0,0], zscl=1e7, /bouncesymb

 FOR it =0,1023,2 DO BEGIN  
   ystart=reform(xyzt[it,0:2])
   tt=reform(xyzt[it,3]);*100./95.0
   yp=ystart
   ODEINT, ystart, eps, h1, string1, yp
   xplot3d, yp(0,*), yp(1,*), yp(2,*), COLOR=[0,0,255], /OVERPLOT
  ENDFOR

END
