;common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 39, /silent
tvlct, r, g, b, /get

npts=4*4*1
nom="Data/"
;nom="../rcode/Data/"
;nom="../teq0/Data/"

;STOP

testparticle=npts/2-sqrt(npts)/2

particletrack, testparticle, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /mm, myxyrange=[-10,10], myzrange=[0,60]
for i=1,npts DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[0,0,0], zscl=1e6, /mm, myxyrange=[-10,10], myzrange=[0,60]
endfor

;particletrack, npts/2, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /bsymb

ds=getndata(testparticle, /ke,/rdotperp, /vpar)
m=9.1093826e-31
aq=1.602176530e-19

myvpar=sqrt(4.0d0*aq/m)/sqrt(2.0d0)



!p.background=255
window, 0, ysize=800
!p.multi=[0,1,3]
plot, ds.t, ds.ke, psym=-2, thick=2, xtitle='t (s)', ytitle='ke (eV)', xr=[0,100], yr=[1,3], charsize=3, $
title=string(testparticle,format='("non-relativistic: kinetic energy vs time, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2
oplot, [0,100], [2,2], col=240, linestyle=2, thick=2
legend,['numerical','analytical'],psym=[-2,0], color=[0,240], linestyle=[0,2],/right, /top, textcolors=0, charsize=2, outline_color=0

plot, ds.t, ds.rdotperp/1000000.d0, psym=-2, thick=2, xtitle='t (s)', ytitle='d|R!d!9x!n!3|/dt (Mms!e-1!n)', $
yr=[-1,1], xr=[0,100], charsize=3, title=string(testparticle,format='("non-rel: d|R!d!9x!n!3|/dt vs time, particle no ",i2)'), col=0,  xthick=2, ythick=2, symsize=2
oplot, [0,100], [0,0], col=240, linestyle=2, thick=2

plot, ds.t, ds.vpar/1e6, psym=-2, thick=2, xtitle='t (s)', ytitle='V!d!9#!n!3 (Mms!e-1!n)', yr=[0.587,0.597],$
xr=[0,100],charsize=3, title=string(testparticle,format='("non-rel: V!d!9#!n!3 vs time, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2
oplot, [0,100], [myvpar,myvpar]/1e6, linestyle=2, col=240, thick=2

WRITE_PNG, "nonrel_p6_b0.png", TVRD(/TRUE)

!p.background=0 
!p.multi=0   

;vpar='V!d!9#!n!3'
;vperp='V!d!9x!n!3'
END
