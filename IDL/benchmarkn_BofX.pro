;common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 39, /silent
tvlct, r, g, b, /get

npts=32*32*1
nom="Data/"
;nom="../rcode/Data/"
;nom="../teq0/Data/"

;STOP

testparticle=npts/2-sqrt(npts)/2

particletrack, testparticle, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /mm, myxyrange=[-10,10], myzrange=[0,60]
for i=1,npts,16 DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[0,0,0], zscl=1e6, /mm, myxyrange=[-10,10], myzrange=[0,60]
endfor

;particletrack, npts/2, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /bsymb

ds=getndata(testparticle, /ke,/rdotperp, /fields, /vpar)
m=9.1093826e-31
aq=1.602176530e-19

myvpar=sqrt(4.0d0*aq/m)/sqrt(2.0d0)

myrdotperp=0.5d0*m*(ds.vpar(0)*ds.vpar(0))/aq/ds.b(2,0)*1e-7

myydisp=0.5d0*m*(ds.vpar(0)*ds.vpar(0))/aq/ds.b(2,0)*ds.t*1e-7

thisLetter = "143B
char_gamma = '!4' + String(thisLetter) + '!X'
thisLetter2 = "154B
char_mu = '!4' + String(thisLetter2) + '!X'


!p.background=255
window, 0, ysize=800
!p.multi=[0,1,3]
plot, ds.t, ds.ke, psym=-2, thick=2, xtitle='t (s)', ytitle='ke (eV)', charsize=3, yr=[1.95,2.05], xr=[0,100],$
title=string(testparticle,format='("non-relativistic: kinetic energy vs time, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2
oplot, [0,100], [2,2], col=240, linestyle=2, thick=2
;legend,['numerical','analytical'],psym=[-2,0], color=[0,240], linestyle=[0,2],/right, /top, textcolors=0, charsize=2, outline_color=0

;plot, ds.t, ds.rdotperp*1e6, psym=-2, thick=2, xtitle='t (s)', ytitle='d|R!d!9x!n!3|/dt ('+char_mu+'ms!e-1!n)', yr=[52.5,53.5],xr=[0,100],$
plot, ds.t, ds.rdotperp*1e6, psym=-2, thick=2, xtitle='t (s)', ytitle='d|R!d!9x!n!3|/dt ('+char_mu+'ms!e-1!n)', yr=[9,11],xr=[0,100],$
charsize=3, title=string(testparticle,format='("non-rel: d|R!d!9x!n!3|/dt vs time, particle no ",i2)'), col=0,  xthick=2, ythick=2, symsize=2
oplot, [0,100], [myrdotperp,myrdotperp]*1e6, col=240, linestyle=2, thick=2

plot, ds.t, ds.vpar/1e6, psym=-2, thick=2, xtitle='t (s)', ytitle='V!d!9#!n!3 (Mms!e-1!n)', yr=[0.587,0.597],xr=[0,100],$
charsize=3, title=string(testparticle,format='("non-rel: V!d!9#!n!3 vs time, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2
oplot, [0,100], [myvpar,myvpar]/1e6, linestyle=2, col=240, thick=2

;WRITE_PNG, "nonrel_p6_b0exp.png", TVRD(/TRUE)

!p.multi=0   
window, 2
plot, ds.t, 1e3*(ds.y-ds.y(0)), psym=-2, thick=2, xtitle='t (s)', ytitle='y-displacement (mm)', charsize=2, xr=[0,100], yr=[0,1],$
title=string(testparticle,format='("non-relativistic: y-displacement, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2
oplot, ds.t, myydisp*1e3, linestyle=2, col=240, thick=2

;WRITE_PNG, "nonrel_p6_b0exp_displ.png", TVRD(/TRUE)

!p.background=0 

;vpar='V!d!9#!n!3'
;vperp='V!d!9x!n!3'
END
