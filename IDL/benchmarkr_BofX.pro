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

particletrack, testparticle, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /mm, myxyrange=[-10,10], myzrange=[0,60], /rel
for i=1,npts,16 DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[0,0,0], zscl=1e6, /mm, /rel
endfor
particletrack, testparticle, /op, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /mm, /rel
;particletrack, npts/2, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /bsymb

ds=getrdata(testparticle, /ke,/rdotperp, /vpar, /rel, /fields)

aq=1.602176530e-19
m=9.1093826e-31
c=2.99792458e8

mygammaminus1=2.0d0*aq/m/c/c
myupar=sqrt(((mygammaminus1+1.0d0)*(mygammaminus1+1.0d0)-1.0d0))*c/sqrt(2.0d0)

thisLetter = "143B
char_gamma = '!4' + String(thisLetter) + '!X'
str1="relativistic: "+char_gamma+'-1 vs time, particle no.'
thisLetter2 = "154B
char_mu = '!4' + String(thisLetter2) + '!X'
;myupar=5.930975493e5


myvperp=sqrt(c*c-c/(mygammaminus1+1.0d0)*c/(mygammaminus1+1.0d0))/sqrt(2)
;myrdotperp=0.5d0*m*(myvperp*myvperp)/aq*(mygammaminus1+1.0d0)/ds.b(2,0)*1e-7
myrdotperp=0.5d0*m*(ds.vpar(0)*ds.vpar(0))/aq*(mygammaminus1+1.0d0)/ds.b(2,0)*1e-7

myydisp=0.5d0*m*(ds.vpar(0)*ds.vpar(0))/aq/(mygammaminus1+1.0d0)/ds.b(2,0)*ds.t*1e-7

;mygammaminus1=3.908489100e-6

!p.background=255
window, 0, ysize=800
!p.multi=[0,1,3]
;plot, ds.t, ds.ek, psym=-2, thick=2, xtitle='t (s)', ytitle='ke (eV)', yr=[1,3], xr=[0,100], charsize=3, title=string(testparticle,format='("rel kinetic energy vs time, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2
plot, ds.t, ds.gamma-1.0d0, psym=-2, thick=2, xtitle='t (s)', ytitle=char_gamma+'-1', xr=[0,100], yr=[3.85e-6,3.95e-6],charsize=3,$
 title=string(str1,testparticle,format='(A,i2)'), col=0, xthick=2, ythick=2, symsize=2, xmargin=[13,3]
oplot, [0,100], [mygammaminus1,mygammaminus1], linestyle=2, thick=2, col=240
legend,['numerical','analytical'],psym=[-2,0], color=[0,240], linestyle=[0,2],/right, /bottom, textcolors=0, charsize=2, outline_color=0

;plot, ds.t, ds.rdotperp, psym=-2, thick=2, xtitle='t (s)', ytitle='d|R!d!9x!n!3|/dt (ms!e-1!n)', xr=[0,100], yr=[-1,1],charsize=3, $
;title=string(testparticle,format='("rel: d|R!d!9x!n!3|/dt vs time, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2, xmargin=[13,3]
;oplot, [0,100], [0.0d0,0.0d0], linestyle=2, thick=2, col=240

plot, ds.t, ds.rdotperp*1e6, psym=-2, thick=2, xtitle='t (s)', ytitle='d|R!d!9x!n!3|/dt ('+char_mu+'ms!e-1!n)', yr=[9.5,10.5], xmargin=[13,3],$
charsize=3, title=string(testparticle,format='("rel: d|R!d!9x!n!3|/dt vs time, particle no ",i2)'), col=0,  xthick=2, ythick=2, symsize=2, xr=[0,100]
oplot, [0,100], [myrdotperp,myrdotperp]*1e6, col=240, linestyle=2, thick=2

plot, ds.t, ds.upar/ds.gamma*1e-6, psym=-2, thick=2, xtitle='t (s)', ytitle='u!d!9#!n!3 (Mms!e-1!n)', xr=[0,100], yr=[0.587,0.597],charsize=3, $
title=string(char_gamma,testparticle,format='("rel: u!d!9#!n!3 (=",A,+"v!d!9#!n!3) vs time, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2, xmargin=[13,3]
oplot, [0,100], [myupar,myupar]*1e-6, linestyle=2, thick=2, col=240

;WRITE_PNG, "rel_p6_b0exp.png", TVRD(/TRUE)


!p.multi=0   
window, 2
plot, ds.t, 1e3*(ds.y-ds.y(0)), psym=-2, thick=2, xtitle='t (s)', ytitle='y-displacement (mm)', charsize=2, yr=[0,1],xr=[0,100],$
title=string(testparticle,format='("relativistic: y-displacement, particle no ",i2)'), col=0, xthick=2, ythick=2, symsize=2
oplot, ds.t, [myydisp,myydisp]*1e3, linestyle=2, col=240, thick=2

;WRITE_PNG, "rel_p6_b0exp_displ.png", TVRD(/TRUE)


!p.background=0 
!p.multi=0 

;vpar='V!d!9#!n!3'
;vperp='V!d!9x!n!3'
END
