;KG: changed to use 23 element RV

n=13475
nc=23
nplot=n-1
A=fltarr(nc,n)
t=fltarr(n)
x=fltarr(n)
y=fltarr(n)
z=fltarr(n)
vpar=fltarr(n)
vpar2=fltarr(n)
vperp2=fltarr(n)
epar=fltarr(n)
eperp=fltarr(n)
edrift=fltarr(n)
openr,lun,'RV.dat',/get_lun
readf,lun,A
free_lun,lun
t=transpose(A(0,*))
x=transpose(A(1,*))
y=transpose(A(2,*))
z=transpose(A(3,*))
vpar=transpose(A(4,*))
vpar2=vpar*vpar
vperp2=transpose(A(6,*))
epar  =transpose(A(10,*))
eperp =transpose(A(11,*))
edrift=transpose(A(12,*))
set_plot,'ps'
device,filename='et_drift_2.eps',/encapsulated
plot,t,(epar(0:nplot)+eperp(0:nplot)+edrift(0:nplot))/1000., psym=3, $
 xrange=[0,100], $
 xtitle='time (s)', ytitle='particle energy (keV)', charsize=2, $
 xmargin=[7.3,2],ymargin=[3.5,1]
oplot,t,epar(0:nplot)+eperp(0:nplot)+edrift(0:nplot)
device,/close
;wait,10
;
;set_plot,'ps'
;device,filename='e_drift_2_new.eps',/encapsulated
;plot,t,epar(0:nplot), psym=4, $
; xrange=[0,200], $
; xtitle='time (s)', ytitle='energy (eV)', charsize=1.5
;oplot,t*100.-105.,epar(0:nplot)
;oplot,t*100.-105.,eperp(0:nplot),psym=2
;oplot,t*100.-105.,eperp(0:nplot)
;device,/close

set_plot,'X'
end
