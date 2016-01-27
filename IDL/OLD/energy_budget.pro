n=8344
n=8876
n=13475 ;n is number of datapoints kept (see extra.dat)
A=fltarr(5,n)
openr,lun,'energy_budget.dat',/get_lun
;openr,lun,'energy_budget1e-7.dat',/get_lun
;openr,lun,'energy_budget1e-8.dat',/get_lun
readf,lun,A
free_lun,lun
t=  transpose(A(0,*))
te= transpose(A(2,*))
RHS1=transpose(A(3,*))
RHS2=transpose(A(4,*))
RHS=RHS1+RHS2
;
dte=deriv(t,te)
;set_plot,'ps'
;device,filename='energy_budget_4_rhs1+rhs2.eps',/encapsulated
plot,t,dte/1000.,psym=4;,xrange=[40,40.5];,yrange=[0.65,0.7] ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;,yrange=[0,3.2], ystyle=1;,  $
;     charsize=1.5,xtitle='time (s)', $
;     ytitle='dE!Dk!N/dt (keV s!U-1!N)'
;oplot,t,dte/1000.
;oplot,t,RHS1/1000.,psym=5
;oplot,t,RHS2/1000.,psym=6
oplot,t,RHS/1000.,psym=5
;oplot,t,RHS/1000.
;oplot,t,replicate(0,n-1)
;device,/close

;set_plot,'ps'
;device,filename='energy_budget_new.eps',/encapsulated
plot,t,dte/1000.,psym=4,xrange=[5,7], $ 
     yrange=[-1,15], ystyle=1,       $
     charsize=1.5,xtitle='time (s)', $
     ytitle='dE!Dk!N/dt (keV s!U-1!N)'
oplot,t,RHS/1000,psym=7
oplot,t,RHS/1000
oplot,t,replicate(0,n)
;device,/close
set_plot,'X'

; Example of ticks and labels
;     ytickv=[ 0,  5*10^3,10^4,1.5*10^4], $
;     ytickname=['0','5!U.!N10!U3!N','10!U4!N','1.5!U.!N10!U4!N']

end
