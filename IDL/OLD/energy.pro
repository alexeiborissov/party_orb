;;;we use this routine to produce plots of particle energy
Tscl=100
n_c=23 ;was 21

;get the length of the file to read in
;probibly a more efficient way to do this, but it works!
spawn,'wc -l RV.dat',res
n_data=1UL
dum='blah'
reads,res,n_data,dum

A=fltarr(n_c,n_data)
openr,lun,'RV.dat',/get_lun
readf,lun,A
Free_lun,lun
t=transpose(A(0,*))

loadct,5
et=transpose(A(13,*)+A(14,*)+A(15,*)) ;new RV file (not RV1) doesn't give all required elements. A(15,*) may not be right
plot,t,et/1000., $ ;xrange=[5,7],yrange=[14.6,16.4], $
xtitle='time (s)',ytitle='particle energy (keV)', charsize=2., ystyle=1, $
xmargin=[7.3,2],ymargin=[3.5,1],color=55

;set_plot,'ps'
;device,filename='energy_5to7sec_new.eps',/encapsulated
;plot,t,et/1000.,xrange=[5,7],yrange=[14.6,16.4], $
;xtitle='time (s)',ytitle='particle energy (keV)', charsize=2., ystyle=1, $
;xmargin=[7.3,2],ymargin=[3.5,1]
;device,/close
set_plot,'X'

; EXAMPLE OF LABELS
; ytickv=[1.5*10^4,1.55*10^4,1.6*10^4],  $
; ytickname=['1.5!U.!N10!U4!N','1.55!U.!N10!U4!N','1.6!U.!N10!U4!N']


; IF YOU WANT TO PLOT PARALLEL AND GYRATIONAL ENERGY
;eparall=transpose(A(10,*))
;erot   =transpose(A(11,*))
;plot,t,erot
;oplot,t,eparall

end
