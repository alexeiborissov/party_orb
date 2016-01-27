;;;we use this routine to produce plots of particle energy
Tscl=100
n_c=29 ;was 23, before that 21

;st1='V!d!9x!n!3'; V_perp
;st2='V!d!9#!n!3'; V_par


npts=9
;filename=string(1,format='("../Data/RV",i8.8,".dat")')
;tstring=string(1,format='("wc -l ../Data/RV",i8.8,".dat")')

;postscript=0

;!p.multi=[0,1,2]
loadct, 39
IF Keyword_Set(postscript) THEN PS_Start, FILENAME='particle_test1.ps',FONT=1, Charsize=3.0 ELSE  window, 0, xsize=1000;, ysize=1000

for i=2,2 DO BEGIN
 filename=string(i,format='("../Data/RV",i8.8,".dat")')
 tstring=string(i,format='("wc -l ../Data/RV",i8.8,".dat")')
 spawn,tstring,res
 n_data=1UL
 dum='blah'
 reads,res,n_data,dum
 A=dblarr(n_c,n_data)
 openr,lun,filename,/get_lun
 readf,lun,A
 Free_lun,lun
 t=Tscl*transpose(A(0,*))
 ;modB=sqrt(transpose(A(1,*))*transpose(A(1,*))+transpose(A(2,*))*transpose(A(2,*))+transpose(A(3,*))*transpose(A(3,*)))
 Vpar=transpose(A(4,*))
  GBdrift=A[23:25,*]
  Cdrift=A[26:28,*]

plot,t,Vpar/max(abs(Vpar)), $
;xtitle='time (s)',ytitle='V!d!9#!n!3/max(|V!d!9#!n!3|) (ms!e-1!n)', charsize=2., ystyle=8, POSITION=[0.15, 0.15, 0.85, 0.95]
xtitle='time (s)',ytitle='V!d!9#!n!3/max(|V!d!9#!n!3|)', charsize=2., POSITION=[0.15, 0.15, 0.85, 0.95]
;oplot, [0,max(t)], [0,0], linestyle=1

;AXIS, YAXIS=1, ytitle='|B|/max(|B|)', YRANGE=[-1, 1], col=220, charsize=2,/SAVE
;AXIS, XAXIS=1, XRANGE=[0, max(t)], col=220, charsize=2,xtickname='',/SAVE

;oplot,t,modgradBdrift/max(abs(modgradBdrift)), col=220


;oplot,t,modB/max(abs(modB)), col=220;, $
;xtitle='time (s)',ytitle='|B| (T)', charsize=2., ystyle=1

;malocs=extrema(modB/max(abs(modB)))

;oplot, [t(malocs),t(malocs)], [(modB/max(abs(modB)))[malocs],(modB/max(abs(modB)))[malocs]]

;FOR j=0,n_elements(malocs)-1 DO oplot, [t[malocs[j]],t[malocs[j]]], [(modB/max(abs(modB)))[malocs[j]],(modB/max(abs(modB)))[malocs[j]]], col=180, psym=2, symsize=2
;FOR j=0,n_elements(malocs)-1 DO oplot, [t[malocs[j]],t[malocs[j]]], [-10,10], col=180, linestyle=1

;WAIT, 1.5


endfor
    IF Keyword_Set(postscript) THEN PS_End, /PNG

!p.multi=0
end
