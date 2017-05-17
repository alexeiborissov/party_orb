;general program to plot some stuff
;datloc='/sraid5v1/alan/twoloops/twoloops/Data/'
datloc='../../laredata/twoloopsunstable/'
OUT=1
od='figs/'
ENERGYGRAPH=1

IF (ENERGYGRAPH) THEN BEGIN
 de=getenergy(datloc) 
 aspect_ratio=1.8
 myxs=30
 IF NOT(OUT) THEN myxs=myxs*50
 myys=myxs/aspect_ratio
 loadct, 39
 !p.background=255


 tote=de.en_b+de.en_ke+de.en_int

 IF (out) THEN BEGIN	; set up plot area
   set_plot,'ps'
   !p.font=0
   device, filename=od+'dEn_twoloops.eps', encapsulated=1, /helvetica, /color, BITS_PER_PIXEL=8
   device, xsize=myxs, ysize=myys
   mycharsize=1
   myth=4
   PRINT, 'output to: '+od+'dEn_twoloops.eps'
  ENDIF ELSE BEGIN
   window, 5, ysize=myys, xsize=myxs
   mycharsize=1
   myth=2
  ENDELSE 
 plot, de.time, (de.en_b-de.en_b(0))/totE*100., title='Change in Energy (relative to initial state)', xtitle='time', ytitle='dE (% of total energy)', yr=[-2,2], col=abs(!p.background-255),  thick=myth, xthick=myth, ythick=myth 
 oplot, de.time, (de.en_ke-de.en_ke(0))/totE*100., linestyle=2, col=abs(!p.background-255), thick=myth
 oplot, de.time, (de.en_int-de.en_int(0))/totE*100., linestyle=1, col=abs(!p.background-255), thick=myth
 legend, ['magnetic energy', 'kinetic energy', 'thermal energy'], linestyle=[0,2,1], color=abs(!p.background-255), /top, /left, textcolor=abs(!p.background-255),  thick=myth
 IF (out) THEN BEGIN
  device, /close
  set_plot, 'x'
  !p.font=-1
 ENDIF
ENDIF

END
