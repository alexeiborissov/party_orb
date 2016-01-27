;PRO SCATTER3D, PostScript=postscript


Tscl=100
n_c=23

;npts=9*9*11
npts=21
;filename=string(1,format='("../Data/RV",i8.8,".dat")')
;tstring=string(1,format='("wc -l ../Data/RV",i8.8,".dat")')

thisPalette = Obj_New('IDLgrPalette')
thisPalette->LoadCT, 34

oOrb = OBJ_NEW('orb', COLOR=[0, 255 ,0]) 
oOrb2 = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 
oOrb->Scale, .04, .04, .25 
oOrb2->Scale, .04, .04, .25 
oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol is green orb for start point
oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 is red orb for stop point
 
 
;postscript=0

;IF Keyword_Set(postscript) THEN PS_Start, FILENAME='particle_test1.ps',FONT=1, Charsize=3.0 ELSE  window, 0
; ; Load a color table and create colors for the scatterplot.
; cgLoadCT, 33

;  cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[-1,1], rotz=60,$
;   YRANGE=[-1,1], ZRANGE=[-6, 6], XSTYLE=1, $
;   YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
;   POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
;   XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1
;   cgAXIS, XAXIS=1, /T3D, CHARSIZE=2.0, title='x'
;   cgAXIS, YAXIS=1, /T3D, CHARSIZE=2.0, title='y'
;
;   phi = Findgen(32) * (!PI * 2 / 32.)
;   phi = [ phi, phi(0) ]
 

;for i=1,npts DO BEGIN
for i=24,24 DO BEGIN
 
 particletrack, i, op=0
 
 
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
 t=100*transpose(A(0,*))
 x=transpose(A(1,*))
 y=transpose(A(2,*))
 z=transpose(A(3,*))

   
  ; Set the 3D coordinate space with axes.
 ; cgSurf, DIST(5), /NODATA, /SAVE, XRANGE=[-1,1], $
 ;  YRANGE=[-1,1], ZRANGE=[-6, 6], XSTYLE=1, $
 ;  YSTYLE=1, ZSTYLE=1, CHARSIZE=2.0, $
 ;  POSITION=[0.1, 0.1, 0.95, 0.95, 0.1, 0.95], $
 ;  XTICKLEN=1, YTICKLEN=1, XGRIDSTYLE=1, YGRIDSTYLE=1
 ;  cgAXIS, XAXIS=1, /T3D, CHARSIZE=2.0
 ;  cgAXIS, YAXIS=1, /T3D, CHARSIZE=2.0
    
   ; Plot the random points in 3D space with a filled circle shape.
   tt=congrid(t,1024)
   xx=congrid(x,1024)
   yy=congrid(y,1024)
   zz=congrid(z,1024)
   tcolors = BYTSCL(tt)
   ;cgPlotS, xx, yy, zz, PSYM=2, COLOR=tcolors, SYMSIZE=1.5, /T3D
   xplot3d, xx, yy, zz, xrange=[-1,1], yrange=[-1,1], zrange=[-6,6], color=[0,0,0] 

;tricking the 3d graphics - xplot3d needs an array, as it puts symbols at the end vertices of lines.
; if I feed it the start and end points as a two element array, line starts and ends at same vertex!
;ta=[ [x[0],x[n_data-1]],[y[0],y[n_data-1]],[z[0],z[n_data-1]]]
ts=[ [x[0],x[0]],[y[0],y[0]],[z[0],z[0]]]
te=[ [x[n_data-1],x[n_data-1]],[y[n_data-1],y[n_data-1]],[z[n_data-1],z[n_data-1]]]

   XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
   XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
 
 wait, 0.5  
   
   ; Connect the data points to the XY plane of the plot.
 ;FOR j=0,n_data-1 DO cgPlotS, [x(j), x(j)], [y(j), y(j)], [-6, z(j)], COLOR=zcolors(j), /T3D   
 ;wait, 0.5    
endfor
    ; Close the PostScript file and clean-up, if required.
;    IF Keyword_Set(postscript) THEN PS_End, /PNG

;a[0] 	    - time
;a[1:3]     - R
;a[4]	    - Vpar
;a[5]	    - muB.B
;a[6]	    - sum((dRdt-Vpar*bb)^2)
;a[7:9]	    - Escl*E
;a[10:12]   - Bscl*B
;a[13:15]   - particle energy?
        
END
