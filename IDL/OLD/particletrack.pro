PRO particletrack, ii, n_c=n_c, tscl=tscl, op=op, floc=floc, xyzt
; plots a 3d particle track using xplot3d
; op flag tells whether to overplot
; floc says where the datafile is..


IF n_elements(op) eq 0 THEN op=0
IF n_elements(floc) eq 0 THEN floc="../Data/"
IF n_elements(n_c) eq 0 THEN n_c=23
IF n_elements(tscl) eq 0 THEN tscl=100

 thisPalette = Obj_New('IDLgrPalette')
 thisPalette->LoadCT, 34

 oOrb = OBJ_NEW('orb', COLOR=[0, 255 ,0]) 
 oOrb2 = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 
 oOrb->Scale, .04, .04, .25 
 oOrb2->Scale, .04, .04, .25 
 oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol is green orb for start point
 oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 is red orb for stop point


; Tscl=100
 ;n_c=23

 filename=string(floc, ii,format='(a,"RV",i8.8,".dat")')
 tstring=string(floc, ii,format='("wc -l ",a,"RV",i8.8,".dat")')
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
   
   IF (op EQ 0) THEN BEGIN
    xplot3d, xx, yy, zz, xrange=[-1,1], yrange=[-1,1], zrange=[-6,6], color=[0,0,0] 
   ENDIF ELSE BEGIN
    xplot3d, xx, yy, zz, color=[0,0,0], /OVERPLOT
   ENDELSE

;tricking the 3d graphics - xplot3d needs an array, as it puts symbols at the end vertices of lines.
 ; if I feed it the start and end points as a two element array, line starts and ends at same vertex!
 ;ta=[ [x[0],x[n_data-1]],[y[0],y[n_data-1]],[z[0],z[n_data-1]]]
 ts=[ [x[0],x[0]],[y[0],y[0]],[z[0],z[0]]]
 te=[ [x[n_data-1],x[n_data-1]],[y[n_data-1],y[n_data-1]],[z[n_data-1],z[n_data-1]]]

   XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
   XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT


xyzt=[[xx],[yy],[zz],[tt]]


END
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
