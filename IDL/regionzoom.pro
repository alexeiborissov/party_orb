PRO regionzoom, ii, n_c=n_c, tscl=tscl, op=op, floc=floc, lscl=lscl, xyzt, symb=symb, linecol=linecol, zscl=zscl, bouncesymb=bouncesymb, mytitle=mytitle, zoomy=zoomy, zoomrange=zoomrange
;@getdata
; plots a 3d particle track using xplot3d THROUGH A GIVEN REGION!
; op flag tells whether to overplot
; floc says where the datafile is..

IF n_elements(op)   	eq 0 THEN op=0    	    	    ; default action - do not overplot
IF n_elements(floc) 	eq 0 THEN floc="../Data/"	    ; default location of data
IF n_elements(n_c)  	eq 0 THEN n_c=32 	    	    ; default no of columns in datafile
IF n_elements(tscl) 	eq 0 THEN tscl=100.0    	    ; default timescale
IF n_elements(lscl) 	eq 0 THEN lscl=1e7	    	    ; default lengthscale
IF n_elements(linecol)	eq 0 THEN linecol=[0,0,0]   	    ; default track color
IF n_elements(zscl) 	eq 0 THEN zscl=lscl     	    ; default z scale
IF n_elements(mytitle)	eq 0 THEN mytitle='zoomy woomy'

 thisPalette = Obj_New('IDLgrPalette')
 thisPalette->LoadCT, 34
 dims=[1.,1.,6.]
 IF keyword_set(zoomy) THEN mysymsize=reform(zoomrange[1,*]-zoomrange[0,*])*0.1*dims ELSE mysymsize=[lscl,lscl,zscl]*0.1*dims  
 
 oOrb = OBJ_NEW('orb', COLOR=[0, 255 ,0]) 
 oOrb2 = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 
 IF keyword_set(zoomy) THEN oOrb->Scale, .02*0.5*(zoomrange[1,0]-zoomrange[0,0]), .02*0.5*(zoomrange[1,1]-zoomrange[0,1]), .125*0.5*(zoomrange[1,2]-zoomrange[0,2]) $
    	    	    	ELSE oOrb->Scale, .02*lscl, .02*lscl, .125*zscl 
 IF keyword_set(zoomy) THEN oOrb2->Scale, .02*0.5*(zoomrange[1,0]-zoomrange[0,0]), .02*0.5*(zoomrange[1,1]-zoomrange[0,1]), .125*0.5*(zoomrange[1,2]-zoomrange[0,2]) $
    	    	    	ELSE oOrb2->Scale, .02*lscl, .02*lscl, .125*zscl 
 oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol is green orb for start point
 oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 is red orb for stop point
 oModel=obj_new('IDLgrModel')
 oModel->add,symbol_obj(2,[0,0,255],[0,0,0])
 oTet = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oTet) THEN oTet->SetProperty, Size=mysymsize 
 ;oOrb3 = OBJ_NEW('orb', COLOR=[0,0,255]) 
 ;oOrb3->Scale, .02*lscl, .02*lscl, .125*zscl 
 ;oSymbol3 = OBJ_NEW('IDLgrSymbol', oOrb3) ;osymbol2 is red orb for stop point
 
 ds=getdata(ii, wkdir=floc, /vpar, /fields)
 newx=ds.x
 newy=ds.y
 newz=ds.z
 newt=ds.t
 vpar=ds.vpar

 IF keyword_set(zoomy) THEN zr=zoomrange ELSE zr=[[[-1,1]*lscl],[[-1,1]*lscl],[[-6,6]*zscl]]
   ;zr=[[-1e6,1e6],[-1e6,1e6],[-1e6,1e6]]
  
  ex=where((newx gt zr[0,0]) and (newx lt zr[1,0]))
  newx=newx(ex)
  newy=newy(ex)
  newz=newz(ex)
  newt=newt(ex)
  newvp=vpar(ex)
 
  ey=where((newy gt zr[0,1]) and (newy lt zr[1,1]))
  newx=newx(ey)
  newy=newy(ey)
  newz=newz(ey)
  newt=newt(ey)
  newvp=newvp(ey)
 
  ez=where((newz gt zr[0,2]) and (newz lt zr[1,2]))
  newx=newx(ez)
  newy=newy(ez)
  newz=newz(ez)
  newt=newt(ez)
  newvp=newvp(ez)
 
  n_data=n_elements(newt)
    
 ; Plot the random points in 3D space with a filled circle shape.
  tt=congrid(newt,1024)
  xx=congrid(newx,1024)
  yy=congrid(newy,1024)
  zz=congrid(newz,1024)
  tcolors = BYTSCL(tt)
  ;cgPlotS, xx, yy, zz, PSYM=2, COLOR=tcolors, SYMSIZE=1.5, /T3D
   
  IF (op EQ 0) THEN BEGIN
   xplot3d, newx, newy, newz, xrange=zr[*,0], yrange=zr[*,1], zrange=zr[*,2], color=linecol, title=mytitle
  ENDIF ELSE BEGIN
   xplot3d, newx, newy, newz, color=linecol, /OVERPLOT
  ENDELSE

 ;tricking the 3d graphics - xplot3d needs an array, as it puts symbols at the end vertices of lines.
 ; if I feed it the start and end points as a two element array, line starts and ends at same vertex!
 ;ta=[ [x[0],x[n_data-1]],[y[0],y[n_data-1]],[z[0],z[n_data-1]]]
 ;ts=[ [x[0],x[0]],[y[0],y[0]],[z[0],z[0]]]
 ;te=[ [x[n_data-1],x[n_data-1]],[y[n_data-1],y[n_data-1]],[z[n_data-1],z[n_data-1]]]

 ;IF keyword_set(symb) THEN XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
 ;IF keyword_set(symb) THEN XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT

 IF keyword_set(bouncesymb) THEN BEGIN
 
  ;malocs=extrema(reform(sqrt(ds.B[0,*]*ds.B[0,*]+ds.B[1,*]*ds.B[1,*]+ds.B[2,*]*ds.B[2,*])),/max_only)
  ;milocs=extrema(abs(ds.vpar),/min_only)   	    	; there can be vpar minima without necessarily being bounce points
  gtz=(where(newvp*shift(newvp,1) lt 0))
  IF (n_elements(gtz) ge 2) then gtz=gtz[1:*]
  
  ;IF n_elements(malocs) ne n_elements(milocs) THEN PRINT, 'ERROR, B(max) and Vpar(min) do not agree!'
  IF n_elements(gtz) ge 2 THEN BEGIN
   FOR j=0,n_elements(gtz)-1 DO BEGIN
    tel=[ [newx[gtz[j]],newx[gtz[j]]], [newy[gtz[j]],newy[gtz[j]]], [newz[gtz[j]],newz[gtz[j]]]]
    ;XPLOT3D, tel[*,0], tel[*,1], tel[*,2], COLOR=[0,0,0], NAME='bounce', SYMBOL=oTet, THICK=2, /OVERPLOT  
   ENDFOR
  ENDIF
 ENDIF
 

xyzt=[[newx],[newy],[newz],[newt]]


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
