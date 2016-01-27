PRO iniorb, ii, n_c=n_c, tscl=tscl, op=op, floc=floc, lscl=lscl, ocol=ocol, zscl=zscl, mytitle=mytitle
;@getdata
; plots a 3d orb at initial particle location using xplot3d, coloured in a specific way
; op flag tells whether to overplot
; floc says where the datafile is..

IF n_elements(op)   	eq 0 THEN op=0    	    	    ; default action - do not overplot
IF n_elements(floc) 	eq 0 THEN floc="../Data/"	    ; default location of data
IF n_elements(n_c)  	eq 0 THEN n_c=32 	    	    ; default no of columns in datafile
IF n_elements(tscl) 	eq 0 THEN tscl=100.0    	    ; default timescale
IF n_elements(lscl) 	eq 0 THEN lscl=1e6	    	    ; default lengthscale
IF n_elements(ocol) 	eq 0 THEN ocol=[0,0,0]      	    ; default track color
IF n_elements(zscl) 	eq 0 THEN zscl=lscl         	    ; default track color
IF n_elements(mytitle)	eq 0 THEN mytitle='particle paths'

 thisPalette = Obj_New('IDLgrPalette')
 thisPalette->LoadCT, 34

 oOrb = OBJ_NEW('orb', COLOR=ocol) 
 ;oOrb->Scale, .02, .02, .125
 oOrb->Scale, 1.0, 1.0, 1.0/3.0
 oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol is green orb for start point

 ds=getdata(ii, wkdir=floc)
 
 x=ds.x
 y=ds.y
 z=ds.z
 t=ds.t
 
 ts=[ [x[0],x[0]],[y[0],y[0]],[z[0],z[0]]]
 
 IF (op EQ 0) THEN BEGIN
  XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/zscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, $ 
  xrange=[-30.0,30.0], yrange=[-30.0,30.0], zrange=[-20.0,0.0], xtitle='x (m)', ytitle='y (m)', ztitle='z (Mm)';, title=mytitle
 ENDIF ELSE BEGIN
  XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/zscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
 ENDELSE
 
 
 
END
