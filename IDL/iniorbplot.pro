PRO iniorbplot, ii, n_c=n_c, tscl=tscl, op=op, floc=floc, lscl=lscl, ocol=ocol, zscl=zscl, mytitle=mytitle, me=me, km=km, Mm=Mm, $
rel=rel, myxyrange=myxyrange, myzrange=myzrange, ax=ax, az=az, filewidth=filewidth, frac=frac
; plots a 3d orb at initial particle location using xplot3d, coloured in a specific way
; op flag tells whether to overplot
; floc says where the datafile is..

IF n_elements(op)   	eq 0 THEN op=0    	    	    ; default action - do not overplot
IF n_elements(floc) 	eq 0 THEN floc="Data/"	    ; default location of data
IF n_elements(n_c)  	eq 0 THEN n_c=32 	    	    ; default no of columns in datafile
IF n_elements(tscl) 	eq 0 THEN tscl=100.0    	    ; default timescale
;IF n_elements(lscl) 	eq 0 THEN lscl=1e6	    	    ; default lengthscale
IF n_elements(lscl) 	eq 0 THEN lscl=1e3	    	    ; default lengthscale
IF n_elements(ocol) 	eq 0 THEN ocol=[0,0,0]      	    ; default track color
IF n_elements(zscl) 	eq 0 THEN zscl=lscl         	    ; default track color
IF n_elements(mytitle)	eq 0 THEN mytitle='particle paths'
IF n_elements(myzrange)	eq 0 THEN myzrange=[-20,20]    	    ; default track color
IF n_elements(myxyrange)eq 0 THEN myxyrange=[-0.4,0.4]    	    ; default track color
IF n_elements(ax)   	eq 0 THEN ax=-60
IF n_elements(az)   	eq 0 THEN az=30
IF n_elements(filewidth) eq 0  THEN filewidth=800
IF n_elements(frac) 	eq 0  THEN frac=0.1

IF keyword_set(me) THEN BEGIN
  lscl=1e0
  rmult=10.0
  xt='x (m)'
  yt='y (m)'
  xyrat=10
ENDIF  
IF keyword_set(km) THEN BEGIN
  lscl=1e3
  rmult=1.0
  xt='x (km)'
  yt='y (km)'
  xyrat=0.1
  ;frac=0.1
ENDIF
IF keyword_set(Mm) THEN BEGIN
  lscl=1e6
  rmult=1.0
  xt='x (Mm)'
  yt='y (Mm)'
  xyrat=1.
  ;frac=0.1
ENDIF

 xylen=myxyrange[1]-myxyrange[0]
 zlen=myzrange[1]-myzrange[0]

 thisPalette = Obj_New('IDLgrPalette')
 thisPalette->LoadCT, 34

 oOrb = OBJ_NEW('orb', COLOR=ocol) 
 ;oOrb->Scale, .02, .02, .125
 oOrb->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol is green orb for start point

 IF keyword_set(rel) THEN BEGIN
  ds=getrdata(ii, wkdir=floc+'/DataR/')
  ;print, floc+'/DataR/'
 ENDIF ELSE BEGIN
  ds=getndata(ii, wkdir=floc+'/DataN/')
  ;print, floc+'/DataN/'
 ENDELSE
 
 x=ds.x
 y=ds.y
 z=ds.z
 t=ds.t
 
 ts=[ [x[0],x[0]],[y[0],y[0]],[z[0],z[0]]]
 
 IF (op EQ 0) THEN BEGIN
  XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/zscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, $ 
  xrange=myxyrange*rmult, yrange=myxyrange*rmult, zrange=myzrange, xtitle=xt, ytitle=yt, ztitle='z (Mm [L/10Mm])',$
  ax=ax, az=az, filewidth=filewidth
  ;, title=mytitle
 ENDIF ELSE BEGIN
  XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/zscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
 ENDELSE
 ;IF (op EQ 0) THEN BEGIN
 ; XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/zscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, $ 
 ; xrange=[-30.0,30.0], yrange=[-30.0,30.0], zrange=[-20.0,0.0], xtitle='x (m)', ytitle='y (m)', ztitle='z (Mm)';, title=mytitle
 ;ENDIF ELSE BEGIN
 ; XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/zscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
 ;ENDELSE
 
 
END
