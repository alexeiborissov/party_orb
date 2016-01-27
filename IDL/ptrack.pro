PRO particletrack, ii, n_c=n_c, tscl=tscl, op=op, floc=floc, lscl=lscl, xyzt, symb=symb, lcol=lcol, zscl=zscl, bsymb=bsymb, mytitle=mytitle, me=me, km=km, Mm=Mm, myxyrange=myxyrange, myzrange=myzrange, rel=rel, myls=myls, myfline=myfline, ax=ax, az=az, filewidth=filewidth
;@getdata
; plots a 3d particle track using xplot3d
; op flag tells whether to overplot
; floc says where the datafile is..
;@ODEINT
; run xplot3dJT to use angle and output options on xplot3d
common somename, tt, FRon, len, an



IF n_elements(op)   	eq 0 THEN op=0    	    	    ; default action - do not overplot
IF n_elements(floc) 	eq 0 THEN floc="Data/"	    	    ; default location of data
IF n_elements(n_c)  	eq 0 THEN n_c=32 	    	    ; default no of columns in datafile
IF n_elements(tscl) 	eq 0 THEN tscl=100.0    	    ; default timescale
IF n_elements(lscl) 	eq 0 THEN lscl=1e7	    	    ; default lengthscale
IF n_elements(lcol)	eq 0 THEN lcol=[0,0,0]   	    ; default track color
IF n_elements(zscl) 	eq 0 THEN zscl=lscl     	    ; default z scale
IF n_elements(mytitle)	eq 0 THEN mytitle=' '
IF n_elements(myxyrange)eq 0 THEN myxyrange=[-100,100]
IF n_elements(myzrange)	eq 0 THEN myzrange=[-60,60]
IF n_elements(myls)	eq 0 THEN myls=0
IF n_elements(ax)   	eq 0 THEN ax=-60
IF n_elements(az)   	eq 0 THEN az=30
IF n_elements(filewidth) eq 0  THEN filewidth=800

zt='z (Mm)'

IF keyword_set(me) THEN BEGIN
  lscl=1.0
  rmult=1.0
  xt='x (m)'
  yt='y (m)'
  xyrat=1.
ENDIF  
IF keyword_set(km) THEN BEGIN
  lscl=1e3
  rmult=1.0
  xt='x (km)'
  yt='y (km)'
  xyrat=1.
ENDIF
IF keyword_set(Mm) THEN BEGIN
  lscl=1e6
  rmult=1.0
  xt='x (Mm)'
  yt='y (Mm)'
  xyrat=1.
ENDIF
 frac=2.7
 xylen=DOUBLE(myxyrange[1]-myxyrange[0])
 zlen=DOUBLE(myzrange[1]-myzrange[0])

 thisPalette = Obj_New('IDLgrPalette')
 thisPalette->LoadCT, 34
 dims=[1.,1.,zlen/xylen]
 mysymsize=[lscl,lscl,zscl]*0.1*dims  
 
 oOrb = OBJ_NEW('orb', COLOR=[0, 255 ,0]) 
 oOrb2 = OBJ_NEW('orb', COLOR=[255, 0 ,0])
 oOrb3 = OBJ_NEW('orb', COLOR=[0, 0 ,255]) 	;blue orb
 ;oOrb3 = OBJ_NEW('orb', COLOR=[75, 0 ,130]) 	;indigo orb
 
 ;oOrb->Scale, 2., 2., 1.2 
 ;oOrb2->Scale, 2. , 2. , 1.2
 ;oOrb3->Scale, .4 , .4 , 1.2
 
 ;oOrb->Scale, 0.1*xyrat, 0.1*xyrat, 20.0/(6.*rmult)
 ;oOrb2->Scale, 0.1*xyrat, 0.1*xyrat, 20.0/(6.*rmult)
 ;oOrb3->Scale, 0.1*xyrat, 0.1*xyrat, 20.0/(6.*rmult)
 oOrb->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 oOrb2->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 ;oOrb3->Scale, frac*xyrat*1.1, frac*xyrat*1.1, zlen/xylen*xyrat*frac*1.1
 
 oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol is green orb for start point
 oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 is red orb for stop point
 
  ;oview = obj_new('IDLgrView')
 ;oview->add, oOrb
 ;idlgrwindow
 ;idlgrview, oOrb
 
 oModel=obj_new('IDLgrModel')
 oModel->add,symbol_obj(2,[0, 0 ,255],[0,0,0])	    ; blue tetrahedon
 ;oModel->add,symbol_obj(2,[75, 0 ,130],[0,0,0])    ; purple tetrahedon
 ;oModel->add,symbol_obj(3,[75, 0 ,130],[0,0,0])    ; purple pyramid
 oTet = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oTet) THEN oTet->SetProperty, Size=[frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac]*9.9
 ;xobjview,oTet
 ;STOP
 ;oOrb3 = OBJ_NEW('orb', COLOR=[0,0,255]) 
 ;oOrb3->Scale, .02*lscl, .02*lscl, .125*zscl 
 ;oTet = OBJ_NEW('IDLgrSymbol', oTet) ;osymbol3 is blue orb for bounce point
 ;oTet = OBJ_NEW('IDLgrSymbol', oOrb3) ;osymbol3 is blue orb for bounce point
 
 ;oTet->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 
 IF keyword_set(rel) THEN BEGIN
  ds=getrdata(ii, wkdir=floc+'/DataR/', /vpar, /fields)
 ENDIF ELSE BEGIN
  ds=getndata(ii, wkdir=floc+'/DataN/', /vpar, /fields)
 ENDELSE
 
 x=ds.x
 y=ds.y
 z=ds.z
 t=ds.t
 n_data=n_elements(t)
    
 ; Plot the random points in 3D space with a filled circle shape.
 ttt=congrid(t,1024)
 xx=congrid(x,1024)
 yy=congrid(y,1024)
 zz=congrid(z,1024)
 tcolors = BYTSCL(ttt)
 ;cgPlotS, xx, yy, zz, PSYM=2, COLOR=tcolors, SYMSIZE=1.5, /T3D
 
 IF keyword_set(myfline) THEN BEGIN ; if field lines are required, draw them first (in black)
   flc=[127,127,127]  	    ; field line and arrow color
   h1=0.01     	    	    ; step size?
   string1='derivs' 	    ; a string
   FRon=1	            ; flux ring ON?
   ystart=[x[0],y[0],z[0]]  ; where to start?
   yp=ystart
   ypb=ystart	    	    
   tt=0     	    	    ; time
   eps=0.0000001    	    ; error
   ODEINT, ystart, eps, h1, string1, yp     ; ODE INTEGRATION SUBROUTINE (FORWARD)
   ODEINT, ystart, eps, -h1, string1, ypb     ; ODE INTEGRATION SUBROUTINE (BACK)
   
   nd=size(yp)
   ; start and end points of forward field lines   
   fs=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
   fe=[ [yp[0,nd[2]-1],yp[0,nd[2]-1]],[yp[1,nd[2]-1],yp[1,nd[2]-1]],[yp[2,nd[2]-1],yp[2,nd[2]-1]]]
   oi=[ [yp[0,nd[2]-1]-yp[0,nd[2]-2]],[yp[1,nd[2]-1]-yp[1,nd[2]-2]],[yp[2,nd[2]-1]-yp[2,nd[2]-2]]] 
   ; start and end points of backwards field lines   
   ndb=size(ypb)
   fsb=[ [ypb[0,0],ypb[0,0]],[ypb[1,0],ypb[1,0]],[ypb[2,0],ypb[2,0]]]
   feb=[ [ypb[0,ndb[2]-1],ypb[0,ndb[2]-1]],[ypb[1,ndb[2]-1],ypb[1,ndb[2]-1]],[ypb[2,ndb[2]-1],ypb[2,ndb[2]-1]]]
   oib=[ [-ypb[0,ndb[2]-1]+ypb[0,ndb[2]-2]],[-ypb[1,ndb[2]-1]+ypb[1,ndb[2]-2]],[-ypb[2,ndb[2]-1]+ypb[2,ndb[2]-2]]] 
   ;orientation of forward arrow
   oModela=obj_new('IDLgrModel')
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=flc)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize*0.0001
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=flc)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize*0.0001
   IF (op EQ 0) THEN BEGIN
    xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=flc, xrange=myxyrange, yrange=myxyrange, $
    zrange=myzrange, xtitle=xt, ytitle=yt, ztitle=zt, title=mytitle, linestyle=ls, ax=ax, az=az, filewidth=filewidth
   ENDIF ELSE BEGIN
    xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=flc, /OVERPLOT
   ENDELSE
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=flc, /OVERPLOT
   xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=flc, NAME='atest', SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
   xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=flc, NAME='atestb', SYMBOL=oARb, THICK=2, /OVERPLOT   ; the initial field line arrow
   xplot3d, x/lscl, y/lscl, z/zscl, color=lcol, /OVERPLOT, linestyle=lsm, THICK=3     	    	    	    	    	    	; the actual particle track!
  ENDIF ELSE BEGIN ; else just draw the field lines (check whether overplot is on though)
   IF (op EQ 0) THEN BEGIN
    xplot3d, x/lscl, y/lscl, z/zscl, xrange=myxyrange, yrange=myxyrange, $
    zrange=myzrange, color=lcol, xtitle=xt, ytitle=yt, ztitle=zt, title=mytitle, linestyle=ls, THICK=3, ax=ax, az=az, filewidth=filewidth
   ENDIF ELSE BEGIN
    xplot3d, x/lscl, y/lscl, z/zscl, color=lcol, /OVERPLOT, linestyle=ls, THICK=3
   ENDELSE
  ENDELSE 

 ;tricking the 3d graphics - xplot3d needs an array, as it puts symbols at the end vertices of lines.
 ; if I feed it the start and end points as a two element array, line starts and ends at same vertex!
 ;ta=[ [x[0],x[n_data-1]],[y[0],y[n_data-1]],[z[0],z[n_data-1]]]
 ts=[ [x[0],x[0]],[y[0],y[0]],[z[0],z[0]]]
 te=[ [x[n_data-1],x[n_data-1]],[y[n_data-1],y[n_data-1]],[z[n_data-1],z[n_data-1]]]

 ;help, te
 IF keyword_set(symb) THEN XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/zscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
 IF keyword_set(symb) THEN XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/zscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT

 IF keyword_set(bsymb) THEN BEGIN
  ns=n_elements(ds.vpar)
  gtz=where(ds.vpar(0:ns-2)*ds.vpar(1:ns-1) le 0)
   IF (n_elements(gtz) eq 1) THEN BEGIN     ; if n_elements=1, there could still be no solutions (gtz=-1)
    IF gtz eq -1 THEN BEGIN
     PRINT, 'no bounce points found' 
    ENDIF ELSE BEGIN
     tel=[ [x[gtz],x[gtz]], [y[gtz],y[gtz]], [z[gtz],z[gtz]]]
     XPLOT3D, tel[*,0]/lscl, tel[*,1]/lscl, tel[*,2]/zscl, COLOR=[0,0,0], NAME='bounce', SYMBOL=oTet, THICK=2, /OVERPLOT
    ENDELSE
   ENDIF ELSE BEGIN
    FOR j=0,n_elements(gtz)-1 DO BEGIN
     tel=[ [x[gtz[j]],x[gtz[j]]], [y[gtz[j]],y[gtz[j]]], [z[gtz[j]],z[gtz[j]]]]
     XPLOT3D, tel[*,0]/lscl, tel[*,1]/lscl, tel[*,2]/zscl, COLOR=[0,0,0], NAME='bounce', SYMBOL=oTet, THICK=2, /OVERPLOT  
    ENDFOR
   ENDELSE
  ENDIF
  
  ;overplotting field lines?
  

xyzt=[[xx],[yy],[zz],[ttt]]


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
