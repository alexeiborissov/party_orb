PRO fwrap, op=op, myt=myt, mya=mya, rad=rad, midrad=midrad, zloc=zloc, intercept=intercept, isocurrent=isocurrent, $
isoccol=isoccol, ax=ax, az=az, filename=filename, transparent=transparent, filewidth=filewidth, frac=frac, centre=centre, arrows=arrows, flinethick=flinethick
;@ODEINT
;@xplot3dJT
;USES HACKED VERSION OF xplot3d, called xplot3dJT

common somename, tt, FRon, len, an

IF n_elements(op) eq 0     	THEN op=0
IF n_elements(myt) eq 0     	THEN myt=0.0d0 
IF n_elements(mya) eq 0     	THEN mya=0.1d0 
IF n_elements(rad) eq 0     	THEN rad=5000000.0d0
IF n_elements(midrad) eq 0   	THEN midrad=3000000.0d0
IF n_elements(zloc) eq 0    	THEN zloc=5e7
IF n_elements(intercept) eq 0 	THEN intercept=5e7
IF n_elements(isoccol) eq 0 	THEN isoccol=[200,0,0]
IF n_elements(ax) eq 0      	THEN ax=-60
IF n_elements(az) eq 0      	THEN az=30
IF n_elements(filewidth) eq 0   THEN filewidth=800
IF n_elements(filename) eq 0  	THEN filename='fieldwrappingtest.png'
IF n_elements(frac) eq 0   	THEN frac=100.0
IF n_elements(centre) eq 0  	THEN centre=0
IF n_elements(flinethick) eq 0	THEN flinethick=1

;type=strsplit(filename,".",/extract)
IF keyword_set(filename) and keyword_set(transparent) THEN BEGIN
 newfilenames=[type[0]+'_a.'+type[1],type[0]+'_b.'+type[1]]
ENDIF

 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]
 purple=[153,51,255]
 grey=[170,170,170]
 poobrown=[153,76,0]
 olive=[128,128,0]
 salmon=[255,160,122]
 lightsteelblue=[176,196,222]
 lavender=[230,230,250]
 rosybrown=[188,143,143]
 burlywood=[222,184,135]
 tann=[210,180,140]
 newpurple=[128,0,128]
 newpink=[255,153,255]
;col1=lightsteelblue
;col2=tann
col1=newpurple
col2=newpink

an=mya
me=0
km=0
Mm=1

myxrange=[-100e6,100e6]
myyrange=[-100e6,100e6]
myzrange=[-75e6,75e6]
lscl=1e6
zscl=1e6 
mytitle=' '

h1=0.00001     	    	    ; step size?
string1='derivs' 	    ; a string
FRon=1	            ; flux ring ON?
tt=myt     	    ; time
eps=0.000000001    	    ; error
zt='z (Mm)'

IF keyword_set(me) THEN BEGIN
  lscl=1.0
  rmult=1.0
  xt='x (m)'
  yt='y (m)'
ENDIF  
IF keyword_set(km) THEN BEGIN
  lscl=1e3
  rmult=1.0
  xt='x (km)'
  yt='y (km)'
ENDIF
IF keyword_set(Mm) THEN BEGIN
  lscl=1e6
  rmult=1.0
  xt='x (Mm)'
  yt='y (Mm)'
ENDIF

 xlen=DOUBLE(myxrange[1]-myxrange[0])
 ylen=DOUBLE(myyrange[1]-myyrange[0])
 zlen=DOUBLE(myzrange[1]-myzrange[0])

 thisPalette = Obj_New('IDLgrPalette')
 thisPalette->LoadCT, 34
 dims=[1.,1.,zlen/xlen]
 mysymsize=[lscl,lscl,zscl]*dims*1.5e-7*frac  
 
 oOrb4 = OBJ_NEW('orb', COLOR=[034, 139 ,034]) 	;Forest Green orb
 oOrb5 = OBJ_NEW('orb', COLOR=[124,252,0]) 	;Other green orb
 ;oOrb4 = OBJ_NEW('orb', COLOR=[255, 255 ,255]) 	;white orb
 ;oOrb5 = OBJ_NEW('orb', COLOR=[0, 0 ,0]) 	;black orb
 oOrb4->Scale, 0.01*frac, 0.01*frac, 0.01*zlen/xlen*frac
 oOrb5->Scale, frac, frac, zlen/xlen*frac
 oSy4 = OBJ_NEW('IDLgrSymbol', oOrb4) ;osymbol4 is black orb for null point 
 oSy5 = OBJ_NEW('IDLgrSymbol', oOrb5) ;osymbol4 is black orb for null point 
 
 nullu=[ [0.0d0,0.0d0],[0.0d0,0.0d0],[5e7,5e7]]
 nulld=[ [0.0d0,0.0d0],[0.0d0,0.0d0],[-5e7,-5e7]]   
 
 oModel=obj_new('IDLgrModel')
 oModel->add,symbol_obj(2,col1,[0,0,0])	    ; steelblue tetrahedon
 oTet = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oTet) THEN oTet->SetProperty, Size=[frac, frac, zlen/xlen*frac]
 oModel=obj_new('IDLgrModel')
 oModel->add,symbol_obj(2,col2,[0,0,0])	    ; grey tetrahedon
 oTet2 = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oTet2) THEN oTet2->SetProperty, Size=[frac, frac, zlen/xlen*frac]
 
print, 'making a movie..'
naz=00


FOR ii=0,naz DO BEGIN
 az=30+ii
 ;filename=type[0]+string(ii,format='(i3.3,".")')+type[1]
 print, ii, naz, format='("frame ", i3.3, "/", i3.3)' 
IF (keyword_set(filename) and keyword_set(transparent)) THEN ns=1 ELSE ns=0
FOR ij=0,ns DO BEGIN  ; for transparency (as outputting closes the xWindow)
 tmax=360
 ;tinc=24
 tinc=11
 ;flc=[127,127,127]  	    ; field line and arrow color (grey)
 IF (op EQ 0) THEN BEGIN
  xplot3d, nullu[*,0]/lscl, nullu[*,1]/lscl, nullu[*,2]/zscl, COLOR=flc, xrange=myxrange/lscl, yrange=myyrange/lscl, zrange=myzrange/lscl, $
  xtitle=xt, ytitle=yt, ztitle=zt, title=mytitle, linestyle=ls, SYMBOL=oTet2,ax=ax,az=az, filewidth=filewidth
 ENDIF ELSE BEGIN
  xplot3d, nullu[*,0]/lscl, nullu[*,1]/lscl, nullu[*,2]/zscl, COLOR=flc, linestyle=ls, SYMBOL=oTet2, /OVERPLOT
 ENDELSE
 
 xplot3d, nulld[*,0]/lscl, nulld[*,1]/lscl, nulld[*,2]/zscl, COLOR=flc, SYMBOL=oTet, /OVERPLOT
; fan rings: upper fan
  FOR i=0,tmax,tinc DO BEGIN 
   theta=DOUBLE(i)/360.*2.*!dpi
   ystart=[rad*cos(theta)-0.0d0,0.0,5.0*rad*sin(theta)+zloc]
   yp=ystart
   ypb=ystart
   ODEINT, ystart, eps, h1, string1, yp     ; ODE INTEGRATION SUBROUTINE (FORWARD)
   ODEINT, ystart, eps, -h1, string1, ypb     ; ODE INTEGRATION SUBROUTINE (BACK)
   iji=where((abs(yp[0,*]) lt myxrange[1]) and (abs(yp[1,*]) lt myyrange[1]) and (abs(yp[2,*]) lt myzrange[1]))
   yp=yp[*,iji]
   iji=where((abs(ypb[0,*]) lt myxrange[1]) and (abs(ypb[1,*]) lt myyrange[1]) and (abs(ypb[2,*]) lt myzrange[1]))
   ypb=ypb[*,iji]
   nd=size(yp)
   ;; start and end points of forward field lines   
   fs=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
   fe=[ [yp[0,nd[2]-1],yp[0,nd[2]-1]],[yp[1,nd[2]-1],yp[1,nd[2]-1]],[yp[2,nd[2]-1],yp[2,nd[2]-1]]]
   oi=[ [yp[0,nd[2]-1]-yp[0,nd[2]-2]],[yp[1,nd[2]-1]-yp[1,nd[2]-2]],[yp[2,nd[2]-1]-yp[2,nd[2]-2]]]
   ;; start and end points of backwards field lines   
   ndb=size(ypb)
   fsb=[ [ypb[0,0],ypb[0,0]],[ypb[1,0],ypb[1,0]],[ypb[2,0],ypb[2,0]]]
   feb=[ [ypb[0,ndb[2]-1],ypb[0,ndb[2]-1]],[ypb[1,ndb[2]-1],ypb[1,ndb[2]-1]],[ypb[2,ndb[2]-1],ypb[2,ndb[2]-1]]]
   oib=[ [-ypb[0,ndb[2]-1]+ypb[0,ndb[2]-2]],[-ypb[1,ndb[2]-1]+ypb[1,ndb[2]-2]],[-ypb[2,ndb[2]-1]+ypb[2,ndb[2]-2]]] 
   ;orientation of forward arrow
   oModela=obj_new('IDLgrModel')
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col2)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col2)
   oARb = obj_new('IDLgrSymbol', oModelb)
   ;xplot3d, temp[*,0]/lscl, temp[*,1]/lscl, temp[*,2]/zscl, COLOR=flc, linestyle=ls, SYMBOL=osy4, /OVERPLOT
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col2, thick=flinethick, /OVERPLOT
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col2, thick=flinethick, /OVERPLOT
   IF keyword_set(arrows) THEN xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col2, SYMBOL=oAR, THICK=2, /OVERPLOT  
   IF keyword_set(arrows) THEN xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col2, SYMBOL=oARb, THICK=2, /OVERPLOT   
  ENDFOR
 ;lower fan:
  FOR i=0,tmax,tinc DO BEGIN 
   theta=DOUBLE(i)/360.*2.*!dpi
   ystart=[0.0,rad*cos(theta)-0.0d0,5.0*rad*sin(theta)-zloc]
   yp=ystart
   ypb=ystart
   ODEINT, ystart, eps, h1, string1, yp     ; ODE INTEGRATION SUBROUTINE (FORWARD)
   ODEINT, ystart, eps, -h1, string1, ypb     ; ODE INTEGRATION SUBROUTINE (BACK)
   iji=where((abs(yp[0,*]) lt myxrange[1]) and (abs(yp[1,*]) lt myyrange[1]) and (abs(yp[2,*]) lt myzrange[1]))
   yp=yp[*,iji]
   iji=where((abs(ypb[0,*]) lt myxrange[1]) and (abs(ypb[1,*]) lt myyrange[1]) and (abs(ypb[2,*]) lt myzrange[1]))
   ypb=ypb[*,iji]
   nd=size(yp)
   ; start and end points of forward field lines   
   fs=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
   fe=[ [yp[0,nd[2]-1],yp[0,nd[2]-1]],[yp[1,nd[2]-1],yp[1,nd[2]-1]],[yp[2,nd[2]-1],yp[2,nd[2]-1]]]
   oi=[ [yp[0,nd[2]-1]-yp[0,nd[2]-2]],[yp[1,nd[2]-1]-yp[1,nd[2]-2]],[yp[2,nd[2]-1]-yp[2,nd[2]-2]]] 
   ;; start and end points of backwards field lines   
   ndb=size(ypb)
   fsb=[ [ypb[0,0],ypb[0,0]],[ypb[1,0],ypb[1,0]],[ypb[2,0],ypb[2,0]]]
   feb=[ [ypb[0,ndb[2]-1],ypb[0,ndb[2]-1]],[ypb[1,ndb[2]-1],ypb[1,ndb[2]-1]],[ypb[2,ndb[2]-1],ypb[2,ndb[2]-1]]]
   oib=[ [-ypb[0,ndb[2]-1]+ypb[0,ndb[2]-2]],[-ypb[1,ndb[2]-1]+ypb[1,ndb[2]-2]],[-ypb[2,ndb[2]-1]+ypb[2,ndb[2]-2]]] 
   ;orientation of forward arrow
   oModela=obj_new('IDLgrModel')
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col1)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col1)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col1, thick=flinethick, /OVERPLOT
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col1, thick=flinethick, /OVERPLOT
   IF keyword_set(arrows) THEN xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col1, NAME='atest', SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
   IF keyword_set(arrows) THEN xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col1, NAME='atestb', SYMBOL=oARb, THICK=2, /OVERPLOT   ; the initial field line arrow
  ENDFOR
   mysymsize=[lscl,lscl,zscl]*dims*1.8e-7*frac
   ystart=[1e6,0,-49e6]
   ypb=ystart
   ODEINT, ystart, eps, -h1, string1, ypb     ; ODE INTEGRATION SUBROUTINE (BACK)
   iji=where((abs(ypb[0,*]) lt myxrange[1]) and (abs(ypb[1,*]) lt myyrange[1]) and (abs(ypb[2,*]) lt myzrange[1]))
   ypb=ypb[*,iji]
   ; start and end points of backwards field lines   
   ndb=size(ypb)
   fsb=[ [ypb[0,0],ypb[0,0]],[ypb[1,0],ypb[1,0]],[ypb[2,0],ypb[2,0]]]
   feb=[ [ypb[0,ndb[2]-1],ypb[0,ndb[2]-1]],[ypb[1,ndb[2]-1],ypb[1,ndb[2]-1]],[ypb[2,ndb[2]-1],ypb[2,ndb[2]-1]]]
   epf1=[ypb[0,ndb[2]-1],ypb[1,ndb[2]-1],ypb[2,ndb[2]-1]]
   oib=[ [-ypb[0,ndb[2]-1]+ypb[0,ndb[2]-2]],[-ypb[1,ndb[2]-1]+ypb[1,ndb[2]-2]],[-ypb[2,ndb[2]-1]+ypb[2,ndb[2]-2]]] 
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col1)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col1, thick=flinethick+4, /OVERPLOT
   IF keyword_set(arrows) THEN xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col1, NAME='atestb', SYMBOL=oARb, THICK=2, /OVERPLOT   ; the initial field line arrow
   ystart=[-1e6,0,-49e6]
   ;yp=ystart
   ypb=ystart
   ODEINT, ystart, eps, -h1, string1, ypb     ; ODE INTEGRATION SUBROUTINE (BACK)
   iji=where((abs(ypb[0,*]) lt myxrange[1]) and (abs(ypb[1,*]) lt myyrange[1]) and (abs(ypb[2,*]) lt myzrange[1]))
   ypb=ypb[*,iji]
   ; start and end points of backwards field lines   
   ndb=size(ypb)
   fsb=[ [ypb[0,0],ypb[0,0]],[ypb[1,0],ypb[1,0]],[ypb[2,0],ypb[2,0]]]
   feb=[ [ypb[0,ndb[2]-1],ypb[0,ndb[2]-1]],[ypb[1,ndb[2]-1],ypb[1,ndb[2]-1]],[ypb[2,ndb[2]-1],ypb[2,ndb[2]-1]]]
   epf2=[ypb[0,ndb[2]-1],ypb[1,ndb[2]-1],ypb[2,ndb[2]-1]]
   oib=[ [-ypb[0,ndb[2]-1]+ypb[0,ndb[2]-2]],[-ypb[1,ndb[2]-1]+ypb[1,ndb[2]-2]],[-ypb[2,ndb[2]-1]+ypb[2,ndb[2]-2]]] 
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col1)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col1, thick=flinethick+4, /OVERPLOT
   IF keyword_set(arrows) THEN xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col1, SYMBOL=oARb, THICK=2, /OVERPLOT   ; the initial field line arrow
   ystart=[0,1e6,49e6]
   yp=ystart
   ;ypb=ystart
   ODEINT, ystart, eps, h1, string1, yp     ; ODE INTEGRATION SUBROUTINE (FORWARD)
   iji=where((abs(yp[0,*]) lt myxrange[1]) and (abs(yp[1,*]) lt myyrange[1]) and (abs(yp[2,*]) lt myzrange[1]))
   yp=yp[*,iji]
   nd=size(yp)
   ; start and end points of forward field lines   
   fs=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
   fe=[ [yp[0,nd[2]-1],yp[0,nd[2]-1]],[yp[1,nd[2]-1],yp[1,nd[2]-1]],[yp[2,nd[2]-1],yp[2,nd[2]-1]]]
   epb1=[yp[0,nd[2]-1],yp[1,nd[2]-1],yp[2,nd[2]-1]]
   oi=[ [yp[0,nd[2]-1]-yp[0,nd[2]-2]],[yp[1,nd[2]-1]-yp[1,nd[2]-2]],[yp[2,nd[2]-1]-yp[2,nd[2]-2]]] 
   ;orientation of forward arrow
   oModela=obj_new('IDLgrModel')
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col2)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col2, thick=flinethick+4,  /OVERPLOT
   ;xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col2, thick=4, /OVERPLOT
   IF keyword_set(arrows) THEN xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col2, SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
   ystart=[0,-1e6,49e6]
   yp=ystart
   ODEINT, ystart, eps, h1, string1, yp     ; ODE INTEGRATION SUBROUTINE (FORWARD)
   iji=where((abs(yp[0,*]) lt myxrange[1]) and (abs(yp[1,*]) lt myyrange[1]) and (abs(yp[2,*]) lt myzrange[1]))
   yp=yp[*,iji]
   nd=size(yp)
   ; start and end points of forward field lines   
   
   fs=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
   fe=[ [yp[0,nd[2]-1],yp[0,nd[2]-1]],[yp[1,nd[2]-1],yp[1,nd[2]-1]],[yp[2,nd[2]-1],yp[2,nd[2]-1]]]
   epb2=[yp[0,nd[2]-1],yp[1,nd[2]-1],yp[2,nd[2]-1]]
   oi=[ [yp[0,nd[2]-1]-yp[0,nd[2]-2]],[yp[1,nd[2]-1]-yp[1,nd[2]-2]],[yp[2,nd[2]-1]-yp[2,nd[2]-2]]]  
   ;orientation of forward arrow
   oModela=obj_new('IDLgrModel')
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col2)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col2, thick=flinethick+4,  /OVERPLOT
   IF keyword_set(arrows) THEN xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col2, SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
   
   rmax=720
   endpts1=dblarr(3,rmax)
   endpts2=dblarr(3,rmax) 
    
  FOR i=0,rmax-1 DO BEGIN 
   theta=DOUBLE(i)/DOUBLE(rmax+1)*2.*!dpi
   ystart=[rad*cos(theta)-0.0d0,0.0,1.0*rad*sin(theta)+5e7]
   ;yp=ystart
   ypb=ystart
   ODEINT, ystart, eps, -h1, string1, ypb     ; ODE INTEGRATION SUBROUTINE (BACK)
   iji=where((abs(ypb[0,*]) lt myxrange[1]) and (abs(ypb[1,*]) lt myyrange[1]) and (abs(ypb[2,*]) lt myzrange[1]))
   ypb=ypb[*,iji]
   ;; start and end points of backwards field lines   
   ndb=size(ypb)
   endpts1[0,i]=ypb[0,ndb[2]-1]
   endpts1[1,i]=ypb[1,ndb[2]-1]
   endpts1[2,i]=ypb[2,ndb[2]-1]
   ystart=[0.0,rad*cos(theta)-0.0d0,1.0*rad*sin(theta)-5e7]
   yp=ystart
   ODEINT, ystart, eps, h1, string1, yp     ; ODE INTEGRATION SUBROUTINE (FORWARD)
   iji=where((abs(yp[0,*]) lt myxrange[1]) and (abs(yp[1,*]) lt myyrange[1]) and (abs(yp[2,*]) lt myzrange[1]))
   yp=yp[*,iji]
   nd=size(yp)
   endpts2[0,i]=yp[0,nd[2]-1]
   endpts2[1,i]=yp[1,nd[2]-1]
   endpts2[2,i]=yp[2,nd[2]-1]
  ENDFOR
   iii= where((endpts1[0,*] gt myxrange[1]/lscl-0.001) and (endpts1[2,*] lt 61.0*zscl))
   nee=n_elements(iii)
   ufp1=dblarr(3,nee+1)
   ufp1[0,1:nee]= myxrange[1]
   ufp1[1,1:nee]=(endpts1[1,iii])[sort(endpts1[2,iii])]
   ufp1[2,1:nee]=(endpts1[2,iii])[sort(endpts1[2,iii])]
   ufp1[0,0]=myxrange[1]
   ufp1[1,0]=epf1[1]
   ufp1[2,0]=epf1[2]
   xplot3d, ufp1[0,*]/lscl, ufp1[1,*]/lscl, ufp1[2,*]/lscl, COLOR=col2, /OVERPLOT, thick=flinethick
   iii= where((endpts1[0,*] lt myxrange[0]/lscl+0.001) and (endpts1[2,*] lt 61*zscl))
   nee=n_elements(iii)
   ufp2=dblarr(3,nee+1)
   ;lfp2(0,0:nee-1)=(endpts1[0,iii])[sort(endpts1[2,iii])]
   ufp2[0,1:nee]=myxrange[0]
   ufp2[1,1:nee]=(endpts1[1,iii])[sort(endpts1[2,iii])]
   ufp2[2,1:nee]=(endpts1[2,iii])[sort(endpts1[2,iii])]
   ufp2[0,0]=myxrange[0]
   ufp2[1,0]=epf2[1]
   ufp2[2,0]=epf2[2]
   xplot3d, ufp2[0,*]/lscl, ufp2[1,*]/lscl, ufp2[2,*]/zscl, COLOR=col2, /OVERPLOT, thick=flinethick
   
   iii= where((endpts2[1,*] gt myyrange[1]/lscl-0.001) and (endpts2[2,*] gt -61.0*zscl))
   nee=n_elements(iii)
   lfp1=dblarr(3,nee+1)
   lfp1[0,0:nee-1]=(endpts2[0,iii])[sort(endpts2[2,iii])]
   lfp1[1,0:nee-1]= myyrange[1]
   lfp1[2,0:nee-1]=(endpts2[2,iii])[sort(endpts2[2,iii])]
   lfp1[0,nee]=epb1[0]
   lfp1[1,nee]=myyrange[1]
   lfp1[2,nee]=epb1[2]
   xplot3d, lfp1[0,*]/lscl, lfp1[1,*]/lscl, lfp1[2,*]/lscl, COLOR=col1, /OVERPLOT, thick=flinethick
   iii= where((endpts2[1,*] lt myyrange[0]/lscl+0.001) and (endpts2[2,*] gt -61*zscl))
   nee=n_elements(iii)
   lfp2=dblarr(3,nee+1)
   ;lfp2(0,0:nee-1)=
   lfp2[0,0:nee-1]=(endpts2[0,iii])[sort(endpts2[2,iii])]
   lfp2[1,0:nee-1]=myyrange[0]
   lfp2[2,0:nee-1]=(endpts2[2,iii])[sort(endpts2[2,iii])]
   lfp2[0,nee]=epb2[0]
   lfp2[1,nee]=myyrange[0]
   lfp2[2,nee]=epb2[2]
   xplot3d, lfp2[0,*]/lscl, lfp2[1,*]/lscl, lfp2[2,*]/zscl, COLOR=col1, /OVERPLOT, thick=flinethick
 
 ;  STOP
   
   ;each output closes the current window, therefore we need to do output and repeat if filename is set
   IF ((ij eq 0) and keyword_set(filename) and keyword_set(transparent)) THEN xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=flc, NAME='atestb', SYMBOL=oARb, THICK=2, /OVERPLOT, filename=newfilenames[0]    
   IF ((ij eq 0) and keyword_set(filename) and keyword_set(transparent)) THEN continue
  
  IF (keyword_set(isocurrent)) THEN BEGIN
   r1j=dblarr(3)
   r2j=dblarr(3)
   rstepsj=intarr(3)
   r1j[0]=-4.0e7	;xstart
   r2j[0]=4.0e7 	;xend
   rstepsj[0]=101 	;nx
   r1j[1]=-4.0e7	;ystart
   r2j[1]=4.0e7 	;yend
   rstepsj[1]=101	;ny
   r1j[2]=-8e7  	;zstart
   r2j[2]=8e7   	;zend
   rstepsj[2]=101 	;nz
   IF (rstepsj[0] EQ 1) THEN  gdxj=1.0d0 ELSE gdxj=1.0d0/FLOAT(rstepsj[0]-1)
   IF (rstepsj[1] EQ 1) THEN  gdyj=1.0d0 ELSE gdyj=1.0d0/FLOAT(rstepsj[1]-1)
   IF (rstepsj[2] EQ 1) THEN  gdzj=1.0d0 ELSE gdzj=1.0d0/FLOAT(rstepsj[2]-1)
   gdsj=[gdxj,gdyj,gdzj]
   lboxj=[r2j(0)-r1j(0),r2j(1)-r1j(1),r2j(2)-r1j(2)]
   npj=rstepsj[0]*rstepsj[1]*rstepsj[2]
   flinecol=[0,255,0]
   j=dblarr(rstepsj[0],rstepsj[1],rstepsj[2],3)
   modj=dblarr(rstepsj[0],rstepsj[1],rstepsj[2])
   FOR ix=0,rstepsj[0]-1 DO BEGIN
    ;print, ix
    FOR iy=0,rstepsj[1]-1 DO BEGIN
     FOR iz=0,rstepsj[2]-1 DO BEGIN
      i=[ix,iy,iz]
      ystart= R1j+lboxj*(i*1.0D0)*[gdxj,gdyj,gdzj]
      ;print, "ystart=", ystart
      j(ix,iy,iz,*)=current(ystart)
     ENDFOR
    ENDFOR
   ENDFOR
   modj=sqrt(j(*,*,*,0)*j(*,*,*,0)+j(*,*,*,1)*j(*,*,*,1)+j(*,*,*,2)*j(*,*,*,2))
   ;STOP
   ;IF NOT(keyword_set(isoccol)) THEN isoccol=[0,240,0]
   IF (tt ne 0) THEN BEGIN ;PLOTTING ONLY WORKS WHEN B IS NON-POTENTIAL
    lsize=1
    theModel = Obj_New('IDLgrModel')
    theModel -> Scale, ((r2j[0]-r1j[0])/rstepsj[0]/1e6),(r2j[1]-r1j[1])/rstepsj[1]/1e6,(r2j[2]-r1j[2])/rstepsj[2]/1e6
    obj_temp=mk_iso(modj(*,*,*),isocurrent,scol=isoccol, /low);, transparency=0.7)
    print, max(modj)
    
    IF (obj_valid(obj_temp)) THEN theModel -> Add, obj_temp
    newsym=obj_new('IDLgrSymbol', theModel)
    IF obj_valid(newsym) THEN newsym->setproperty, size=lsize
    IF obj_valid(newsym) THEN newsym->SetProperty, COLOR = isocol
    IF keyword_set(filename) THEN BEGIN
     IF keyword_set(transparent) THEN BEGIN
      xplot3d, [r1j[0],r1j[0]]/1e6, [r1j[1],r1j[1]]/1e6, [r1j[2],r1j[2]]/1e6, symbol=newsym, /overplot, filename=newfilenames[1] 
     ENDIF ELSE BEGIN
      xplot3d, [r1j[0],r1j[0]]/1e6, [r1j[1],r1j[1]]/1e6, [r1j[2],r1j[2]]/1e6, symbol=newsym, /overplot, filename=filename
     ENDELSE
    ENDIF ELSE BEGIN
     xplot3d, [r1j[0],r1j[0]]/1e6, [r1j[1],r1j[1]]/1e6, [r1j[2],r1j[2]]/1e6, symbol=newsym, /overplot
    ENDELSE 
   ENDIF
  ENDIF
ENDFOR  

;if transparency then combine both temp images to be transparent
IF (keyword_set(filename) AND keyword_set(transparent)) THEN BEGIN
 im1=read_image(newfilenames[0])
 im2=read_image(newfilenames[1])
 window, 0, xsize=filewidth, ysize=filewidth
 tv, im1, true=1
 tv, byte((1-transparent)*im1+transparent*im2), true=1
 write_image,filename,type[1],byte((1-transparent)*im1+transparent*im2)
 wdelete, 0
ENDIF
ENDFOR
END
;------------------------------------;
;---OLD VERSION, FIELD WRAPPING------;
;------------------------------------;
PRO fieldwrapping, op=op, myt=myt, mya=mya, rad=rad, midrad=midrad, zloc=zloc, intercept=intercept, isocurrent=isocurrent, $
isoccol=isoccol, ax=ax, az=az, filename=filename, transparent=transparent, filewidth=filewidth, frac=frac, centre=centre
;@ODEINT
;@xplot3dJT
;USES HACKED VERSION OF xplot3d, called xplot3dJT

common somename, tt, FRon, len, an

IF n_elements(op) eq 0     	THEN op=0
IF n_elements(myt) eq 0     	THEN myt=0.0d0 
IF n_elements(mya) eq 0     	THEN mya=0.1d0 
IF n_elements(rad) eq 0     	THEN rad=5000000.0d0
IF n_elements(midrad) eq 0   	THEN midrad=3000000.0d0
IF n_elements(zloc) eq 0    	THEN zloc=5e7
IF n_elements(intercept) eq 0 	THEN intercept=5e7
IF n_elements(isoccol) eq 0 	THEN isoccol=[200,0,0]
IF n_elements(ax) eq 0      	THEN ax=-60
IF n_elements(az) eq 0      	THEN az=30
IF n_elements(filewidth) eq 0   THEN filewidth=800
IF n_elements(filename) eq 0  	THEN filename='fieldwrappingtest.png'
IF n_elements(frac) eq 0   	THEN frac=100.0
IF n_elements(centre) eq 0  	THEN centre=0

;type=strsplit(filename,".",/extract)
IF keyword_set(filename) and keyword_set(transparent) THEN BEGIN
 newfilenames=[type[0]+'_a.'+type[1],type[0]+'_b.'+type[1]]
ENDIF

 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]
 purple=[153,51,255]
 grey=[170,170,170]
 poobrown=[153,76,0]
 olive=[128,128,0]
 salmon=[255,160,122]
 lightsteelblue=[176,196,222]
 lavender=[230,230,250]
 rosybrown=[188,143,143]
 burlywood=[222,184,135]
 tann=[210,180,140]
 
col1=lightsteelblue
col2=tann


an=mya
me=0
km=0
Mm=1

myxrange=[-100e6,100e6]
myyrange=[-100e6,100e6]
myzrange=[-75e6,75e6]
lscl=1e6
zscl=1e6 
mytitle=' '

h1=0.01     	    	    ; step size?
string1='derivs' 	    ; a string
FRon=1	            ; flux ring ON?
tt=myt     	    ; time
eps=0.0000001    	    ; error
zt='z (Mm)'

IF keyword_set(me) THEN BEGIN
  lscl=1.0
  rmult=1.0
  xt='x (m)'
  yt='y (m)'
ENDIF  
IF keyword_set(km) THEN BEGIN
  lscl=1e3
  rmult=1.0
  xt='x (km)'
  yt='y (km)'
ENDIF
IF keyword_set(Mm) THEN BEGIN
  lscl=1e6
  rmult=1.0
  xt='x (Mm)'
  yt='y (Mm)'
ENDIF

 xlen=DOUBLE(myxrange[1]-myxrange[0])
 ylen=DOUBLE(myyrange[1]-myyrange[0])
 zlen=DOUBLE(myzrange[1]-myzrange[0])

 thisPalette = Obj_New('IDLgrPalette')
 thisPalette->LoadCT, 34
 dims=[1.,1.,zlen/xlen]
 mysymsize=[lscl,lscl,zscl]*dims*1.5e-7*frac  
 
 oOrb4 = OBJ_NEW('orb', COLOR=[034, 139 ,034]) 	;Forest Green orb
 oOrb5 = OBJ_NEW('orb', COLOR=[124,252,0]) 	;Other green orb
 ;oOrb4 = OBJ_NEW('orb', COLOR=[255, 255 ,255]) 	;white orb
 ;oOrb5 = OBJ_NEW('orb', COLOR=[0, 0 ,0]) 	;black orb
 oOrb4->Scale, frac, frac, zlen/xlen*frac
 oOrb5->Scale, frac, frac, zlen/xlen*frac
 oSy4 = OBJ_NEW('IDLgrSymbol', oOrb4) ;osymbol4 is black orb for null point 
 oSy5 = OBJ_NEW('IDLgrSymbol', oOrb5) ;osymbol4 is black orb for null point 
 
 nullu=[ [0.0d0,0.0d0],[0.0d0,0.0d0],[5e7,5e7]]
 nulld=[ [0.0d0,0.0d0],[0.0d0,0.0d0],[-5e7,-5e7]]   
 
 oModel=obj_new('IDLgrModel')
 oModel->add,symbol_obj(2,col1,[0,0,0])	    ; steelblue tetrahedon
 oTet = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oTet) THEN oTet->SetProperty, Size=[frac, frac, zlen/xlen*frac]
 oModel=obj_new('IDLgrModel')
 oModel->add,symbol_obj(2,col2,[0,0,0])	    ; grey tetrahedon
 oTet2 = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oTet2) THEN oTet2->SetProperty, Size=[frac, frac, zlen/xlen*frac]
 
print, 'making a movie..'
naz=00


FOR ii=0,naz DO BEGIN
 az=30+ii
 ;filename=type[0]+string(ii,format='(i3.3,".")')+type[1]
 print, ii, naz, format='("frame ", i3.3, "/", i3.3)' 
IF (keyword_set(filename) and keyword_set(transparent)) THEN ns=1 ELSE ns=0
FOR ij=0,ns DO BEGIN  ; for transparency (as outputting closes the xWindow)
 tmax=360
 tinc=24
 ;flc=[127,127,127]  	    ; field line and arrow color (grey)
 IF (op EQ 0) THEN BEGIN
  xplot3d, nullu[*,0]/lscl, nullu[*,1]/lscl, nullu[*,2]/zscl, COLOR=flc, xrange=myxrange/lscl, yrange=myyrange/lscl, zrange=myzrange/lscl, $
  xtitle=xt, ytitle=yt, ztitle=zt, title=mytitle, linestyle=ls, SYMBOL=oTet,ax=ax,az=az, filewidth=filewidth
 ENDIF ELSE BEGIN
  xplot3d, nullu[*,0]/lscl, nullu[*,1]/lscl, nullu[*,2]/zscl, COLOR=flc, linestyle=ls, SYMBOL=oTet, /OVERPLOT
 ENDELSE
 
 xplot3d, nulld[*,0]/lscl, nulld[*,1]/lscl, nulld[*,2]/zscl, COLOR=flc, SYMBOL=oTet2, /OVERPLOT
;multiple spine rings 
  FOR i=0,tmax,tinc DO BEGIN 
   theta=DOUBLE(i)
   ystart=[intercept,rad*cos(theta)-0.0d0, rad*sin(theta)-zloc]
   yp=ystart
   ypb=ystart
   ;flc=[theta/DOUBLE(tmax)*255,0,0]  	    ; field line and arrow color
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
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col2)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col2)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col2, /OVERPLOT
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col2, /OVERPLOT
   xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col2, SYMBOL=oAR, THICK=2, /OVERPLOT  
   xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col2, SYMBOL=oARb, THICK=2, /OVERPLOT   
  ENDFOR
  ;STOP
  FOR i=0,tmax,tinc DO BEGIN 
   theta=DOUBLE(i)
   ystart=[-intercept,rad*cos(theta)-0.0d0, rad*sin(theta)-zloc]
   yp=ystart
   ypb=ystart
   ;flc=[theta/DOUBLE(tmax)*255,0,0]  	    ; field line and arrow color
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
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col2)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col2)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col2, /OVERPLOT
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col2, /OVERPLOT
   xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col2, SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
   xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col2, SYMBOL=oARb, THICK=2, /OVERPLOT   ; the initial field line arrow
  ENDFOR
  FOR i=0,tmax,tinc DO BEGIN 
   theta=DOUBLE(i)
   ystart=[rad*cos(theta)-0.0d0,intercept,rad*sin(theta)+zloc]
   yp=ystart
   ypb=ystart
   ;flc=[theta/DOUBLE(tmax)*255,0,0]  	    ; field line and arrow color
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
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col1)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col1)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col1, /OVERPLOT
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col1, /OVERPLOT
   xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col1, NAME='atest', SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
   xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col1, NAME='atestb', SYMBOL=oARb, THICK=2, /OVERPLOT   ; the initial field line arrow
  ENDFOR
  FOR i=0,tmax,tinc DO BEGIN 
   theta=DOUBLE(i)
   ystart=[rad*cos(theta)-0.0d0,-intercept,rad*sin(theta)+zloc]
   yp=ystart
   ypb=ystart
   ;flc=[theta/DOUBLE(tmax)*255,0,0]  	    ; field line and arrow color
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
   oModela->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=col1)
   oAR = obj_new('IDLgrSymbol', oModela)
   IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
   ;orientation of reverse arrow
   oModelb=obj_new('IDLgrModel')
   oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=col1)
   oARb = obj_new('IDLgrSymbol', oModelb)
   IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
   xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=col1, /OVERPLOT
   xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=col1, /OVERPLOT
   xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=col1, SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
   xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=col1, SYMBOL=oARb, THICK=2, /OVERPLOT   ; the initial field line arrow
  ENDFOR
  IF centre THEN BEGIN
   flc=[0,0,205] ; change field line color to blue
   tinc=45
   FOR i=0,tmax,tinc DO BEGIN 
    theta=DOUBLE(i)+1.0d0
    ystart=[midrad*cos(theta), midrad*sin(theta),0.0d0]
    yp=ystart
    ypb=ystart
    ;flc=[theta/DOUBLE(tmax)*255,0,0]  	    ; field line and arrow color
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
    IF Obj_Valid(oAR) THEN oAR->SetProperty, Size=mysymsize
    ;orientation of reverse arrow
    oModelb=obj_new('IDLgrModel')
    oModelb->add,mk_vector([oib[0],oib[1],oib[2]]/sqrt(oib[0]*oib[0]+oib[1]*oib[1]+oib[2]*oib[2])/dims,color=flc)
     oARb = obj_new('IDLgrSymbol', oModelb)
    IF Obj_Valid(oARb) THEN oARb->SetProperty, Size=mysymsize
    xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=flc, thick=2, /OVERPLOT
    xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=flc, thick=2, /OVERPLOT
    xplot3d, fe[*,0]/lscl, fe[*,1]/lscl, fe[*,2]/zscl, COLOR=flc, NAME='atest', SYMBOL=oAR, THICK=2, /OVERPLOT   	; the final field line arrow
    xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=flc, NAME='atestb', SYMBOL=oARb, THICK=2, /OVERPLOT    ; the initial field line arrow
   ENDFOR
  ENDIF
   ;each output closes the current window, therefore we need to do output and repeat if filename is set
   IF ((ij eq 0) and keyword_set(filename) and keyword_set(transparent)) THEN xplot3d, feb[*,0]/lscl, feb[*,1]/lscl, feb[*,2]/zscl, COLOR=flc, NAME='atestb', SYMBOL=oARb, THICK=2, /OVERPLOT, filename=newfilenames[0]    
   IF ((ij eq 0) and keyword_set(filename) and keyword_set(transparent)) THEN continue
  
  IF (keyword_set(isocurrent)) THEN BEGIN
   r1j=dblarr(3)
   r2j=dblarr(3)
   rstepsj=intarr(3)
   r1j[0]=-4.0e7	;xstart
   r2j[0]=4.0e7 	;xend
   rstepsj[0]=101 	;nx
   r1j[1]=-4.0e7	;ystart
   r2j[1]=4.0e7 	;yend
   rstepsj[1]=101	;ny
   r1j[2]=-8e7  	;zstart
   r2j[2]=8e7   	;zend
   rstepsj[2]=101 	;nz
   IF (rstepsj[0] EQ 1) THEN  gdxj=1.0d0 ELSE gdxj=1.0d0/FLOAT(rstepsj[0]-1)
   IF (rstepsj[1] EQ 1) THEN  gdyj=1.0d0 ELSE gdyj=1.0d0/FLOAT(rstepsj[1]-1)
   IF (rstepsj[2] EQ 1) THEN  gdzj=1.0d0 ELSE gdzj=1.0d0/FLOAT(rstepsj[2]-1)
   gdsj=[gdxj,gdyj,gdzj]
   lboxj=[r2j(0)-r1j(0),r2j(1)-r1j(1),r2j(2)-r1j(2)]
   npj=rstepsj[0]*rstepsj[1]*rstepsj[2]
   flinecol=[0,255,0]
   j=dblarr(rstepsj[0],rstepsj[1],rstepsj[2],3)
   modj=dblarr(rstepsj[0],rstepsj[1],rstepsj[2])
   FOR ix=0,rstepsj[0]-1 DO BEGIN
    ;print, ix
    FOR iy=0,rstepsj[1]-1 DO BEGIN
     FOR iz=0,rstepsj[2]-1 DO BEGIN
      i=[ix,iy,iz]
      ystart= R1j+lboxj*(i*1.0D0)*[gdxj,gdyj,gdzj]
      ;print, "ystart=", ystart
      j(ix,iy,iz,*)=current(ystart)
     ENDFOR
    ENDFOR
   ENDFOR
   modj=sqrt(j(*,*,*,0)*j(*,*,*,0)+j(*,*,*,1)*j(*,*,*,1)+j(*,*,*,2)*j(*,*,*,2))
   ;STOP
   ;IF NOT(keyword_set(isoccol)) THEN isoccol=[0,240,0]
   IF (tt ne 0) THEN BEGIN ;PLOTTING ONLY WORKS WHEN B IS NON-POTENTIAL
    lsize=1
    theModel = Obj_New('IDLgrModel')
    theModel -> Scale, ((r2j[0]-r1j[0])/rstepsj[0]/1e6),(r2j[1]-r1j[1])/rstepsj[1]/1e6,(r2j[2]-r1j[2])/rstepsj[2]/1e6
    obj_temp=mk_iso(modj(*,*,*),isocurrent,scol=isoccol, /low);, transparency=0.7)
    print, max(modj)
    
    IF (obj_valid(obj_temp)) THEN theModel -> Add, obj_temp
    newsym=obj_new('IDLgrSymbol', theModel)
    IF obj_valid(newsym) THEN newsym->setproperty, size=lsize
    IF obj_valid(newsym) THEN newsym->SetProperty, COLOR = isocol
    IF keyword_set(filename) THEN BEGIN
     IF keyword_set(transparent) THEN BEGIN
      xplot3d, [r1j[0],r1j[0]]/1e6, [r1j[1],r1j[1]]/1e6, [r1j[2],r1j[2]]/1e6, symbol=newsym, /overplot, filename=newfilenames[1] 
     ENDIF ELSE BEGIN
      xplot3d, [r1j[0],r1j[0]]/1e6, [r1j[1],r1j[1]]/1e6, [r1j[2],r1j[2]]/1e6, symbol=newsym, /overplot, filename=filename
     ENDELSE
    ENDIF ELSE BEGIN
     xplot3d, [r1j[0],r1j[0]]/1e6, [r1j[1],r1j[1]]/1e6, [r1j[2],r1j[2]]/1e6, symbol=newsym, /overplot
    ENDELSE 
   ENDIF
  ENDIF
ENDFOR  

;if transparency then combine both temp images to be transparent
IF (keyword_set(filename) AND keyword_set(transparent)) THEN BEGIN
 im1=read_image(newfilenames[0])
 im2=read_image(newfilenames[1])
 window, 0, xsize=filewidth, ysize=filewidth
 tv, im1, true=1
 tv, byte((1-transparent)*im1+transparent*im2), true=1
 write_image,filename,type[1],byte((1-transparent)*im1+transparent*im2)
 wdelete, 0
ENDIF
ENDFOR

  
END
