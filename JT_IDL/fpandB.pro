;plots final position (coloured by peak energy gain) against field lines in global Hesse model scaled to be AR scale
@xplot3dJT
@ODEINT
common somename, epsilon

 npts=32*32*8
 dims=[2.,2.,4.]
 xlen=dims[0]
 ylen=dims[1]
 zlen=dims[2]
 mysymsize=[1.0,1.0,1.0]*0.08*dims  

tempdir='../epsilon'
newoutdir='./allepsilon/finalpos/'

;;----- define variables ----------------
eps=0.0000001
h1=0.0001     	    ; step size?
string1='derivs2'

ery=[0,10,19]

FOR e=0, 2 DO BEGIN


es=string(ery(e),format='(i2.2)')
print, 'sorting epsilon '+es


result = FILE_TEST(tempdir+es+'/kestore.sav') 
IF (result eq 0) THEN BEGIN
;stop
; npts=16*16*2
 nom=tempdir+es+'/Data/'
 print, 'loading: ', nom
 ds=getrdata(1,/ke, /gyro, wkdir=nom+'DataR/')
 inike=ds.ke(0)
 ;loadct, 2
 estore=0
 kestore=dblarr(npts)
 maxgyro=dblarr(npts)
 tstop=dblarr(npts)
 startpt=dblarr(3,npts)
 endpt=dblarr(3,npts)
 ns=n_elements(ds.t)
 kestore[0]=inike
 maxg=max(ds.gyror)
 maxgyro[0]=maxg
 startpt[*,0]=[ds.x[0],ds.y[0],ds.z[0]]
 endpt[*,0]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
 tstop[0]=ds.t[ns-1]
 for i=1,npts DO BEGIN
  tflag=0
  print, i, npts, format='(i4,"/",i4)'
  ds=getrdata(i,/ke,/gyro, wkdir=nom+'DataR/')
  ns=n_elements(ds.t)
  startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
  endpt[*,i-1]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
  tstop[i-1]=ds.t[ns-1]
  maxke=max(ds.ke)
  kestore[i-1]=maxke
  ;print, startpt[*,i-1], 'max ke=', maxke
  maxg=max(ds.gyror)
  maxgyro[i-1]=maxg
 endfor
 save, npts, kestore, startpt, endpt, tstop, maxgyro, filename=tempdir+es+'/kestore.sav'
ENDIF ELSE BEGIN
;STOP
restore, filename=tempdir+es+'/kestore.sav', /verbose
ENDELSE


;;----- pick initial conditions ------------
ystart=dblarr(3)
r1=dblarr(3)
r2=dblarr(3)
rsteps=intarr(3)

r1[0]=-4.0	;xstart
r2[0]=4.0 	;xend
rsteps[0]=5 	;nx
r1[1]=0.0	;ystart
r2[1]=0.0 	;yend
rsteps[1]=1	;ny
r1[2]=4.0  	;zstart
r2[2]=15.0   	;zend
rsteps[2]=11 	;nz

lscl=1.0
zscl=1.0
flinecol=[170,170,170]

flcol=[183,183,183]
IF (rsteps[0] EQ 1) THEN  gdx=1.0d0 ELSE gdx=1.0d0/FLOAT(rsteps[0]-1)
IF (rsteps[1] EQ 1) THEN  gdy=1.0d0 ELSE gdy=1.0d0/FLOAT(rsteps[1]-1)
IF (rsteps[2] EQ 1) THEN  gdz=1.0d0 ELSE gdz=1.0d0/FLOAT(rsteps[2]-1)
gds=[gdx,gdy,gdz]
lbox=[r2(0)-r1(0),r2(1)-r1(1),r2(2)-r1(2)]

np=rsteps[0]*rsteps[1]*rsteps[2]
flinecol=[170,170,170]
purple=[148,0,211]
pink=[255,0,255]
flcol=[183,183,183]

floc='mov/'
esteps=30
estart=0
eend=15
;FOR ie=0,esteps-1 DO BEGIN
 pno=0
; epsilon=DOUBLE(eend-estart)*ie/double(esteps)
 epsilon=e
; FOR ix=0,rsteps[0]-1 DO BEGIN
;  FOR iy=0,rsteps[1]-1 DO BEGIN
;   FOR iz=0,rsteps[2]-1 DO BEGIN
;    i=[ix,iy,iz]
;    ystart= R1+lbox*(i*1.0D0)*[gdx,gdy,gdz]
;    yp=ystart
;    ODEINT, ystart, eps, h1, string1, yp
;    n_data=size(yp)
;    ts=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
;    te=[ [yp[0,n_data[2]-1],yp[0,n_data[2]-1]],[yp[1,n_data[2]-1],yp[1,n_data[2]-1]],[yp[2,n_data[2]-1],yp[2,n_data[2]-1]]]
;    oi=[ [yp[0,n_data[2]-1]-yp[0,n_data[2]-2]],[yp[1,n_data[2]-1]-yp[1,n_data[2]-2]],[yp[2,n_data[2]-1]-yp[2,n_data[2]-2]]]  
;    oModel=obj_new('IDLgrModel')
;    oModel->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=flinecol)
;    oAR3 = obj_new('IDLgrSymbol', oModel)
;    IF Obj_Valid(oAR3) THEN oAR3->SetProperty, Size=mysymsize
;    IF pno eq 0 THEN BEGIN
;     xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, xrange=[-10,10], yrange=[-10,10], zrange=[0,20], COLOR=flinecol, $
;     ztitle='z', xtitle='x',ytitle='y', title=string(epsilon, format='("epsilon=",f4.1)'), filewidth=800
;    ENDIF ELSE BEGIN
;     xplot3d, yp(0,*)/lscl, yp(1,*)/lscl, yp(2,*)/zscl, COLOR=[0,0,255], /OVERPLOT
;    ENDELSE
;    ;xplot3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/zscl, COLOR=[0,0,0], NAME='atest', SYMBOL=oAR3, THICK=2, /OVERPLOT
;    pno=pno+1
;    ypb=ystart
;    ODEINT, ystart, eps, -h1, string1, ypb 
;    ;IF ((ix eq rsteps[0]-1) AND (iy eq rsteps[1]-1) AND (iz eq rsteps[2]-1)) THEN BEGIN
;    ; xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=[0,0,255], /OVERPLOT, filename=floc+string(ie,format='("flux",i2.2,".png")')
;    ;ENDIF ELSE BEGIN
;     xplot3d, ypb(0,*)/lscl, ypb(1,*)/lscl, ypb(2,*)/zscl, COLOR=[0,0,255], /OVERPLOT
;    ;ENDELSE
;   ENDFOR
;  ENDFOR
; ENDFOR

arbscl=10. 
lscl=1e6
zscl=1e6

frac=0.12/arbscl

oOrb = OBJ_NEW('orb', COLOR=[170,170,170]) 	    	    	    	    	;initial orb - grey
oOrb->Scale, frac, frac, zlen/xlen*frac
oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol for start point
 
loadct, 4
;stretch, 30,270
tvlct, RR, GG, BB, /GET
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

t1=3
t2=6.3
alke=alog10(kestore-20.0d0)
astore=fix((alke-t1)/(t2-t1)*254)
astore(where(astore lt 0))=0

 ystart=startpt[*,0]/lscl*arbscl
 yp=ystart
 ODEINT, ystart, eps, h1, string1, yp
 extremely1=extrema(reform(yp[2,*])) ;returns the location(s) of extrema in the z component of yp
 n_data=size(yp)
 ts=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
 te=[ [yp[0,n_data[2]-1],yp[0,n_data[2]-1]],[yp[1,n_data[2]-1],yp[1,n_data[2]-1]],[yp[2,n_data[2]-1],yp[2,n_data[2]-1]]]
 oi=[ [yp[0,n_data[2]-1]-yp[0,n_data[2]-2]],[yp[1,n_data[2]-1]-yp[1,n_data[2]-2]],[yp[2,n_data[2]-1]-yp[2,n_data[2]-2]]]
 
 xplot3d, [0,0], [0,0], [0,0], xrange=[-1,1], yrange=[-1,1], zrange=[0,4], COLOR=[255,255,255],$
     ztitle='z (Mm [L/100km])', xtitle='x (Mm [L/100km])',ytitle='y (Mm [L/100km])',$; title=string(epsilon, format='("epsilon=",f4.1)'), 
     filewidth=900, az=40, AX=-65

 ypb=ystart
 ODEINT, ystart, eps, -h1, string1, ypb 
 
 extremely2=extrema(reform(ypb[2,*])) ;returns the location(s) of extrema in the z component of yp
 totextremely=[extremely1,extremely2]
 IF (n_elements(totextremely) gt 2) THEN BEGIN
  nflcol=purple
 ENDIF ELSE BEGIN
  IF (((max(yp[2,*])/arbscl) gt 0.75) or (max(ypb[2,*])/arbscl gt 0.75)) THEN BEGIN
   nflcol=flcol
  ENDIF ELSE BEGIN
   nflcol=[0,0,0]
  ENDELSE
 ENDELSE
 
 oModel=obj_new('IDLgrModel')
 oModel->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=nflcol)
 oAR3 = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oAR3) THEN oAR3->SetProperty, Size=mysymsize
 xplot3d, te[*,0]/arbscl,te[*,1]/arbscl, te[*,2]/arbscl, COLOR=nflcol, sym=oAR3, /OVERPLOT  
 xplot3d, yp(0,*)/arbscl, yp(1,*)/arbscl, yp(2,*)/arbscl, COLOR=nflcol, /OVERPLOT  
 xplot3d, ypb(0,*)/arbscl, ypb(1,*)/arbscl, ypb(2,*)/arbscl, COLOR=nflcol, /OVERPLOT
 
 oOrb2 = OBJ_NEW('orb', COLOR=[rr(astore(0)),gg(astore(0)),bb(astore(0))]) 	;final orb - coloured with ke
 oOrb2->Scale, frac, frac, zlen/xlen*frac
 oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 for stop point
 ts=[ [startpt[0,0],startpt[0,0]],[startpt[1,0],startpt[1,0]],[startpt[2,0],startpt[2,0]]]
 te=[ [endpt[0,0],endpt[0,0]],[endpt[1,0],endpt[1,0]],[endpt[2,0],endpt[2,0]]]
 XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/lscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
 
 XPLOT3D, [-1,1], [0,0], [0,0], COLOR=[0,0,0], NAME='end', THICK=3, /OVERPLOT, linestyle=2
 
 fieldflag=0
FOR i=2,npts DO BEGIN
 IF ((i MOD 133) eq 0) THEN fieldflag=1 ELSE fieldflag=0
 oOrb2 = OBJ_NEW('orb', COLOR=[rr(astore(i-1)),gg(astore(i-1)),bb(astore(i-1))]) 	;final orb - coloured with ke
 oOrb2->Scale, frac, frac, zlen/xlen*frac
 oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 for stop point
 ts=[ [startpt[0,i-1],startpt[0,i-1]],[startpt[1,i-1],startpt[1,i-1]],[startpt[2,i-1],startpt[2,i-1]]]
 te=[ [endpt[0,i-1],endpt[0,i-1]],[endpt[1,i-1],endpt[1,i-1]],[endpt[2,i-1],endpt[2,i-1]]]

 XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/lscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
 ystart=startpt[*,i-1]/lscl*arbscl
 yp=ystart
 ODEINT, ystart, eps, h1, string1, yp
 
 extremely1=extrema(reform(yp[2,*])) ;returns the location(s) of extrema in the z component of yp
 ;STOP
 
 n_data=size(yp)
 ts=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
 te=[ [yp[0,n_data[2]-1],yp[0,n_data[2]-1]],[yp[1,n_data[2]-1],yp[1,n_data[2]-1]],[yp[2,n_data[2]-1],yp[2,n_data[2]-1]]]
 oi=[ [yp[0,n_data[2]-1]-yp[0,n_data[2]-2]],[yp[1,n_data[2]-1]-yp[1,n_data[2]-2]],[yp[2,n_data[2]-1]-yp[2,n_data[2]-2]]]  
 oModel=obj_new('IDLgrModel')
 oModel->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=nflcol)
 oAR3 = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oAR3) THEN oAR3->SetProperty, Size=mysymsize
 
 
 pno=pno+1
 ypb=startpt[*,i-1]/lscl*arbscl
 ODEINT, ystart, eps, -h1, string1, ypb
 
 extremely2=extrema(reform(ypb[2,*])) ;returns the location(s) of extrema in the z component of yp
 totextremely=[extremely1,extremely2]
 IF (n_elements(totextremely) gt 2) THEN BEGIN	; if there are multiple maxima colour purple
  nflcol=purple
 ENDIF ELSE BEGIN
  IF (((max(yp[2,*])/arbscl) gt 0.75) or (max(ypb[2,*])/arbscl gt 0.75)) THEN BEGIN ; if higher loops, colour pink
   nflcol=flcol
  ENDIF ELSE BEGIN  	    	; if low loops, colour grey
   nflcol=[0,0,0]
  ENDELSE
 ENDELSE
 IF (fieldflag) THEN xplot3d, yp(0,*)/arbscl, yp(1,*)/arbscl, yp(2,*)/arbscl, COLOR=nflcol, /OVERPLOT
 IF (fieldflag) THEN xplot3d, ypb(0,*)/arbscl, ypb(1,*)/arbscl, ypb(2,*)/arbscl, COLOR=nflcol, /OVERPLOT
 oi=[ [yp[0,n_data[2]-1]-yp[0,n_data[2]-2]],[yp[1,n_data[2]-1]-yp[1,n_data[2]-2]],[yp[2,n_data[2]-1]-yp[2,n_data[2]-2]]]  
 oModel=obj_new('IDLgrModel')
 oModel->add,mk_vector([oi[0],oi[1],oi[2]]/sqrt(oi[0]*oi[0]+oi[1]*oi[1]+oi[2]*oi[2])/dims,color=nflcol)
 oAR3 = obj_new('IDLgrSymbol', oModel)
 IF Obj_Valid(oAR3) THEN oAR3->SetProperty, Size=mysymsize
 IF (fieldflag) THEN xplot3d, te[*,0]/arbscl,te[*,1]/arbscl, te[*,2]/arbscl, COLOR=nflcol, sym=oAR3, /OVERPLOT
 
ENDFOR
 XPLOT3D, [0,0], [0,0], [0,0], COLOR=[255,255,255], NAME='end', /OVERPLOT, filename=newoutdir+string(e,format='("fp",i2.2,".png")')
ENDFOR
;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
;title='Particle peak energy gain (eV)', /vertical, /ylog


end
