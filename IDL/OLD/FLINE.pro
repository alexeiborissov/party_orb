;PRO FLINE
; FLINE
;  numerically integrates field lines according to field specified in derivs
;  then displays
;  REQUIRES - ODEINT, derivs, RKQC
@ODEINT
;common path, kount;, yp
common somename, tt, FRon
 ;thisPalette = Obj_New('IDLgrPalette')
 ;thisPalette->LoadCT, 34
FRon=1
; oOrb = OBJ_NEW('orb', COLOR=[0, 255 ,0]) 
; oOrb2 = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 
; oOrb->Scale, .04, .04, .25 
; oOrb2->Scale, .04, .04, .25 
; oSymbol = OBJ_NEW('IDLgrSymbol', oOrb)     ;oSymbol is green orb for start point
; oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2)   ;osymbol2 is red orb for stop point

tt=0.3d

;;----- define variables ----------------
eps=0.0000001
h1=0.01
string1='derivs'
tim=0.3

;;----- pick initial conditions ------------
ystart=dblarr(3)
r1=dblarr(3)
r2=dblarr(3)
rsteps=intarr(3)


r1[0]=-0.5d0
r2[0]=0.5d0
rsteps[0]=7
r1[1]=-0.5d0
r2[1]=0.5d0
rsteps[1]=7
r1[2]=-5.5d0
r2[2]=5.5d0
rsteps[2]=15

IF (rsteps[0] EQ 1) THEN  gdx=1.0d0 ELSE gdx=1.0d0/FLOAT(rsteps[0]-1)
IF (rsteps[1] EQ 1) THEN  gdy=1.0d0 ELSE gdy=1.0d0/FLOAT(rsteps[1]-1)
IF (rsteps[2] EQ 1) THEN  gdz=1.0d0 ELSE gdz=1.0d0/FLOAT(rsteps[2]-1)
gds=[gdx,gdy,gdz]
lbox=[r2(0)-r1(0),r2(1)-r1(1),r2(2)-r1(2)]

pno=0
np=rsteps[0]*rsteps[1]*rsteps[2]

   ;x_start=0.15
   ;y_start=0.0
   ;z_start=0.0
   ;ystart=[, y_start, z_start]

FOR ix=0,rsteps[0]-1 DO BEGIN
 FOR iy=0,rsteps[1]-1 DO BEGIN
  FOR iz=0,rsteps[2]-1 DO BEGIN

   i=[ix,iy,iz]
   ystart= R1+lbox*(i*1.0D0)*gds
   yp=ystart
   ODEINT, ystart, eps, h1, string1, yp
   n_data=size(yp)
      
     ts=[ [yp[0,0],yp[0,0]],[yp[1,0],yp[1,0]],[yp[2,0],yp[2,0]]]
     te=[ [yp[0,n_data[2]-1],yp[0,n_data[2]-1]],[yp[1,n_data[2]-1],yp[1,n_data[2]-1]],[yp[2,n_data[2]-1],yp[2,n_data[2]-1]]]

   
    IF pno eq 0 THEN BEGIN
     xplot3d, yp(0,*), yp(1,*), yp(2,*), xrange=[-1,1], yrange=[-1,1], zrange=[-6,6], COLOR=[0,0,255]
    ENDIF ELSE BEGIN
     xplot3d, yp(0,*), yp(1,*), yp(2,*), COLOR=[0,0,255], /OVERPLOT
    ENDELSE
   ; XPLOT3D, ts[*,0], ts[*,1], ts[*,2], COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol, THICK=2, /OVERPLOT
   ; XPLOT3D, te[*,0], te[*,1], te[*,2], COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
    pno=pno+1
  ENDFOR
 ENDFOR
ENDFOR


;;----3D Plot of Field Lines----------------------------------------------

;xplot3d, yp(0,*), yp(1,*), yp(2,*), xrange=[-1,1], yrange=[-1,1], zrange=[-6,6] 


end
