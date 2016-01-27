PRO pvels, ii, wkdir=wkdir, silent=silent, mvals=mvals, lscl=lscl, tscl=tscl, bscl=bscl, stride=stride
;pvels takes a given RV file number, reads in various parameters and calculates some velocity/force vectors
; it then plots the particle track for that particle, and the vectors at various times, determined by stride.


;@ODEINT
;@mk_vector
;heap_gc
;common somename, tt, FRon, an

IF n_elements(wkdir) eq 0 THEN wkdir="../Data/"
IF n_elements(mvals) eq 0 THEN mvals=0
IF n_elements(lscl) eq 0 THEN lscl=1e7
IF n_elements(bscl) eq 0 THEN bscl=0.01
IF n_elements(tscl) eq 0 THEN tscl=100.0
IF n_elements(stride) eq 0 THEN stride=15

vscl=lscl/tscl
escl=vscl*bscl
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

dims=[1.,1.,6.]
mysymsize=[lscl,lscl,lscl]*0.1*dims  

files=FINDFILE(wkdir+'RV*.dat',count=count)
IF count eq 0 THEN BEGIN
 PRINT, "ERROR: no RV*.dat files found"
 RETURN
ENDIF
IF ii lt 1 THEN BEGIN
 print, "ERROR: data files from ii=1.."+STRING(COUNT,FORMAT='(i3.3)')
 RETURN
ENDIF

ds=getdata(ii,wkdir=wkdir,/drifts,/fields, /vpar)
nt=n_elements(ds.t)

; if we haven't specified a scale, make one.
IF keyword_set(silent) THEN BEGIN

moduestore=dblarr(nt)
modEstore=dblarr(nt)
modBstore=dblarr(nt)
Vparstore=dblarr(nt)
;print, "calc of max quantities (for scaling 3d arrows)"
FOR j=0,nt-1 DO BEGIN
 moduestore[j]=sqrt(total(ds.ue[*,j]*ds.ue[*,j]))
 modEstore[j]=sqrt(total(ds.E[*,j]*ds.E[*,j]))
 modbstore[j]=sqrt(total(ds.B[*,j]*ds.B[*,j]))
 Vparstore[j]=ds.vpar[j]
ENDFOR
;print, "..done."
maxue=max(moduestore)
maxE=max(modEstore)
maxB=max(modBstore)
maxVpar=max(Vparstore)
mvals=[maxue,maxE,maxB,maxVpar]

ENDIF ELSE BEGIN

maxue=mvals[0]
maxE=mvals[1]
maxB=mvals[2]
maxVpar=mvals[3]

particletrack, ii, n_c=32, tscl=0.0, op=1, floc=nom, lscl=lscl, xyzt

FOR j=0,nt-1,stride DO BEGIN

   tp=[[ds.x[j],ds.x[j]],[ds.y[j],ds.y[j]],[ds.z[j],ds.z[j]]]
   
   tt=ds.t[j]
   ystart=reform(tp[0,*])/lscl
   yp=reform(tp[0,*])/lscl
   ODEINT, ystart, eps, h1, string1, yp ; ODE mod works in NORMALISED COORDS
   xplot3d, yp(0,*)*lscl, yp(1,*)*lscl, yp(2,*)*lscl, COLOR=[0,0,255], /OVERPLOT

  ;xplot3D, tp[*,0], tp[*,1], tp[*,2], COLOR=[0,0,0], NAME='atest', SYMBOL=oSymbol, THICK=2, /OVERPLOT
  
  IF sqrt(total(ds.UE[*,j]*ds.UE[*,j])) gt 0 THEN BEGIN
   oModel=obj_new('IDLgrModel')
   oModel->add,mk_vector([ds.UE[0,j],ds.UE[1,j],ds.UE[2,j]]/maxue/dims,color=[255,0,0])
   ;oModel->add,mk_vector([ds.UE[0,j],ds.UE[1,j],ds.UE[2,j]]/maxvpar/dims,color=[255,0,0])
   oAR1 = obj_new('IDLgrSymbol', oModel)
   IF Obj_Valid(oAR1) THEN oAR1->SetProperty, Size=mysymsize
   xplot3D, tp[*,0], tp[*,1], tp[*,2], COLOR=[0,0,0], NAME='atest', SYMBOL=oAR1, THICK=2, /OVERPLOT
  ENDIF
  IF sqrt(total(ds.E[*,j]*ds.E[*,j])) gt 0 THEN BEGIN
   oModel=obj_new('IDLgrModel')
   oModel->add,mk_vector([ds.E[0,j],ds.E[1,j],ds.E[2,j]]/maxE/dims,color=[0,0,255])
   ;oModel->add,mk_vector([ds.E[0,j],ds.E[1,j],ds.E[2,j]]/sqrt(total(ds.E[*,j]*ds.E[*,j]))/dims,color=[0,0,255])
   oAR2 = obj_new('IDLgrSymbol', oModel)
   IF Obj_Valid(oAR2) THEN oAR2->SetProperty, Size=mysymsize*5.
   xplot3D, tp[*,0], tp[*,1], tp[*,2], COLOR=[0,0,0], NAME='atest', SYMBOL=oAR2, THICK=2, /OVERPLOT
  ENDIF
  IF sqrt(total(ds.B[*,j]*ds.B[*,j])) gt 0 THEN BEGIN
   oModel=obj_new('IDLgrModel')
   oModel->add,mk_vector([ds.B[0,j],ds.B[1,j],ds.B[2,j]]/sqrt(total(ds.B[*,j]*ds.B[*,j]))*ds.vpar[j]/maxvpar/dims,color=[0,255,0])
   ;oModel->add,mk_vector([ds.B[0,j],ds.B[1,j],ds.B[2,j]]/sqrt(total(ds.B[*,j]*ds.B[*,j]))/dims,color=[0,255,0])
   ;oModel->add,mk_vector([ds.B[0,j],ds.B[1,j],ds.B[2,j]]/maxB/dims,color=[0,255,0])
   oAR3 = obj_new('IDLgrSymbol', oModel)
   IF Obj_Valid(oAR3) THEN oAR3->SetProperty, Size=mysymsize
   xplot3D, tp[*,0], tp[*,1], tp[*,2], COLOR=[0,0,0], NAME='atest', SYMBOL=oAR3, THICK=2, /OVERPLOT
  ENDIF

ENDFOR
ENDELSE



END
