function getndata, ii, n_c=n_c, wkdir=wkdir, muBsquared=muBsquared, fields=fields, $
drifts=drifts, rdotperp=rdotperp, ke=ke, vpar=vpar,all=all, gyro=gyro, Epar=Epar
on_error, 2

IF n_elements(wkdir) eq 0 THEN wkdir="Data/DataN/"
IF n_elements(n_c) eq 0 THEN n_c=29

files=FINDFILE(wkdir+'RV*.dat',count=count)
IF count eq 0 THEN BEGIN
 PRINT, "ERROR: no RV*.dat files found"
 return, 0
ENDIF
IF (ii lt 1) or (ii gt count) THEN BEGIN
 print, "ERROR: data files from ii=1->"+STRING(COUNT,FORMAT='(i4.4)')$
 +string(ii,format='("; ",i4.4," outside range!")')
 return, 0
ENDIF

 filename=string(wkdir, ii,format='(a,"RV",i8.8,".dat")')
 tstring=string(wkdir, ii,format='("wc -l ",a,"RV",i8.8,".dat")')
 spawn,tstring,res
 n_data=1UL
 dum='blah'
 reads,res,n_data,dum
 A=dblarr(n_c,n_data)
 openr,lun,filename,/get_lun
 readf,lun,A
 Free_lun,lun
 t=transpose(A(0,*))
 x=transpose(A(1,*))
 y=transpose(A(2,*))
 z=transpose(A(3,*))

thedata={name: filename, t:t,x:x,y:y,z:z}

IF keyword_set(vpar) or keyword_set(all) THEN BEGIN
 vpar=transpose(A(4,*))
 thedata=create_struct(thedata,{vpar:vpar})
ENDIF
IF keyword_set(muBsquared) or keyword_set(all) THEN BEGIN
 muBsquared=transpose(A(5,*))
 thedata=create_struct(thedata,{muBsquared:muBsquared})
ENDIF
IF keyword_set(rdotperp) or keyword_set(all) THEN BEGIN
 rdotperp=transpose(A(6,*))
 thedata=create_struct(thedata,{rdotperp:rdotperp})
ENDIF
IF keyword_set(fields) or keyword_set(all) THEN BEGIN
 E=transpose([[transpose(A(7,*))],[transpose(A(8,*))],[transpose(A(9,*))]])
 B=transpose([[transpose(A(10,*))],[transpose(A(11,*))],[transpose(A(12,*))]])
 thedata=create_struct(thedata,{B:B,E:E})
ENDIF
IF keyword_set(ke) or keyword_set(all) THEN BEGIN
 e2=transpose(A(14,*))
 e3=transpose(A(15,*))
 ke=transpose(A(16,*))
 thedata=create_struct(thedata,{e2:e2,e3:e3,ke:ke})
ENDIF
IF keyword_set(drifts) or keyword_set(all) THEN BEGIN
 Ue=transpose([[transpose(A(23,*))],[transpose(A(24,*))],[transpose(A(25,*))]]) 
ENDIF
IF keyword_set(gyro) or keyword_set(all) THEN BEGIN
 gyrofreq=transpose(A(26,*))
 gyroperiod=transpose(A(27,*)) 
 gyrorad=transpose(A(28,*))
 thedata=create_struct(thedata,{gyrof:gyrofreq,gyrop:gyroperiod,gyror:gyrorad})
ENDIF
IF keyword_set(Epar) or keyword_set(all) THEN BEGIN
 Epar=transpose(A(22,*))
 thedata=create_struct(thedata,{Epar:Epar})
ENDIF

return, thedata

END
