function getnstpdata, ii, n_c=n_c, wkdir=wkdir
on_error, 2

IF n_elements(wkdir) eq 0 THEN wkdir="Data/DataR/"
IF n_elements(n_c) eq 0 THEN n_c=8

files=FINDFILE(wkdir+'d*e.tmp',count=count)
IF count eq 0 THEN BEGIN
 PRINT, "ERROR: no d*.tmp files found"
 return, 0
ENDIF
IF (ii lt 1) or (ii gt count) THEN BEGIN
 print, "ERROR: data files from ii=1->"+STRING(COUNT,FORMAT='(i4.4)')$
 +string(ii,format='("; ",i4.4," outside range!")')
 return, 0
ENDIF

 filename=string(wkdir, ii,format='(a,"d",i8.8,"e.tmp")')
 tstring=string(wkdir, ii,format='("wc -l ",a,"d",i8.8,"e.tmp")')
 spawn,tstring,res
 n_data=1UL
 dum='blah'
 reads,res,n_data,dum
 A=dblarr(n_c,n_data)
 openr,lun,filename,/get_lun
 readf,lun,A
 Free_lun,lun
 nstp=transpose(A(0,*))
 t=transpose(A(1,*))
 h=transpose(A(2,*))
 dudt=transpose(A(3,*))
 dxdt=transpose(A(4,*))
 dydt=transpose(A(5,*))
 dzdt=transpose(A(6,*))
 dgammadt=transpose(A(7,*))
thedata={name:filename, t:t,nstp:nstp,h:h,dudt:dudt,dxdt:dxdt,dydt:dydt,dzdt:dzdt,dgammadt:dgammadt}

return, thedata

END
