function gettempdata, ii, n_c=n_c, wkdir=wkdir, sub=sub
on_error, 2

IF n_elements(wkdir) eq 0 THEN wkdir="Data/DataR/"
IF n_elements(n_c) eq 0 THEN n_c=38
IF n_elements(sub) eq 0 THEN sub='e'


efiles=FINDFILE(wkdir+'f*e.tmp',count=count)
IF count eq 0 THEN BEGIN
 PRINT, "ERROR: no f*.tmp files found"
 return, 0
ENDIF
IF (ii lt 1) or (ii gt count) THEN BEGIN
 print, "ERROR: data files from ii=1->"+STRING(COUNT,FORMAT='(i4.4)')$
 +string(ii,format='("; ",i4.4," outside range!")')
 return, 0
ENDIF

 filename=string(wkdir, ii,sub,format='(a,"f",i8.8,a,".tmp")')
 tstring=string(wkdir, ii,sub,format='("wc -l ",a,"f",i8.8,a,".tmp")')
 spawn,tstring,res
 n_data=1UL
 dum='blah'
 reads,res,n_data,dum
 A=dblarr(n_c,n_data)
 openr,lun,filename,/get_lun
 readf,lun,A
 Free_lun,lun

 NSTP=transpose(A(0,*))
 B=transpose([[transpose(A(1,*))],[transpose(A(2,*))],[transpose(A(3,*))]])
 DBDX=transpose([[transpose(A(4,*))],[transpose(A(5,*))],[transpose(A(6,*))]])
 DBDY=transpose([[transpose(A(7,*))],[transpose(A(8,*))],[transpose(A(9,*))]])
 DBDZ=transpose([[transpose(A(10,*))],[transpose(A(11,*))],[transpose(A(12,*))]])
 E=transpose([[transpose(A(13,*))],[transpose(A(14,*))],[transpose(A(15,*))]])
 DEDX=transpose([[transpose(A(16,*))],[transpose(A(17,*))],[transpose(A(18,*))]])
 DEDY=transpose([[transpose(A(19,*))],[transpose(A(20,*))],[transpose(A(21,*))]])
 DEDZ=transpose([[transpose(A(22,*))],[transpose(A(23,*))],[transpose(A(24,*))]])
 DRDT=transpose([[transpose(A(25,*))],[transpose(A(26,*))],[transpose(A(27,*))]])
 DUDT=transpose(A(28,*))
 DGAMMADT=transpose(A(29,*))
 Epar=transpose(A(30,*))
 t=transpose(A(31,*))
 R=transpose([[transpose(A(32,*))],[transpose(A(33,*))],[transpose(A(34,*))]])
 U=transpose(A(35,*))
 ;gamma=transpose(A(36,*))
; DBDS=transpose(A(37,*))
; Vf=transpose([[transpose(A(36,*))],[transpose(A(37,*))],[transpose(A(38,*))]])
 modgradB=transpose(A(36,*)) 
 rg=transpose(A(37,*))

thedata={name: filename, NSTP:NSTP,B:B,DBDX:DBDX,DBDY:DBDY,DBDZ:DBDZ,E:E,DEDX:DEDX,DEDY:DEDY,DEDZ:DEDZ, DRDT:DRDT, DUDT:DUDT, DGAMMADT:DGAMMADT, $
;Epar:Epar, t:t, R:R, U:U, gamma:gamma, DBDS:DBDS, rg:rg}
;Epar:Epar, t:t, R:R, U:U, Vf:Vf}
Epar:Epar, t:t, R:R, U:U,modgradB:modgradB, rg:rg}
return, thedata

END
