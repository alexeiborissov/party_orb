function getsdata, ii, n_c=n_c, wkdir=wkdir
on_error, 2

IF n_elements(wkdir) eq 0 THEN wkdir="Data/DataR/"
IF n_elements(n_c) eq 0 THEN n_c=10


files=FINDFILE(wkdir+'sum?.dat',count=count)
IF count eq 0 THEN BEGIN
 PRINT, "ERROR: no sum*.dat files found"
 return, 0
ENDIF

t0=0.0d0
x0=0.0d0
y0=0.0d0
z0=0.0d0
ek0=0.0d0
t1=0.0d0
x1=0.0d0
y1=0.0d0
z1=0.0d0
ek1=0.0d0
 
  filename=string(wkdir, ii,format='(a,"sum",i1,".dat")')
  tstring=string(wkdir, ii,format='("wc -l ",a,"sum",i1,".dat")')
  spawn,tstring,res
  n_data=1UL
  dum='blah'
  reads,res,n_data,dum
  A=dblarr(n_c,n_data)
  openr,lun,filename,/get_lun
  readf,lun,A
  Free_lun,lun
  t0=transpose(A[0,*])
  x0=transpose(A[1,*])
  y0=transpose(A[2,*])
  z0=transpose(A[3,*])
  ek0=transpose(A[4,*])
  t1=transpose(A[5,*])
  x1=transpose(A[6,*])
  y1=transpose(A[7,*])
  z1=transpose(A[8,*])
  ek1=transpose(A[9,*])
   
thedata={name: filename, t0:t0,x0:x0,y0:y0,z0:z0,ek0:ek0,t1:t1,x1:x1,y1:y1,z1:z1,ek1:ek1}

return, thedata

END
