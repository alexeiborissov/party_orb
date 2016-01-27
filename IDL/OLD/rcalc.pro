R1=fltarr(3)
R2=fltarr(3)
rsteps=intarr(3)

R1[0]=-0.5
R2[0]=0.5
RSTEPS[0]=8

R1[1]=0.0
R2[1]=0.0
RSTEPS[1]=1

R1[2]=0.0d0
R2[2]=0.0d0
RSTEPS[2]=1


NPOS=RSTEPS[0]*RSTEPS[1]*RSTEPS[2]

IF RSTEPS[0] EQ 1 then dx=1.0 ELSE dx=1.0/float(RSTEPS[0]-1)
IF RSTEPS[1] EQ 1 then dy=1.0 ELSE dy=1.0/float(RSTEPS[1]-1)
IF RSTEPS[2] EQ 1 then dz=1.0 ELSE dz=1.0/float(RSTEPS[2]-1)
ds=[dx,dy,dz]

lgrid=[[R2[0]-R1[0]],[R2[1]-R1[1]],[R2[2]-R1[2]]]


count=0
FOR i=0,RSTEPS[0]-1 DO BEGIN
 FOR j=0,RSTEPS[1]-1 DO BEGIN
  FOR k=0,RSTEPS[2]-1 DO BEGIN
   posno=[i,j,k]
   count=count+1
   ;RSTART   = R1+(R2-R1)*((posno*1.0d0)/(RSteps-1))
   ;print, count,i,j,k, RSTART, format='(i2,",",3i1 ,:,": ",3(f6.3, :, ", "))'
   
   RSTART=r1+lgrid*posno*ds 
   print, count,i,j,k, RSTART, format='(i2,",",3i1 ,:,": ",3(f6.3, :, ", "))'
  ENDFOR
 ENDFOR
ENDFOR
  


END
