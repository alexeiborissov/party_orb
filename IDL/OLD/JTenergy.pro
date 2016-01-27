
common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

npts=16*16*5
;nom="../Data/"
;nom="../../SeptNoEfield/Data/"

window, 0

i=10
;for i=2,npts-2,2 DO BEGIN
 dsE=getdata(i,/ke,wkdir="../Data/")
 ds=getdata(i,/ke,wkdir="../../SeptNoEfield/Data/")
 
 plot, dsE.t, dsE.ek

 
;endfor


END
