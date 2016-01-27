;changed to use double precision

;get the length of the file to read in
;probibly a more efficient way to do this, but it works!
spawn,'wc -l RV.dat',res
n=1UL
dum='blah'
reads,res,n,dum

A=dblarr(23,n)
openr,lun,'RV.dat',/get_lun
readf,lun,A
free_lun,lun
t=dblarr(n,/nozero)
t=transpose(A(0,*))
x=dblarr(n,/nozero)
x=transpose(A(1,*))
y=dblarr(n,/nozero)
y=transpose(A(2,*))
z=dblarr(n,/nozero)
z=transpose(A(3,*))
vpar=dblarr(n,/nozero)
vpar=transpose(A(4,*))
vperp2=dblarr(n,/nozero)
vperp2=transpose(A(5,*))
vdrift2=dblarr(n,/nozero)
vdrift2=transpose(A(6,*))
el1=dblarr(n,/nozero)
el1=transpose(A(7,*))
el2=dblarr(n,/nozero)
el2=transpose(A(8,*))
el3=dblarr(n,/nozero)
el3=transpose(A(9,*))
b1=dblarr(n,/nozero)
b1=transpose(A(10,*))
b2=dblarr(n,/nozero)
b2=transpose(A(11,*))
b3=dblarr(n,/nozero)
b3=transpose(A(12,*))
e1=dblarr(n,/nozero)
e1=transpose(A(13,*))
e2=dblarr(n,/nozero)
e2=transpose(A(14,*))
e3=dblarr(n,/nozero)
e3=transpose(A(15,*))
et=dblarr(n,/nozero)
et=transpose(A(16,*))
rhs1=dblarr(n,/nozero)
rhs1=transpose(A(17,*))
rhs2=dblarr(n,/nozero)
rhs2=transpose(A(18,*))
Vfx=dblarr(n,/nozero)
Vfx=transpose(A(19,*))
Vfy=dblarr(n,/nozero)
Vfy=transpose(A(20,*))
Vfz=dblarr(n,/nozero)
Vfz=transpose(A(21,*))
H=dblarr(n,/nozero)
H=transpose(A(22,*))
end
