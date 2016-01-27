function read_RVn,rvfilenumber



;changed to use double precision

;get the length of the file to read in
;probibly a more efficient way to do this, but it works!

spawn,string('wc -l ',rvfilenumber,'.dat'),res
n=1UL
dum='blah'
reads,res,n,dum


A=dblarr(23,n)
openr,lun,string(rvfilenumber,'.dat'),/get_lun
readf,lun,A
free_lun,lun

p = {t:dblarr(n,/nozero), x:dblarr(n,/nozero), y:dblarr(n,/nozero), z:dblarr(n,/nozero), vpar:dblarr(n,/nozero), vperp2:dblarr(n,/nozero), vdrift2:dblarr(n,/nozero), el1:dblarr(n,/nozero), el2:dblarr(n,/nozero), el3:dblarr(n,/nozero), b1:dblarr(n,/nozero), b2:dblarr(n,/nozero), b3:dblarr(n,/nozero), e1:dblarr(n,/nozero), e2:dblarr(n,/nozero), e3:dblarr(n,/nozero), et:dblarr(n,/nozero), rhs1:dblarr(n,/nozero), rhs2:dblarr(n,/nozero), Vfx:dblarr(n,/nozero), Vfy:dblarr(n,/nozero), Vfz:dblarr(n,/nozero), H:dblarr(n,/nozero),n:1UL, alpha:dblarr(n,/nozero) , b:dblarr(n,/nozero)}

p.t=transpose(A(0,*))
p.x=transpose(A(1,*))
p.y=transpose(A(2,*))
p.z=transpose(A(3,*))
p.vpar=transpose(A(4,*))
p.vperp2=transpose(A(5,*))
p.vdrift2=transpose(A(6,*))
p.el1=transpose(A(7,*))
p.el2=transpose(A(8,*))
p.el3=transpose(A(9,*))
p.b1=transpose(A(10,*))
p.b2=transpose(A(11,*))
p.b3=transpose(A(12,*))
p.e1=transpose(A(13,*))
p.e2=transpose(A(14,*))
p.e3=transpose(A(15,*))
p.et=transpose(A(16,*))
p.rhs1=transpose(A(17,*))
p.rhs2=transpose(A(18,*))
p.Vfx=transpose(A(19,*))
p.Vfy=transpose(A(20,*))
p.Vfz=transpose(A(21,*))
p.H=transpose(A(22,*))
p.n=n
;alphasigns=1-2*(p.vpar lt 0)
;p.alpha=atan(sqrt(p.e2),sqrt(p.e1)*alphasigns)*!radeg
p.alpha=acos(p.vpar/p.e1)*!radeg   ;e1 used to be the perp energy
p.b=sqrt(p.b1^2+p.b2^2+p.b3^2)

return, p

end
