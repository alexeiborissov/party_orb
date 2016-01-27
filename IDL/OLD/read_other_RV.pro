;changed to use double precision

;get the length of the file to read in
;probibly a more efficient way to do this, but it works!
spawn,'wc -l ~/non_rel_0_4.2_5.5_alpha_160.4/RV.dat',res
n=1UL
dum='blah'
reads,res,n,dum

A=dblarr(23,n)
openr,lun,'~/non_rel_0_4.2_5.5_alpha_160.4/RV.dat',/get_lun
readf,lun,A
free_lun,lun
o_t=dblarr(n,/nozero)
o_t=transpose(A(0,*))
o_x=dblarr(n,/nozero)
o_x=transpose(A(1,*))
o_y=dblarr(n,/nozero)
o_y=transpose(A(2,*))
o_z=dblarr(n,/nozero)
o_z=transpose(A(3,*))
o_vpar=dblarr(n,/nozero)
o_vpar=transpose(A(4,*))
o_vperp2=dblarr(n,/nozero)
o_vperp2=transpose(A(5,*))
o_vdrift2=dblarr(n,/nozero)
o_vdrift2=transpose(A(6,*))
o_el1=dblarr(n,/nozero)
o_el1=transpose(A(7,*))
o_el2=dblarr(n,/nozero)
o_el2=transpose(A(8,*))
o_el3=dblarr(n,/nozero)
o_el3=transpose(A(9,*))
o_b1=dblarr(n,/nozero)
o_b1=transpose(A(10,*))
o_b2=dblarr(n,/nozero)
o_b2=transpose(A(11,*))
o_b3=dblarr(n,/nozero)
o_b3=transpose(A(12,*))
o_e1=dblarr(n,/nozero)
o_e1=transpose(A(13,*))
o_e2=dblarr(n,/nozero)
o_e2=transpose(A(14,*))
o_e3=dblarr(n,/nozero)
o_e3=transpose(A(15,*))
o_et=dblarr(n,/nozero)
o_et=transpose(A(16,*))
o_rhs1=dblarr(n,/nozero)
o_rhs1=transpose(A(17,*))
o_rhs2=dblarr(n,/nozero)
o_rhs2=transpose(A(18,*))
o_Vfx=dblarr(n,/nozero)
o_Vfx=transpose(A(19,*))
o_Vfy=dblarr(n,/nozero)
o_Vfy=transpose(A(20,*))
o_Vfz=dblarr(n,/nozero)
o_Vfz=transpose(A(21,*))
o_H=dblarr(n,/nozero)
o_H=transpose(A(22,*))
end
