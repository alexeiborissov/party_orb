;n=8344
n=500
n1=0
A=fltarr(5,n)
openr,lun,'mbt.dat',/get_lun
readf,lun,A
free_lun,lun
t=transpose(A(0,n1:n-1))
B=transpose(A(1,n1:n-1))
B1=transpose(A(2,n1:n-1))  ;dpart_b_dt
B2=transpose(A(3,n1:n-1))  ;vpdbds
B3=transpose(A(4,n1:n-1))     ;uegradb
;dtot_b_dt=deriv(t,B)

;set_plot,'ps'
;device,filename='mbt.eps',/encapsulated
;plot,t,dtot_b_dt,psym=4;,xrange=[5,7],yrange=[-10^3,15*10^3], $
    ; charsize=1.5,xtitle='t (s)', $
    ; ytitle='dE!Dk!N/dt (eV/s)',ytickv=[ 0,  5*10^3,10^4,1.5*10^4], $
    ; ytickname=['0','5!U.!N10!U3!N','10!U4!N','1.5!U.!N10!U4!N'],ystyle=1
;plot,t,dpart_B_dt
;oplot,t,dpart_B_dpart_t+vpdbds+uegradb
;oplot,t,replicate(0,n)
;device,/close
;set_plot,'X'

close,/all
end
