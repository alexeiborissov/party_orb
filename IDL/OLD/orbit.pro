;KG (29/1/08): changed to use double precision variables

;;;wee use this routine to produce plots of 2d x-y orbits
Tscl=100
n_c=23
n_data=13475
A=dblarr(n_c,n_data)
openr,lun,'../Data/RV00000001.dat',/get_lun
;openr,lun,'RV.dat',/get_lun
readf,lun,A
Free_lun,lun
t=100*transpose(A(0,*))
x=transpose(A(1,*))
y=transpose(A(2,*))
z=transpose(A(3,*))
;set_plot,'ps'
;device,filename='xy.eps',/encapsulated
window, 0
plot,x(0:1000),y(0:1000),yrange=[1.4,2.5],xtitle='x',ytitle='y',charsize=1.5
;device,/close
;set_plot,'X'
end
