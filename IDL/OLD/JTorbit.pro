;KG (29/1/08): changed to use double precision variables

;;;wee use this routine to produce plots of 2d x-y orbits
Tscl=100
n_c=23

;i=1
for i=1,20 DO BEGIN
 
 ;  filename='../Data/RV00000001.dat'
filename=string(i,format='("../Data/RV",i8.8,".dat")')
tstring=string(i,format='("wc -l ../Data/RV",i8.8,".dat")')
;spawn,'wc -l ../Data/RV00000001.dat',res
spawn,tstring,res

n_data=1UL
dum='blah'
reads,res,n_data,dum

;n_data=13475
A=dblarr(n_c,n_data)
openr,lun,filename,/get_lun
readf,lun,A
Free_lun,lun
t=100*transpose(A(0,*))
x=transpose(A(1,*))
y=transpose(A(2,*))
z=transpose(A(3,*))
;set_plot,'ps'
;device,filename='xy.eps',/encapsulated

help, A
;window, 0
;plot,x,y, title='i', xr=[-1e-7,1e-7]
;plot,x(0:1000),y(0:1000),yrange=[1.4,2.5],xtitle='x',ytitle='y',charsize=1.5
;device,/close
;set_plot,'X'

;wait, 0.5
endfor

end
