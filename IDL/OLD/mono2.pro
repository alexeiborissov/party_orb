c1=-0.15*10.^(6)
cc=1.65; 0.4
;;;;;;;;; the transformation is in dimensionless variables ;;;
L=10.*10^6 ; 10^7 m = 10 Mm
Lvert=10.*10^6
d=L
Tscl=10; 100.
Vscl=L/Tscl
B0=0.01
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
a1=0.9
a2=0.9
esp=1.0
;;;;;;;  Grid ;;;;;;;;;;;;;;;;;;;;;;;;
n=500
n1=n/2
n2=n/2
xmin=0.01
xmax=2.
ymin=0.01
ymax=5. 
x1= xmin+findgen(n1)*(xmax-xmin)/(n1-1)
x2=-xmax+findgen(n2)*(-xmin+xmax)/(n2-1)
x=[ x2, x1]
y1= ymin+findgen(n1)*(ymax-ymin)/(n1-1)
y2=-ymax+findgen(n2)*(-ymin+ymax)/(n2-1)
y=[ y2, y1]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
r=fltarr(n,n)
mu=fltarr(n,n)
A=fltarr(n,n)
el=fltarr(n,n)
erase
t1=1.
t2=10.
tp1=1.05 ;100  ; 1.05
tp2=2.0  ;1000.1
;loadCT,3
deltat=1
deltatp=((tp2-tp1)/(t2-t1))*deltat
print,'tp1=',tp1
print,'tp2=',tp2
print,'deltatp=',deltatp
;openr,lun,'/user/paolo/simulations/part_orb/paulwood/monopole/dataold.dat',$
;   /get_lun
;openr,lun,'/user/paolo/simulations/part_orb/paulwood/monopole/drift/RV.dat',/get_lun
;dati=fltarr(13,5857)
;ndati=5857
;readf,lun,dati
;free_lun,lun
;int=0
;;;;;;;;;;;;;; TEMPORAL LOOP ;;;;;;;
;for k=1,ndati do begin
for t=t1,t2,deltat do begin
tp=tp1+((tp2-tp1)/(t2-t1))*(t-t1)
;tp=dati(0,k-1)
;;;; tanh transformation ;;;;;
;y1a=cc*tp*tanh(y1/cc/tp)
;ya=[y2, y1a]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; log transformation
;y1a=(cc*tp)^esp*alog(1.+(1./((cc*tp)^esp))*y1)
;ya=[y2, y1a]
;;;;;;; log transformation more complicated ;;;;;;
y1a=(cc*tp)^esp*alog(1.+(1./((cc*tp)^esp))*y1)*(1+tanh((y1-Lvert/L)*a1))*0.5 + $
     (1-tanh((y1-Lvert/L)*a2))*0.5*y1
ya=[y2, y1a]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;dipole;;;;;;;;;;;;;;
;for i=0,n-1 do begin
; for j=0,n-1 do begin
;  r(i,j)=sqrt(x(i)^2+ya(j)^2)
;  mu(i,j)=ya(j)/r(i,j)
; endfor
;endfor
;A=mu^2/r
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;monopole;;;;;;;;;;;;
for i=0,n-1 do begin
 for j=0,n-1 do begin
scra1_y = ya(j)+d/L
scra1_x = x(i)+1./2.
scra2_y = ya(j)+d/L
scra2_x = x(i)-1./2.
A(i,j)=c1*(atan(scra1_y,scra1_x)-atan(scra2_y,scra2_x))
endfor
endfor
;;;;;;;;;;;
;for i=0,n-1 do begin

; for j=0,n/2-1 do begin
; el(i,j)=0.
; endfor

; for j=n/2,n-1 do begin
;el(i,j)= -1.*(c1*(1./2.*(cc*tp)^esp*esp*alog(1.+ya(j)/(cc*tp)^esp)* $
;(1.+tanh((ya(j)-Lvert/L)*a1))/tp-1./2.*ya(j)*esp* $
;(1.+tanh((ya(j)-Lvert/L)*a1))/(tp*(1.+ya(j)/(cc*tp)^esp)))/ $
;((x(i)+1./2.)*(1.+(1./2.*(cc*tp)^esp*alog(1.+ya(j)/(cc*tp)^esp)* $
;(1.+tanh((ya(j)-Lvert/L)*a1))+1./2.*(1.-tanh((ya(j)-Lvert/L)*a2))*ya(j)+ $
;d/L)^2/(x(i)+1./2.)^2))-c1*(1./2.*(cc*tp)^esp*esp*alog(1.+ya(j)/ $
;(cc*tp)^esp)*(1.+tanh((ya(j)-Lvert/L)*a1))/tp-1./2.*ya(j)*esp* $
;(1.+tanh((ya(j)-Lvert/L)*a1))/(tp*(1.+ya(j)/(cc*tp)^esp)))/ $
;((x(i)-1./2.)*(1.+(1./2.*(cc*tp)^esp*alog(1.+ya(j)/ $
;(cc*tp)^esp)*(1.+tanh((ya(j)-Lvert/L)*a1))+1./2.* $
;(1.-tanh((ya(j)-Lvert/L)*a2))*ya(j)+d/L)^2/(x(i)-1./2.)^2)))/Tscl

;el(i,j)= (1./Vscl/B0)*el(i,j)
;endfor

;endfor
;;;;;;;;;;;;; PLOTTING ;;;;;;;;;;;;
;print,max(el),min(el)
erase
valmin=1.5e4 ;0.44 ; 0.48, c1=-1
valmax=7.90e4 ;0.55 ; 0.49, c1=-1 
ncont=20
vals=valmin+findgen(ncont)*(valmax-valmin)/(ncont-1.)
;surface,Ey,x,z,ax=30,az=45,xrange=[-500,500],zrange=[-1,0.5],xstyle=1,zstyle=1
;if (tp eq 1.05) then begin
;LoadCt,5
;set_plot,'ps'
;device,filename='cmt1.ps',/encapsulated
;shade_surf,el,x,y,ax=90,az=0,xrange=[-2,2],yrange=[0,5],xstyle=1,ystyle=1
contour,A,x,y,level=vals,xrange = [-2,2], yrange =[0,5]
; xrange = [-1,1], yrange=[1,3.5] ;,c_annotation=;vals  ;xrange=[-2,2],yrange=[0,5] ,c_annotation=vals
;oplot,dati(1,*),dati(2,*),psym=4
;plots,dati(1,int),dati(2,int),psym=4
;int=int+1
;device,/close
;endif
;set_plot,'X'
wait,0.
endfor
end
