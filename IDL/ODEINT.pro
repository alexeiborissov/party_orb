FUNCTION derivs2, x, y
common somename, epsilon
;------------------------------------------------------------------
;******************************************************************
; right hand side of field line ODE
;******************************************************************
!except=0
;
;JT B-field

Lscl=5.0d0
dydx = dblarr(3)
Lx=sqrt(2.0d0)*Lscl
Lz=Lscl
Ly=1.0d0

;Bx=-1.0d0-epsilon*(1.0d0-y[2]*y[2]/Lscl/Lscl)/((1.0d0+y[2]*y[2]/Lscl/Lscl)*(1.0d0+y[1]*y[1]/Lscl/Lscl))
;By=0.2d0
;Bz=y[0]
; Bx=-1.0d0-epsilon*(1.0d0-y[1]*y[1]/Lscl/Lscl)/((1.0d0+y[1]*y[1]/Lscl/Lscl)^2*(1.0d0+y[2]*y[2]/Lscl/Lscl))
; By=y[0]
; Bz=0.2d0
 Bx=0.2d0
 By=-1.0d0-epsilon*(1.0d0-y[2]*y[2]/Lz/Lz)/((1.0d0+y[2]*y[2]/Lz/Lz)^2*(1.0d0+y[0]*y[0]/Lx/Lx))
 Bz=y[1]/Ly
 
 dydx[0]=Bx
 dydx[1]=By
 dydx[2]=Bz

 
return, dydx
end

;##################################################################
; Runge-Kutta routines with adaptive step size control taken
; from Numerical Recipes
; rk4 from idl lib
;##################################################################
PRO RKQC, y, dydx,x,htry,eps,yscal,hdid,hnext,string1

common somename, epsilon
;----------------------------------------------------
; *****************************************************************
;  RKQC is the stepper routine which makes one adaptive step 
;  forward
; *****************************************************************
;-----------------------------------------------------------------
; assign some fixed values to variables
;-----------------------------------------------------------------

pgrow = -0.2
pshrnk= -0.25
fcor  = 1.0/15.0
safety = 0.9
errcon = 6.0e-4
ytemp = dblarr(3)
ysav  = ytemp
dysav = ysav

;-----------------------------------------------------------------
; save initial values for step since it is needed twice
;-----------------------------------------------------------------

xsav=x
ysav=y
dysav=dydx

;-----------------------------------------------------------------
; try stepsize htry first and define hh as stepsize for the 
; two half-steps
;-----------------------------------------------------------------

h = htry
jump1: hh= 0.5 *h


;-----------------------------------------------------------------
; take two half-steps first
;-----------------------------------------------------------------

ytemp=rk4(ysav,dysav,xsav,hh,string1)

x=xsav+hh

dydx=derivs2(x,ytemp)
y=rk4(ytemp,dydx,x,hh,string1)

;-----------------------------------------------------------------
; one step next
;-----------------------------------------------------------------

x=xsav+h
ytemp=rk4(ysav,dysav,xsav,h,string1)

;-----------------------------------------------------------------
; now evaluate the accuracy
;-----------------------------------------------------------------

errmax = 0.
ytemp = y-ytemp

errmax= max(abs(ytemp/yscal))

;-----------------------------------------------------------------
; scale error to required accuracy eps and reduce stepsize if 
; necessary
;-----------------------------------------------------------------

errmax = errmax/eps
if errmax gt 1.0 then begin
   h = safety*double(h)*(errmax^pshrnk)
   goto, jump1
endif else begin            
;-----------------------------------------------------------------
; if stepsize ok, compute stepsize for next step
;-----------------------------------------------------------------

   hdid = h                 
   if errmax gt errcon then begin
       hnext = safety*h*(errmax^pgrow)
   endif else begin
       hnext = 4.0*h
   endelse
endelse



;-----------------------------------------------------------------
; reduce truncation error to fifth order
;-----------------------------------------------------------------

y = y + ytemp*fcor

hnext=hnext

end ;_____of RKQC
;*****************************************************************
; ODEINT is the driver routine for the ODE solver
;*****************************************************************

PRO ODEINT, ystart,eps,h1,string1, yp
!except=0
common somename, epsilon
;common path, kount;, yp
;maxsteps = 1e4-2
maxsteps=5000
tiny = 1.0e-20

;-----------------------------------------------------------------
; we always start at arclength(=x) zero; set initial steplength to
; h1; and step counter to zero
;-----------------------------------------------------------------
x=0
h=h1
kount = 0
;-----------------------------------------------------------------
; assign initial values and store them in solution vector
;-----------------------------------------------------------------
y = ystart
yp=dblarr(3)
yp(*,0) = ystart(*)
;-----------------------------------------------------------------
; start integrating ODE
;-----------------------------------------------------------------
for iit=1,maxsteps do begin
	dydx = derivs2(x,y)
;-----------------------------------------------------------------
; for monitoring accuracy; see Numerical Recipes for explanation
;-----------------------------------------------------------------
  	yscal = abs(y) + abs(h*dydx) + tiny
;-----------------------------------------------------------------
; make one adaptive step forward and write result to yp
;-----------------------------------------------------------------
  RKQC,y, dydx,x,h,eps,yscal,hdid,hnext,string1
  kount = kount + 1
  yp=[[yp],[dblarr(3)]]
  yp(*,kount) = y(*)
  IF ((ABS(y(0)) GT 10.0) OR (ABS(y(1)) GT 10.0) OR (y(2) GT 40.0) OR (y(2) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (ABS(y(2)) GT 10.0) OR (y(1) GT 20.0) OR (y(1) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (ABS(y(1)) GT 10.0) OR (y(2) GT 20.0) OR (y(2) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (ABS(y(2)) GT 10.0) OR (y(1) GT 20.0) OR (y(1) LT 0.0)) THEN BEGIN
  ;IF ((ABS(y(0)) GT 10.0) OR (y(2) LT 0.0) OR (y(2) GT 20.0) OR (y(1) GT 20.0) OR (y(1) LT 0.0)) THEN BEGIN
   return
  ENDIF
  ;print, y(*)
  ;STOP

  ; limiting step size growth  
  IF abs(hnext) lt 0.1 THEN h = hnext
  
  
endfor
end;____________of ODEINT
