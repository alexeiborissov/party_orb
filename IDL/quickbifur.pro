FUNCTION quickbifur, x1, x2, y1, y2

FOR i=0,10 DO BEGIN

 dx=x2-x1
 xnew=x1+0.5d0*dx
 ynew=y1+(y2-y1)/(x2-x1)*0.5*dx

 IF (y1*ynew le 0) THEN BEGIN
  x1=x1
  x2=xnew
  y1=y1
  y2=ynew
 ENDIF ELSE BEGIN
  x1=xnew
  x2=x2
  y1=ynew
  y2=y2
 ENDELSE 

;oplot, [xnew, xnew], [-1e10,1e10], col=150, thick=2

ENDFOR

return, xnew

END
