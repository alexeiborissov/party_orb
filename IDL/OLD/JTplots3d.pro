;Pro PlotS3D

   ; Create some data.
   
;time = Findgen(200)



;x = Sin(time/2.5) / Exp(time/100) ; Damped Sine Curve
;y = Sin(time/2.5) / Exp(time/100) ; Damped Sine Curve
;z = Sin(time/5.0 * !DtoR) ; Function of the time
;x = Reverse(x) ; Reverse the x vector
loadct, 2


;minx = Min(xx, Max=maxx)
;miny = Min(yy, Max=maxy)

   ; Open a window. Load a yellow color for the plot.
   
Window, /Free, XSize=450, YSize=450
TVLCT, 255, 255, 0, 1

   ; Set up the 3D space. Draw axes for the plot.
   ; Save the 3D transformation in !P.T.
   
Surface, Dist(5), XRange=[-1, 1], YRange=[-1, 1], $
   ZRange=[-6,6], XStyle=1, YStyle=1, ZStyle=1, /NoData, /Save

   ; Plot the data in 3D space (use T3D keyword).
   
PlotS, xx, yy, zz, /T3D, Color=1

   ; Create an animation of the data over time.
   ; Let's do 50 frames.
   
XInterAnimate, Set=[250,250,250], /Showload
FOR j=1,249 DO BEGIN
   Surface, Dist(5), XRange=[-1, 1], YRange=[-1, 1], $
      ZRange=[-6,6], XStyle=1, YStyle=1, ZStyle=1, /NoData, /Save
   PlotS, xx((j-1)*4:j*4), yy((j-1)*4:j*4), zz((j-1)*4:j*4), /T3D, Color = 1
   ;PlotS, xx, yy, zz, /T3D, Color = 100
   XInterAnimate, Frame=j, Window=!D.Window
ENDFOR
XInterAnimate

END
