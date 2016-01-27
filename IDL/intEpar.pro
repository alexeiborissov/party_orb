; attempting to manually integrate Epar and compare with the change in KE

eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

rnrt='non-rel: '
rnr='nrel_'
writeout=1
imloc="images/"+rnr

thisLetter = "143B
char_gamma = '!4' + String(thisLetter) + '!X'
thisLetter2 = "154B

char_mu = '!4' + String(thisLetter2) + '!X'
thisLetter3 = "104B
char_DELTA = '!4' + String(thisLetter3) + '!X'
char_rperp='R!d!9x!n!3'
char_vpar='V!d!9#!n!3'


npts=16*16*5	; no of particles
nom="Data/"	; location of data
nomN="Data/DataN"; location of data
nomR="Data/DataR"; location of data

;i=43: close to y=0, far in x, benchmark for bouncing
;i=95: far in x, y and positive in z, travel upwards
;i=363: y=0, closer in x (see how rapidly electric field increases?
;i=603: closest particle to separator (see how rapidly electric field increases?

loadct, 39, /silent
 
is=43	; starting particle number
ie=43	; ending particle number
 
FOR i=is,ie DO BEGIN
 ds=getndata(i,/fields, /ke,/vpar, /gyro)
 ;ds=getrdata(i,/fields, /ke,/vpar, /gyro)
 ns=n_elements(ds.x)
 print, ds.x(0), ds.y(0), ds.z(0), format='("inipos=[",e9.2,",",e9.2,",",e9.2,"]")'
 tmax=max(ds.t)

 IF (i eq is) THEN BEGIN    ; on first go, do a 3d plot
  ;particletrack, i, n_c=32, tscl=0.0, lscl=1e6, op=0, floc=nom, xyzt, linecol=[255,0,0], zscl=1e6, /bouncesymb
  particletrack, i, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /mm, myxyrange=[-100,100], myzrange=[-60,60], /bsymb, mytitle=string(rnrt,i,format='(a,"particle path ", i3)')
 ENDIF ELSE BEGIN   	    ; else overplot
  ;particletrack, I, n_c=32, tscl=0.0, lscl=1e6, op=1, floc=nom, xyzt, linecol=[255,0,0], zscl=1e6, /bouncesymb
  particletrack, i, /op, floc=nom, xyzt, lcol=[255,0,0], zscl=1e6, /mm, myxyrange=[-100,100], myzrange=[-60,60], /bsymb
 ENDELSE

 modB=sqrt((ds.b[0,*]*ds.b[0,*]+ds.b[1,*]*ds.b[1,*]+ds.b[2,*]*ds.b[2,*]))
 Epar=reform((ds.E[0,*]*ds.B[0,*]+ds.E[1,*]*ds.B[1,*]+ds.E[2,*]*ds.B[2,*])/modB)

 bps=where(ds.vpar(0:ns-2)*ds.vpar(1:ns-1) lt 0)    ; bounce points where vpar changes sign
 nb=n_elements(bps)

 lb=0
 
 tr=0.5*(ds.t(lb:ns-2)+ds.t(lb+1:ns-1))
 nEpar=Epar(lb:ns-2)

 delx=ds.x(lb:ns-2)-ds.x(lb+1:ns-1) 	    	; delta_x
 dely=ds.y(lb:ns-2)-ds.y(lb+1:ns-1) 	    	; delta_y
 delz=ds.z(lb:ns-2)-ds.z(lb+1:ns-1)  	    	; delta_z

 delr=sqrt(delx*delx+dely*dely+delz*delz)   	; radial distance from one position to next
 
 trap=0.5d*(Epar(lb:ns-2)+Epar(lb+1:ns-1))  	; calculate the trapezoidal thing (to be integrated)
   
 deltake=ds.ke[lb+1:ns-1]-ds.ke[lb:ns-2]
 window, 0, ysize=350, xsize=850, xpos=1000, ypos=1000
cs=1.8
 !p.background=255
 ;!p.multi=[0,1,3]
 ;----------------------------------------------------------------------------------;  
 kemax=max(ds.ke)
 Eparmax=max(abs(Epar-Epar(0)))
  
 plot, ds.t, ds.ke, col=0, psym=-1, xr=[0,tmax], yr=[1.,ROUND(1.5*kemax)], charsize=cs, $
 title=string(rnrT,i,format='(a,"KE of particle",i4, " vs E!d!9#!n!3")'), xtitle='t (s)', $
 XSTYLE=8, YSTYLE=8, XMARGIN=[10, 12], ytitle='KE [eV]'
 plot, ds.t, Epar-Epar(0), col=240, linestyle=2, psym=-1, xstyle=4, ystyle=4, XMARGIN=[10, 12], $
 charsize=cs, yr=[-1.5*Eparmax,1.5*Eparmax], xr=[0,tmax], /noerase
 AXIS, YAXIS=1, YSTYLE = 1, YRANGE = [-1.5*Eparmax,1.5*Eparmax], YTITLE = 'E!d!9#!n!3 [Vm!e-1!n]', col=240, charsize=cs
 AXIS, XAXIS=1, XSTYLE = 1, col=240, charsize=cs, XTICKFORMAT="(A1)", XR=[0,round(tmax)]
  
 legend, ['KE','E!d!9#!n!3','bounce point loc'], linestyle=[0,2,2], psym=[-1,-1,0],col=[0,240,75], /top, /right, textcolors=0, charsize=cs-0.25
 
 IF (nb gt 1) THEN BEGIN    ; if there are more than one:
  FOR j=0,nb-1 DO BEGIN
   tnew=quickbifur(ds.t[bps[j]], ds.t[bps[j]+1], ds.vpar[bps[j]], ds.vpar[bps[j]+1])
   oplot, [tnew,tnew], [-1,1.5*kemax], col=75, linestyle=2
  ENDFOR
 ENDIF
 IF (nb eq 1) THEN BEGIN    ; if there is only one (and bps isn't -1)
  IF (bps ge 0) THEN BEGIN
    tnew=quickbifur(ds.t[bps], ds.t[bps+1], ds.vpar[bps], ds.vpar[bps+1])
    oplot, [tnew,tnew], [-1,1.5*kemax], col=75, linestyle=2
  ENDIF
 ENDIF
xyouts, 0.15, 0.82, 'peak(KE)='+string(max(ds.ke),format='(g9.3, "eV")'), /normal, col=0, charsize=cs-0.25
IF writeout THEN WRITE_PNG, string(imloc,i,format='(a,"p",i3.3,"_KEEpar.png")'), TVRD(/TRUE) 
 
 ;----------------------------------------------------------------------------------; 
; window, 1, ysize=350, xsize=850, xpos=1000, ypos=450
; maxdke=max(abs(deltake))
; largevel=max(abs(ds.vpar),mloc)
; plot, tr, abs(deltake), yr=[0,1.5*maxdke], col=0, xr=[0,tmax], charsize=cs, XMARGIN=[10, 12], xstyle=8,$
; title=string(rnrT,char_DELTA,i,format='(A,"|",A,"ke| of particle ",i4, " vs integrated E!d!9#!n!3")'), xtitle='t (s)', ytitle="|"+char_DELTA+'ke| [eV]'
; AXIS, XAXIS=1, XSTYLE = 1, col=0, charsize=cs, XTICKFORMAT="(A1)", XR=[0,round(tmax)]
; oplot, tr, delr*trap*(-1), col=240, linestyle=2
; legend, ['|'+char_DELTA+'ke|','integrated Epar','bounce point loc'], linestyle=[0,2,2], col=[0,240,75],/top, /right, textcolors=0, charsize=cs-0.25
; 
;  ; can we overplot the bounce points?
; IF (nb gt 1) THEN BEGIN    ; if there are more than one:
;  FOR j=0,nb-1 DO BEGIN
;   tnew=quickbifur(ds.t[bps[j]],ds.t[bps[j]+1],ds.vpar[bps[j]],ds.vpar[bps[j]+1])
;   oplot, [tnew,tnew], [-1,1], col=75, linestyle=2
;  ENDFOR
; ENDIF
; IF (nb eq 1) THEN BEGIN    ; if there is only one (and bps isn't -1)
;  IF (bps ge 0) THEN BEGIN
;     tnew=quickbifur(ds.t[bps],ds.t[bps+1],ds.vpar[bps],ds.vpar[bps+1])
;     oplot, [tnew,tnew], [-1,1], col=75, linestyle=2
;     ;oplot, [0.5*(tr[bps]+tr[bps+1]),0.5*(tr[bps]+tr[bps+1])], [0,1], col=75, linestyle=2
;  ENDIF
; ENDIF
; ;xyouts, 0.15, 0.82, 'final ke='+string(ds.vpar(mloc),format='(g10.4, "ms!e-1!n")'), /normal, col=0, charsize=1.5
; IF writeout THEN WRITE_PNG, string(imloc,i,format='(a,"p",i3.3,"_dKEintEpar.png")'), TVRD(/TRUE)

 dsn=getndata(i,/ke,/vpar, /gyro)
 dsr=getrdata(i,/ke,/vpar, /gyro)


 window, 1, ysize=350, xsize=850, xpos=1000, ypos=450
 ;maxdke=max(abs(deltake))
 ;largevel=max(abs(ds.vpar),mloc)
 maxgrn=max(dsn.gyror)
 maxgrr=max(dsr.gyror)
 maxgr=1e3*max([maxgrn,maxgrr])
 maxtn=max(dsn.t)
 maxtr=max(dsr.t)
 tmax2=max([maxtn,maxtr])
 plot, dsn.t, 1e3*dsn.gyror, yr=[0,1.5*maxgr], col=0, xr=[0,tmax2], charsize=cs, XMARGIN=[10, 12], xstyle=8,$
 title=string(i,format='(" gyroradius of particle ",i4)'), xtitle='t (s)', ytitle='r!dg!n [mm]'
 AXIS, XAXIS=1, XSTYLE = 1, col=0, charsize=cs, XTICKFORMAT="(A1)", XR=[0,tmax2]
 oplot, dsr.t, 1e3*dsr.gyror, col=240, linestyle=2
 legend, ['non-rel r!dg!n','rel r!dg!n', 'mirror point(s)'], linestyle=[0,2,2], col=[0,240,75],/top, /right, textcolors=0, charsize=cs-0.25
 

  ; can we overplot the bounce points?
 IF (nb gt 1) THEN BEGIN    ; if there are more than one:
  FOR j=0,nb-1 DO BEGIN
   tnew=quickbifur(ds.t[bps[j]],ds.t[bps[j]+1],ds.vpar[bps[j]],ds.vpar[bps[j]+1])
   oplot, [tnew,tnew], [-1,1], col=75, linestyle=2
  ENDFOR
 ENDIF
 IF (nb eq 1) THEN BEGIN    ; if there is only one (and bps isn't -1)
  IF (bps ge 0) THEN BEGIN
     tnew=quickbifur(ds.t[bps],ds.t[bps+1],ds.vpar[bps],ds.vpar[bps+1])
     oplot, [tnew,tnew], [-1,1], col=75, linestyle=2
     ;oplot, [0.5*(tr[bps]+tr[bps+1]),0.5*(tr[bps]+tr[bps+1])], [0,1], col=75, linestyle=2
  ENDIF
 ENDIF
 ;xyouts, 0.15, 0.82, 'final ke='+string(ds.vpar(mloc),format='(g10.4, "ms!e-1!n")'), /normal, col=0, charsize=1.5
 IF writeout THEN WRITE_PNG, string(imloc,i,format='(a,"p",i3.3,"_mirror.png")'), TVRD(/TRUE)
 
 ;----------------------------------------------------------------------------------; 
window, 2, ysize=350, xsize=850, xpos=1000, ypos=0
 bmax=max(modB)
 plot, ds.t, modB, col=240, linestyle=0, psym=-1, xstyle=4, ystyle=4, XMARGIN=[10, 12], charsize=cs, yr=[0,1.5*bmax], xr=[0,tmax]
 plot, ds.t, ds.vpar/max(abs(ds.vpar)), col=0, psym=-1, xr=[0,tmax], yr=[-1,1],charsize=cs, /noerase, $
 title=string(rnrt,i,format='(a,"v!d!9#!n!3 of particle",i4, " vs |B|")'), xtitle='t (s)', ytitle='v!d!9#!n!3/max(v!d!9#!n!3)',  $
 XSTYLE=8, YSTYLE=8, XMARGIN=[10, 12]
 AXIS, YAXIS=1, YSTYLE = 1, YRANGE = [0,1.5*bmax], YTITLE = '|B| [T]', col=240, charsize=cs
 AXIS, XAXIS=1, XSTYLE = 1, col=240, charsize=cs, XTICKFORMAT="(A1)", XR=[0,100]
 oplot, [0,100], [0,0], linestyle=1, col=0
 legend, ['v!d!9#!n!3/max(v!d!9#!n!3)','|B|','bounce point loc'], linestyle=[0,0,2], psym=[-1,-1,0], col=[0,240,75], /top, /right, textcolors=0, charsize=cs-0.25

 IF (nb gt 1) THEN BEGIN    ; if there are more than one:
  FOR j=0,nb-1 DO BEGIN
   tnew=quickbifur(ds.t[bps[j]], ds.t[bps[j]+1], ds.vpar[bps[j]], ds.vpar[bps[j]+1])
   oplot, [tnew,tnew], [-1,1], col=75, linestyle=2
  ENDFOR
 ENDIF
 IF (nb eq 1) THEN BEGIN    ; if there is only one (and bps isn't -1)
  IF (bps ge 0) THEN BEGIN
    tnew=quickbifur(ds.t[bps], ds.t[bps+1], ds.vpar[bps], ds.vpar[bps+1])
    oplot, [tnew,tnew], [-1,1], col=75, linestyle=2
  ENDIF
 ENDIF
 largevel=max(abs(ds.vpar),mloc)
 xyouts, 0.15, 0.82, 'peak(v!d!9#!n!3)='+string(ds.vpar(mloc),format='(g10.4, "ms!e-1!n")'), /normal, col=0, charsize=cs-0.25
 IF writeout THEN WRITE_PNG, string(imloc,i,format='(a,"p",i3.3,"_modBvpar.png")'), TVRD(/TRUE)
 
 WAIT, 0.4
print, 'max(ke) diff (ken-ker):', max(dsn.ke)-max(dsr.ke)

print, 'max(r_g) diff (rg_r-rg_n):', maxgrr-maxgrn
ENDFOR



!p.multi=0
END
