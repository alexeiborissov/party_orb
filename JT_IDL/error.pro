;plotting some quantities from particle runs across a separatrix
@gettempdata
LOADCT, 1

viewtime=1
viewmodB=0
viewmodE=0
viewgamma=0
viewU=0
viewRg=1
viewDUDT=0
viewDGAMMADT=0
viewDBXDX=0
viewDBYDX=0
viewDBZDX=0
viewDBXDY=0
viewDBYDY=0
viewDBZDY=0
viewDBXDZ=0
viewDBYDZ=0
viewDBZDZ=0
viewRX=0
viewRY=0
viewRZ=0
viewDRDTX=0
viewDRDTY=0
viewDRDTZ=0
viewEpar=0
viewDBDS=0
viewUEx=0
viewUEy=0
viewUEz=0
viewmodVf=0
viewVfx=0
viewVfy=0
viewVfz=0
viewmodgradB=0
viewBomodgradB=1

ds1=gettempdata(1, sub='e')
count=0
pmin=50
pmax=60
pstep=1
q=-1.60217653E-19
m=9.1093826E-31
bscl=0.001d0
omscl=q*bscl/m

;----EPS=1E-13-----;
;y1a=9.9e3
;y1b=9.95e3
;y2a=1.1725e4
;y2b=1.1775e4
;----EPS=1E-16-----;
y1a=1.23e4
y1b=1.25e4
y2a=1.46e4
y2b=1.48e4
tm=2.5e4

IF viewtime THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.t, xr=[0,tm], title='ds.t versus NSTP', xtitle='NSTP',  ytitle='ds.t' 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  print, i
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) OR (i eq 51) or (i eq 53) THEN mycol=200
  oplot, ds.NSTP, ds.t, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
; WRITE_PNG, "time.png", TVRD(/TRUE)
ENDIF

IF viewmodB THEN BEGIN
 modB1=sqrt(ds1.B[0,*]*ds1.B[0,*]+ds1.B[1,*]*ds1.B[1,*]+ds1.B[2,*]*ds1.B[2,*])
 window, count
 plot, ds1.NSTP, modB1, xr=[0,tm], title='|B| versus NSTP', xtitle='NSTP',  ytitle='|B|', yr=[0.01,1], /ylog
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  modB=sqrt(ds.B[0,*]*ds.B[0,*]+ds.B[1,*]*ds.B[1,*]+ds.B[2,*]*ds.B[2,*])
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, modB, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewmodE THEN BEGIN
 modE1=sqrt(ds1.E[0,*]*ds1.E[0,*]+ds1.E[1,*]*ds1.E[1,*]+ds1.E[2,*]*ds1.E[2,*])
 window, count
 plot, ds1.NSTP, modE1, xr=[0,tm], title='|E| versus NSTP', xtitle='NSTP', ytitle='|E|', yr=[1E-8,1E-4]
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  modE=sqrt(ds.E[0,*]*ds.E[0,*]+ds.E[1,*]*ds.E[1,*]+ds.E[2,*]*ds.E[2,*])
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, modE, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewgamma THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.gamma-1.0d0, title='gamma versus NSTP', xtitle='NSTP',  ytitle='gamma-1', xmargin=[10,14], xr=[0,tm] ;, yr=[3.9138e-8,3.9140e-8] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.gamma-1.0, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewU THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.u, title='U versus NSTP', xtitle='NSTP',  ytitle='U', xmargin=[10,14], yr=[-100,100], xr=[0,tm]
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.u, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewRG THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.rg, title='Rg versus NSTP', xtitle='NSTP',  ytitle='Rg',$
;  xmargin=[10,14], yr=[0.001,0.02]*0.001,  xr=[0,tm]
  xmargin=[10,14], yr=[1e-7,1e-5],  xr=[0,tm], /ylog
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) OR (i eq 51) or (i eq 53) THEN mycol=200
  oplot, ds.NSTP, ds.rg, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDGAMMADT THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DGAMMADT, title='DGAMMADT versus NSTP', xtitle='NSTP',  ytitle='DGAMMADT', xmargin=[10,14],xr=[0,tm], yr=[-3E-9,3E-9] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DGAMMADT, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDUDT THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DUDT, title='DUDT versus NSTP', xtitle='NSTP',  ytitle='DUDT', xmargin=[10,14], xr=[0,tm], yr=[1000,1000] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DUDT, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBXDX THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDX[0,*], title='DBXDX versus NSTP', xtitle='NSTP',  ytitle='DBXDX', xmargin=[10,14], xr=[0,tm], yr=[0,3] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDX[0,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBYDX THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDX[1,*], title='DBYDX versus NSTP', xtitle='NSTP',  ytitle='DBYDX', xmargin=[10,14], xr=[0,tm], yr=[0,3] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDX[1,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBZDX THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDX[2,*], title='DBZDX versus NSTP', xtitle='NSTP',  ytitle='DBZDX', xmargin=[10,14], xr=[0,tm], yr=[-4,0] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDX[2,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBXDY THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDY[0,*], title='DBXDY versus NSTP', xtitle='NSTP',  ytitle='DBXDY', xmargin=[10,14], xr=[0,tm], yr=[-4,-1] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDY[0,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBYDY THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDY[1,*], title='DBYDY versus NSTP', xtitle='NSTP',  ytitle='DBYDY', xmargin=[10,14], xr=[0,tm], yr=[-3,0] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDY[1,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBZDY THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDY[2,*], title='DBZDY versus NSTP', xtitle='NSTP',  ytitle='DBZDY', xmargin=[10,14], xr=[0,tm], yr=[-2,1] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDY[2,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBXDZ THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDZ[0,*], title='DBXDZ versus NSTP', xtitle='NSTP',  ytitle='DBXDZ', xmargin=[10,14], xr=[0,tm], yr=[-3,0] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDZ[0,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBYDZ THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDZ[1,*], title='DBYDZ versus NSTP', xtitle='NSTP',  ytitle='DBYDZ', xmargin=[10,14], xr=[0,tm], yr=[-1,0] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDZ[1,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBZDZ THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DBDZ[2,*], title='DBZDZ versus NSTP', xtitle='NSTP',  ytitle='DBZDZ', xmargin=[10,14], xr=[0,tm], yr=[-1,1] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DBDZ[2,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewRX THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.R[0,*], title='Rx versus NSTP', xtitle='NSTP',  ytitle='Rx', xmargin=[10,14], xr=[0,tm], yr=[-1,0] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.R[0,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewRY THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.R[1,*], title='Ry versus NSTP', xtitle='NSTP',  ytitle='Ry', xmargin=[10,14], xr=[0,tm], yr=[-1,0] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.R[1,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewRZ THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.R[2,*], title='Rz versus NSTP', xtitle='NSTP',  ytitle='Rz', xmargin=[10,14], xr=[0,tm], yr=[0,1] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.R[2,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDRDTX THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DRDT[0,*], title='DRDTx versus NSTP', xtitle='NSTP',  ytitle='DRDTx', xmargin=[10,14], xr=[0,tm], yr=[-50,50] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DRDT[0,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDRDTY THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DRDT[1,*], title='Ry versus NSTP', xtitle='NSTP',  ytitle='DRDTy', xmargin=[10,14], xr=[0,tm], yr=[-50,50] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DRDT[1,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDRDTZ THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.DRDT[2,*], title='DRDT|z versus NSTP', xtitle='NSTP',  ytitle='DRDTz', xmargin=[10,14], xr=[0,tm], yr=[-100,100] 
 loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.DRDT[2,*], col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewEpar THEN BEGIN
 window, count
 plot, ds1.NSTP, omscl*ds1.epar, xr=[0,tm], title='ds.t versus Epar', xtitle='NSTP',  ytitle='omscl*epar', yr=[-10,10] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, omscl*ds.epar, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewDBDS THEN BEGIN
 window, count
 plot, ds1.NSTP, ds1.dbds, xr=[0,tm], title='ds.t versus DBDS', xtitle='NSTP',  ytitle='DBDS', yr=[-3,3] 
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.dbds, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewUEx THEN BEGIN
 window, count
 modB1=sqrt(ds1.B[0,*]*ds1.B[0,*]+ds1.B[1,*]*ds1.B[1,*]+ds1.B[2,*]*ds1.B[2,*])
 UEx1=(ds1.E[1,*]*ds1.B[2,*]-ds1.E[2,*]*ds1.B[1,*])/modB1/modB1
 UEy1=(ds1.E[2,*]*ds1.B[0,*]-ds1.E[0,*]*ds1.B[2,*])/modB1/modB1
 UEz1=(ds1.E[0,*]*ds1.B[1,*]-ds1.E[1,*]*ds1.B[0,*])/modB1/modB1
 
 plot, ds1.NSTP, UEx1, xr=[0,tm], title='ds.t versus UEx', xtitle='NSTP',  ytitle='UEx', yr=[-1e-4, 1e-4]
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  modB=sqrt(ds.B[0,*]*ds.B[0,*]+ds.B[1,*]*ds.B[1,*]+ds.B[2,*]*ds.B[2,*])
  UEx=(ds.E[1,*]*ds.B[2,*]-ds.E[2,*]*ds.B[1,*])/modB/modB
  UEy=(ds.E[2,*]*ds.B[0,*]-ds.E[0,*]*ds.B[2,*])/modB/modB
  UEz=(ds.E[0,*]*ds.B[1,*]-ds.E[1,*]*ds.B[0,*])/modB/modB
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, UEx, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewUEy THEN BEGIN
 window, count
 modB1=sqrt(ds1.B[0,*]*ds1.B[0,*]+ds1.B[1,*]*ds1.B[1,*]+ds1.B[2,*]*ds1.B[2,*])
 UEx1=(ds1.E[1,*]*ds1.B[2,*]-ds1.E[2,*]*ds1.B[1,*])/modB1/modB1
 UEy1=(ds1.E[2,*]*ds1.B[0,*]-ds1.E[0,*]*ds1.B[2,*])/modB1/modB1
 UEz1=(ds1.E[0,*]*ds1.B[1,*]-ds1.E[1,*]*ds1.B[0,*])/modB1/modB1
 
 plot, ds1.NSTP, UEy1, xr=[0,tm], title='ds.t versus UEy', xtitle='NSTP',  ytitle='UEy', yr=[-1e-4, 1e-4]
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  modB=sqrt(ds.B[0,*]*ds.B[0,*]+ds.B[1,*]*ds.B[1,*]+ds.B[2,*]*ds.B[2,*])
  UEx=(ds.E[1,*]*ds.B[2,*]-ds.E[2,*]*ds.B[1,*])/modB/modB
  UEy=(ds.E[2,*]*ds.B[0,*]-ds.E[0,*]*ds.B[2,*])/modB/modB
  UEz=(ds.E[0,*]*ds.B[1,*]-ds.E[1,*]*ds.B[0,*])/modB/modB
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, UEy, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewUEz THEN BEGIN
 window, count
 modB1=sqrt(ds1.B[0,*]*ds1.B[0,*]+ds1.B[1,*]*ds1.B[1,*]+ds1.B[2,*]*ds1.B[2,*])
 UEx1=(ds1.E[1,*]*ds1.B[2,*]-ds1.E[2,*]*ds1.B[1,*])/modB1/modB1
 UEy1=(ds1.E[2,*]*ds1.B[0,*]-ds1.E[0,*]*ds1.B[2,*])/modB1/modB1
 UEz1=(ds1.E[0,*]*ds1.B[1,*]-ds1.E[1,*]*ds1.B[0,*])/modB1/modB1
 
 plot, ds1.NSTP, UEz1, xr=[0,tm], title='ds.t versus UEz', xtitle='NSTP',  ytitle='UEz', yr=[-1e-4, 1e-4]
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  modB=sqrt(ds.B[0,*]*ds.B[0,*]+ds.B[1,*]*ds.B[1,*]+ds.B[2,*]*ds.B[2,*])
  UEx=(ds.E[1,*]*ds.B[2,*]-ds.E[2,*]*ds.B[1,*])/modB/modB
  UEy=(ds.E[2,*]*ds.B[0,*]-ds.E[0,*]*ds.B[2,*])/modB/modB
  UEz=(ds.E[0,*]*ds.B[1,*]-ds.E[1,*]*ds.B[0,*])/modB/modB
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, UEz, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewmodVf THEN BEGIN
 modV1=sqrt(ds1.Vf[0,*]*ds1.Vf[0,*]+ds1.Vf[1,*]*ds1.Vf[1,*]+ds1.Vf[2,*]*ds1.Vf[2,*])
 window, count
 plot, ds1.NSTP, modV1, xr=[0,tm], title='|Vf| versus NSTP', xtitle='NSTP',  ytitle='|Vf|', yr=[1e-6,1e-2], /ylog
 ;loadct, 1
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  modV=sqrt(ds.Vf[0,*]*ds.Vf[0,*]+ds.Vf[1,*]*ds.Vf[1,*]+ds.Vf[2,*]*ds.Vf[2,*])
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, modV, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewmodgradB THEN BEGIN
 ;modV1=sqrt(ds1.Vf[0,*]*ds1.Vf[0,*]+ds1.Vf[1,*]*ds1.Vf[1,*]+ds1.Vf[2,*]*ds1.Vf[2,*])
 window, count
 plot, ds1.NSTP, ds1.modgradB, xr=[0,tm], title='|gradB| versus NSTP', xtitle='NSTP',  ytitle='|gradB|', yr=[0.1,100], /ylog
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  ;modV=sqrt(ds.Vf[0,*]*ds.Vf[0,*]+ds.Vf[1,*]*ds.Vf[1,*]+ds.Vf[2,*]*ds.Vf[2,*])
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) THEN mycol=220
  oplot, ds.NSTP, ds.modgradB, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
ENDIF

IF viewBomodgradB THEN BEGIN
 modB1=sqrt(ds1.B[0,*]*ds1.B[0,*]+ds1.B[1,*]*ds1.B[1,*]+ds1.B[2,*]*ds1.B[2,*])
 window, count
 plot, ds1.NSTP, modB1/ds1.modgradB, xr=[0,tm], title='|B|/|gradB| versus NSTP', xtitle='NSTP',  ytitle='|B|/|gradB|', yr=[0.001,1], /ylog
 FOR i=pmin,pmax,pstep DO BEGIN
  ds=gettempdata(i, sub='e')
  modB=sqrt(ds.B[0,*]*ds.B[0,*]+ds.B[1,*]*ds.B[1,*]+ds.B[2,*]*ds.B[2,*])
  mycol=254
  ;IF (i eq 50) THEN mycol=160
  IF (i eq 50) OR (i eq 51) or (i eq 53) THEN mycol=200
  oplot, ds.NSTP, modB/ds.modgradB, col=mycol
 ENDFOR
 ;oplot, [y1a,y1a], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y1b,y1b], [-1e8,1e8], col=220, linestyle=2
 ;oplot, [y2a,y2a], [-1e8,1e8], col=160, linestyle=2
 ;oplot, [y2b,y2b], [-1e8,1e8], col=160, linestyle=2
 count=count+1
; WRITE_PNG, "MODBoverMODGRADB.png", TVRD(/TRUE)
ENDIF

END
