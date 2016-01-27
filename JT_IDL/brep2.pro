;making Bourdin report (20 Nov 2014)
;checking individual particle properties.
party=25


;@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1
;@xplot3dJT
loadct, 4, /silent
tvlct, r, g, b, /get

f=32
filename='test.png'
ax=-60
az=30
filewidth=800
transparent=0.5

notbored=0

IF notbored THEN BEGIN
;IF keyword_set(filename) and keyword_set(transparent) THEN BEGIN
; type=strsplit(filename,".",/extract)
; newfilenames=[type[0]+'_a.'+type[1],type[0]+'_b.'+type[1]]
;ENDIF

print, 'getting data..'
;blah='../../laredata/seplength4/jcrit25_beta09/'
;blah='../../laredata/seplength4/jcrit20/'
blah='../../../../2014/bourdin/data/VAR378_B.xdr'
blah2='../../../../2014/bourdin/data/VAR378_J_ABS.xdr'
blah3='../../../../2014/bourdin/data/VAR378_E_PARALLEL.xdr'
;jnom=blah+string(f, format='("newJ", I2,".sav")')
;blah2=blah+'Data/'

;ds=getdata(41,/vx,/vy,/vz,/bx,/by,/bz)
;ds=getdata(f,/grid, wkdir="../../julie_dam_files/Data")
;ds=getdata(f,/grid, wkdir='../../laredata/seplength4/jcrit20/Data/')
;ds=getdata(f,/grid, wkdir=blah2)
restore, /verbose, filename=blah3

;grid=ds.grid
delx=x[1]-x[0]
dely=y[1]-y[0]
;delz=z[1]-z[0]
delz=z-shift(z,-1)
;x=grid.x
nx = n_elements(x)-1
dx = x(1)-x(0)
;y=grid.y
ny = n_elements(y)-1
dy = y(1)-y(0)
;z=grid.z
nz = n_elements(z)-1
dz=z-shift(z,-1)

lscl=1e6

xscale=(max(x)-min(x))/DOUBLE(nx)/lscl
yscale=(max(y)-min(y))/DOUBLE(ny)/lscl
zscale=(max(z)-min(z))/DOUBLE(nz)/lscl
zscale2=(max(z)-z(40))/DOUBLE(nz)/lscl
xt='x(Mm)'
yt='y(Mm)'
zt='z(Mm)'


xregion=where((x le 100e6) and (x ge 99e6))
yregion=where((y le 100.5e6) and (y ge 99.5e6))
zregion=where((z le 75.2e6) and (z ge 74.2e6))

;make isosurface of current
;theModel = Obj_New('IDLgrModel')
;theModel -> Scale, xscale,yscale,zscale	    	    	
;theModel -> Add, mk_iso(J_ABS,0.005,scol=[255,255,255], /low)	
;test=obj_new('IDLgrSymbol', theModel)
;if obj_valid(test) THEN test->setproperty, size=1
;if obj_valid(test) THEN test->SetProperty, COLOR = [0, 240, 0]
;xplot3d, [0,0], [0,0], [0,0], symbol=test, $
;xr=[0,max(x)]/lscl, yr=[0,max(y)]/lscl, zr=[0,max(z)]/lscl, ax=ax, az=az, filewidth=filewidth,$
;xtitle=xt, ytitle=yt, ztitle=zt

;J_ABS=J_ABS(*,*,40:*)

lee=60	;lowest z-gridpoint included

;study currents above 
theModel = Obj_New('IDLgrModel')
theModel -> Scale, xscale,yscale,zscale2	    	    	
theModel -> Add, mk_iso(E_PARALLEL(*,*,lee:250),0.75,scol=[0,0,255], /low, minval=0.05)	
test=obj_new('IDLgrSymbol', theModel)
if obj_valid(test) THEN test->setproperty, size=1
if obj_valid(test) THEN test->SetProperty, COLOR = [0, 0, 255]
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test, $
xr=[0,max(x)]/lscl, yr=[0,max(y)]/lscl, zr=[0,max(z)]/lscl, ax=ax, az=az, filewidth=filewidth,$
xtitle=xt, ytitle=yt, ztitle=zt

theModel2 = Obj_New('IDLgrModel')
theModel2 -> Scale, xscale,yscale,zscale2	    	    	
theModel2 -> Add, mk_iso(E_PARALLEL(*,*,lee:250),-0.75,scol=[255,0,0], maxval=-0.05)	
test2=obj_new('IDLgrSymbol', theModel2)
if obj_valid(test2) THEN test2->setproperty, size=1
if obj_valid(test2) THEN test2->SetProperty, COLOR = [255, 0, 0]
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test2, /overplot

theModel3 = Obj_New('IDLgrModel')
theModel3 -> Scale, xscale,yscale,zscale2	    	    	
theModel3 -> Add, mk_iso(E_PARALLEL(*,*,lee:250),0.05,scol=[158, 182, 255], /low, minval=0.0005d, transparency=0.3)	
test3=obj_new('IDLgrSymbol', theModel3)
if obj_valid(test3) THEN test3->setproperty, size=1
if obj_valid(test3) THEN test3->SetProperty, COLOR = [158, 182, 255]
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test3, /overplot

theModel4 = Obj_New('IDLgrModel')
theModel4 -> Scale, xscale,yscale,zscale2	    	    	
theModel4 -> Add, mk_iso(E_PARALLEL(*,*,lee:250),-0.05,scol=[255, 139, 100], maxval=-0.0005d, transparency=0.3)	
test4=obj_new('IDLgrSymbol', theModel4)
if obj_valid(test4) THEN test4->setproperty, size=1
if obj_valid(test4) THEN test4->SetProperty, COLOR = [255, 139, 100]
xplot3d, [0,0], [0,0], [z(lee),z(lee)]/lscl, symbol=test4, /overplot

;STOP

;npts=1
;nom="Data/"
;print, 'loading: ', nom
;ds=getrdata(1,/ke)
;inike=ds.ke(0)
;;loadct, 2
;estore=0
;kestore=dblarr(npts)
;startpt=dblarr(3,npts)
;for i=1,npts DO BEGIN
; tflag=0
; print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(party)
 startpt=[ds.x(0),ds.y(0),ds.z(0)]
; maxke=max(ds.ke)
; print, startpt[*,i-1]
;endfor
;cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

;ref=alog10(1.0e-12)
;alke=alog10(kestore-2.0d0)
;astore=fix((alke-ref)/max(alke-ref)*254)
;t1=ref
;t2=alog10(max(kestore)-2.0d0)

 ch = 0.01*dx 	    	    ;step size
 chmin = 0.001*dx	    	    ;min step size
 chmax = 0.5*dx	    	    ; max step size
 cepsilon = 1.0d-5  	    ; tolerance
 cmxline = 100000    	    ; max nsteps
 csz_stp = size(startpt)

 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]

 particletrack, party, floc=nom, xyzt, lcol=[0,0,0], /op, zscl=mylscl, /mm, lscl=lscl, /rel, /symb, /bsymb
 
restore, /verbose, filename=blah

 
   linef = tracefield_3D(startpt[*,0],b,x,y,z,ch,chmin,chmax,cepsilon,cmxline)
   ii = where(linef[0,*] gt -8000,nii)
   if (ii[0] ne -1) then xplot3d, linef(0,ii)/lscl, linef(1,ii)/lscl, linef(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT
   ;trace field line backward
   lineb = tracefield_3D(startpt[*,0],b,x,y,z,-ch,chmin,chmax,cepsilon,cmxline)
   ii = where(lineb[0,*] gt -8000,nii)
   if (ii[0] ne -1) then xplot3d, lineb(0,ii)/lscl, lineb(1,ii)/lscl, lineb(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT

 

; for i=2,npts DO BEGIN
;  particletrack, i, floc=nom, xyzt, lcol=[0,0,0], /op, zscl=mylscl, /mm, lscl=lscl, /rel
;  linef = tracefield_3D(startpt[*,i-1],b,x,y,z,ch,chmin,chmax,cepsilon,cmxline)
;  ii = where(linef[0,*] gt -8000,nii)
;  if (ii[0] ne -1) then xplot3d, linef(0,ii)/lscl, linef(1,ii)/lscl, linef(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT
;  ;trace field line backwards
;  lineb = tracefield_3D(startpt[*,i-1],b,x,y,z,-ch,chmin,chmax,cepsilon,cmxline)
;  ii = where(lineb[0,*] gt -8000,nii)
;  if (ii[0] ne -1) then xplot3d, lineb(0,ii)/lscl, lineb(1,ii)/lscl, lineb(2,ii)/lscl, COLOR=brightgreen, /OVERPLOT
; endfor
ENDIF
 ds=getrdata(party,/gyro, /fields, /vpar, /ke,/rel)
 ds0=ds
 ns=n_elements(ds.x)
 print, ds.x(0), ds.y(0), ds.z(0), format='("inipos=[",e9.2,",",e9.2,",",e9.2,"]")'
 tmax=max(ds.t)
 lb=0

 modB0=sqrt((ds.b[0,*]*ds.b[0,*]+ds.b[1,*]*ds.b[1,*]+ds.b[2,*]*ds.b[2,*]))
 ;Epar0=reform((ds.E[0,*]*ds.B[0,*]+ds.E[1,*]*ds.B[1,*]+ds.E[2,*]*ds.B[2,*])/modB0)
 tr0=0.5*(ds.t(lb:ns-2)+ds.t(lb+1:ns-1))
 ;nEpar0=Epar0(lb:ns-2)
 nvpar0=ds.vpar/max(abs(ds.vpar))
;STOP
 cs=2
 ss=2
 !p.background=255
 ;!p.multi=[0,1,3]
 ;----------------------------------------------------------------------------------;  
 kemax=max(ds0.ke)
 ;Eparmax=max(abs(Epar0))

;set_plot, 'ps'
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2

mypath='./img/'

 ;window, 2, ysize=900, xsize=1600
 ;psa, file='~/Dropbox/JT_separators/bourdin/report/figs/e_fullfield_EparKE.eps', xs=13.66, ys=7.68, /inches, /encaps, bits_per_pixel=24, /color
 psa, file=mypath+string(party, format='("e_",i2.2,"_EparKE.eps")'), xs=13.66, ys=7.68, /inches, /encaps, bits_per_pixel=24, /color 
 i=3
 loadct, 39, /silent
  
 plot, ds0.t, ds0.ke, col=0, psym=-4, xr=[1e-3,100], yr=[1E0,2E4], charsize=cs, $
 xtitle='t (s)', symsize=ss, XSTYLE=8, YSTYLE=8, XMARGIN=[10, 12], ytitle='KE [eV]', /xlog, /ylog, position=[.11,.18,.89,.9]
 
 ;plot, ds.t, ds.epar, col=240, psym=-4, xstyle=4, ystyle=4, XMARGIN=[10, 12], symsize=ss,$
 ;charsize=cs, yr=[-5E-2,5E-2], xr=[1e-3,100], /noerase, /xlog
 
 plot, ds.t, ds.Epar, col=240, psym=-4, xstyle=4, ystyle=4, XMARGIN=[10, 12], symsize=ss,$
 charsize=cs, yr=[-1e-2,1e-2], xr=[1e-3,100], /noerase, /xlog, position=[.11,.18,.89,.9]
 
 AXIS, YAXIS=1, YSTYLE = 1, YRANGE = [-2e-2,2e-2], YTITLE = 'E!d||!n!3 [Vm!e-1!n]', col=240, charsize=cs, /xlog
 AXIS, XAXIS=1, XSTYLE = 1, col=240, charsize=cs, XTICKFORMAT="(A1)", xr=[1e-3,100], /xlog
  
 ;legend, ["A","B","C","D"], linestyle=[0,0,0,0], psym=[-4,-7,-8,-1],col=[0,0,0,0], /bottom, /left, textcolors=0, charsize=cs-0.5, symsize=[ss,ss,ss-0.5,ss]
 pse

 
 ;----------------------------------------------------------------------------------; 
 ;window, 4, ysize=900, xsize=1600

 bmax=max(modB0)
 psa, file=mypath+string(party, format='("e_",i2.2,"_VB.eps")'), xs=13.66, ys=7.68, /inches, /encaps, bits_per_pixel=24, /color 
 plot, ds0.t, nvpar0, col=0, psym=-4, xr=[1e-3,100], yr=[-1,1], charsize=cs, $
 xtitle='t (s)', symsize=ss, $
 XSTYLE=8, YSTYLE=8, XMARGIN=[10, 12], ytitle='v!d||!n!3/max(v!d||!n!3)', /xlog, position=[.11,.18,.89,.9]
 
 plot, ds0.t, modB0, col=58, psym=-4, xstyle=4, ystyle=4, XMARGIN=[10, 12], symsize=ss, $
 charsize=cs, yr=[1.5E-4,3.5e-4], xr=[1e-3,100], /noerase, /xlog, position=[.11,.18,.89,.9]

 
 AXIS, YAXIS=1, YSTYLE = 1, YRANGE = [1.5E-4,3.5e-4], YTITLE = '|B| [ T ]', col=58, charsize=cs
 AXIS, XAXIS=1, XSTYLE = 1, col=58, charsize=cs, XTICKFORMAT="(A1)", xr=[1e-3,100], /xlog

 ;legend, ['A','B','C','D'], linestyle=[0,0,0,0], psym=[-4,-7,-8,-1],col=[0,0,0,0], /bottom, /left, textcolors=0, charsize=cs-0.5, symsize=[ss,ss,ss-0.5,ss]
 pse

; legend, ['A','B'], linestyle=[0,0], psym=[-4,-8],col=[0,0], /bottom, /left, textcolors=0, charsize=cs-0.5, symsize=[ss,ss-0.5]
 
; IF writeout THEN WRITE_PNG, string(imloc,format='(a,"VvB_nf_2eV_1deg.png")'), TVRD(/TRUE)
 ;IF writeout THEN WRITE_PNG, imloc+nom2, TVRD(/TRUE)

;-------------------------------------------------------------------------------------------------------;
 loadct, 39, /silent

 dsr0=getrdata(party,/ke,/vpar, /gyro)

;
; window, 9, ysize=900, xsize=1600
 psa, file=mypath+string(party, format='("e_",i2.2,"_gyro.eps")'), xs=13.66, ys=7.68, /inches, /encaps, bits_per_pixel=24, /color 
 
 maxgrr=max(dsr0.gyror)
 maxgr=1e3*maxgrr
 maxtr=max(max(dsr0.t))
 mintr=min(max(dsr0.t))
 ;tmax=max([maxtn,maxtr])
 tmax=maxtr
 plot, dsr0.t, dsr0.gyror*1e3, yr=[33,43], col=0, xr=[1e-3,100], charsize=cs, xstyle=8, psym=-4, symsize=ss-0.5, $
 ;title=tit3, xtitle='t (s)', ytitle='r!dg!n [mm]', position=[.11,.18,.93,.88], /xlog
 xtitle='t (s)', ytitle='r!dg!n [mm]', position=[.11,.18,.89,.9], /xlog, XMARGIN=[10, 12]
 AXIS, XAXIS=1, XSTYLE = 1, col=0, charsize=cs, XTICKFORMAT="(A1)", XR=[1e-3,tmax], /xlog

 pse

 ;print, "[A,B,C,D]"
 ;print, max(dsr0.ke), max(dsr2.ke),max(dsr3.ke),max(dsr1.ke), format='("MAX KE=[",e9.2,",",e9.2,",",e9.2,",",e9.2,"]")'
 ;print, max(abs(dsr0.vpar)), max(abs(dsr2.vpar)),max(abs(dsr3.vpar)),max(abs(dsr1.vpar)), format='("MAX vpar=[",e9.2,",",e9.2,",",e9.2,",",e9.2,"]")'

!p.multi=0

END
