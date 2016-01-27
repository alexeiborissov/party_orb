;can we visualise the sign of epar at each position in the particle path?

party=10

ax=-50
az=10
;ax=0
;az=0
filewidth=800
lscl=1e6
xt='x(Mm)'
yt='y(Mm)'
zt='z(Mm)'


;e&p 10
xl=55e6
xu=165e6
yl=100e6
yu=210e6
zl=0.5e6
zu=90e6

;electron 25
;xl=141e6
;xu=151e6
;yl=146e6
;yu=156e6
;zl=74e6
;zu=90e6

;proton 16
;xl=70e6
;xu=100e6
;yl=180e6
;yu=210e6
;zl=73e6
;zu=77e6

;proton 17
;xl=124.9e6
;xu=125.1e6
;yl=174.9e6
;yu=175.1e6
;zl=74.9e6
;zu=75.1e6

;proton 14
;xl=124.5e6
;xu=126.5e6
;yl=124.5e6
;yu=126.5e6
;zl=79.5e6
;zu=80.0e6


xplot3d, [xl,xl]/lscl, [yl,yl]/lscl, [zl,zl]/lscl, ax=ax, az=az, filewidth=filewidth,$
xr=[xl,xu]/lscl, yr=[yl,yu]/lscl, zr=[zl,zu]/lscl,xtitle=xt, ytitle=yt, ztitle=zt

nom="Data/"

particletrack, party, floc=nom, xyzt, lcol=[0,0,0], /op, zscl=lscl, /mm, lscl=lscl, /rel, /symb, /bsymb, frac=1.0, myxrange=[xl,xu], myyrange=[yl,yu], myzrange=[zl,zu]

ds=getrdata(party,/ke,/rel,/fields, wkdir=nom+'DataR/')

nom="../quicksurveyP/Data/"

particletrack, party, floc=nom, xyzt, lcol=[0,0,255], /op, zscl=lscl, /mm, lscl=lscl, /rel, /symb, /bsymb, frac=1.0, myxrange=[xl,xu], myyrange=[yl,yu], myzrange=[zl,zu]

dsp=getrdata(party,/ke,/rel,/fields, wkdir=nom+'DataR/')


;STOP

restore=1
;electric and magnetic field stored? known?
BF='../../../../2014/bourdin/data/VAR378_B.xdr'
EF='../../../../2014/bourdin/data/VAR378_E.xdr'
EPF='../../../../2014/bourdin/data/VAR378_E_PARALLEL.xdr'
IF restore THEN BEGIN
print, 'restoring data from Bourdin snapshots'
 RESTORE, filename=BF, /verbose
ENDIF


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
stop
 ch = 0.01*dx 	    	    ;step size
 chmin = 0.001*dx	    	    ;min step size
 chmax = 0.5*dx	    	    ; max step size
 cepsilon = 1.0d-5  	    ; tolerance
 cmxline = 100000    	    ; max nsteps
 csz_stp = size([ds.x(0),ds.y(0),ds.z(0)])

 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]

jj=0
  ;for jj=0,n_elements(ds.t)-1 do begin
  linef = tracefield_3D([ds.x(jj),ds.y(jj),ds.z(jj)],b,x,y,z,ch,chmin,chmax,cepsilon,cmxline)
  ii = where(linef[0,*] gt -8000,nii)
  ;if (ii[0] ne -1) then linef=linef[*,ii]
  ;jjf = where(linef[0,*] ge xl-1e6 and linef[0,*] le xu+1e6 and linef[1,*] ge yl-1e6 and linef[1,*] le yu+1e6 and linef[2,*] ge zl-1e6 and linef[2,*] le zu+1e6)
  ;if (jjf[0] ne -1) then xplot3d, linef(0,jjf)/lscl, linef(1,jjf)/lscl, linef(2,jjf)/lscl, COLOR=steelblue, /OVERPLOT
  if (ii[0] ne -1) then xplot3d, linef(0,ii)/lscl, linef(1,ii)/lscl, linef(2,ii)/lscl, COLOR=steelblue, /OVERPLOT
  ;trace field line backwards
  lineb = tracefield_3D([ds.x(jj),ds.y(jj),ds.z(jj)],b,x,y,z,-ch,chmin,chmax,cepsilon,cmxline)
  ii = where(lineb[0,*] gt -8000,nii)
  ;if (ii[0] ne -1) then lineb=lineb[*,ii]
  ;jjb = where(lineb[0,*] ge xl-1e6 and lineb[0,*] le xu+1e6 and lineb[1,*] ge yl-1e6 and lineb[1,*] le yu+1e6 and lineb[2,*] ge zl-1e6 and lineb[2,*] le zu+1e6)
  ;if (jjb[0] ne -1) then xplot3d, lineb(0,jjb)/lscl, lineb(1,jjb)/lscl, lineb(2,jjb)/lscl, COLOR=steelblue, /OVERPLOT
  if (ii[0] ne -1) then xplot3d, lineb(0,ii)/lscl, lineb(1,ii)/lscl, lineb(2,ii)/lscl, COLOR=steelblue, /OVERPLOT
  ;endfor
 
  STOP
  


xr=where((x ge xl) and (x le xu))
yr=where((y ge yl) and (y le yu))
zr=where((z ge zl) and (z le zu))

nx=n_elements(xr)
ny=n_elements(yr)
nz=n_elements(zr)

sx=x[xr[0]-1:xr[nx-1]+1]	; define new x, y, z grids.
sy=y[yr[0]-1:yr[ny-1]+1]
sz=z[zr[0]-1:zr[nz-1]+1]

IF restore THEN BEGIN
 ;PRINT, '-B[0,0,0]-
 ;PRINT, 'IDL:', b[1,1,1,*], format='(A," ",E17.10," ",E17.10," ",E17.10)'
 sgB=B[xr[0]-1:xr[nx-1]+1,yr[0]-1:yr[ny-1]+1,zr[0]-1:zr[nz-1]+1,*]
 ;sgBr=B[xrr[0]:xrr[nxr-1],yrr[0]:yrr[nyr-1],zrr[0]:zrr[nzr-1],*]
 B=0
 RESTORE, filename=EPF, /verbose
 ;PRINT, '-E[0,0,0]-
 ;PRINT, 'IDL:', E[1,1,1,*], format='(A," ",E17.10," ",E17.10," ",E17.10)'
; sgE=E(xr[0]-1:xr[nx-1]+1,yr[0]-1:yr[ny-1]+1,zr[0]-1:zr[nz-1]+1,*)
 ;sgEr=E[xrr[0]:xrr[nxr-1],yrr[0]:yrr[nyr-1],zrr[0]:zrr[nzr-1],*]
 ;E=0
ENDIF
;STOP






dx=x[1]-x[0]
dy=y[1]-y[0]
dz=z[200]-z[199]

inipos=[ds.x[0],ds.y[0],ds.z[0]]

;npx=20
;npy=20
;npz=20

;newnx=nx*(npx)-nx+1
;newny=ny*(npy)-ny+1
;newnz=nz*(npz)-nz+1

;tgx=x[xr[0]]+dindgen(nx*(npx)-nx+1)/float(npx-1)*dx
;tgy=y[yr[0]]+dindgen(ny*(npy)-ny+1)/float(npy-1)*dx
;tgz=z[zr[0]]+dindgen(nz*(npz)-nz+1)/float(npz-1)*dz;

;Epargrid=dblarr(nx*npx-(nx-1),ny*npy-(ny-1),nz*npz-(nz-1));

;FOR iz=1,nz DO BEGIN
; FOR iy=1,ny DO BEGIN
;  FOR ix=1,nx DO BEGIN
;    tgE=create_smbgrid(sgE,ix,iy,iz,npx,npy,npz,1)
;    ;tgEy=create_smbgrid(sgE[*,*,*,1],ix,iy,iz,npx,npy,npz,1)
;    ;tgEz=create_smbgrid(sgE[*,*,*,2],ix,iy,iz,npx,npy,npz,1)
;    tgB=create_smbgrid(sgB,ix,iy,iz,npx,npy,npz,1)
;    ;tgBy=create_smbgrid(sgB[*,*,*,2],ix,iy,iz,npx,npy,npz,1)
;    ;tgBz=create_smbgrid(sgB[*,*,*,3],ix,iy,iz,npx,npy,npz,1)
;    Epargrid[(ix-1)*(npx-1):(ix-1)*(npx-1)+npx-1,(iy-1)*(npy-1):(iy-1)*(npy-1)+npy-1,(iz-1)*(npz-1):(iz-1)*(npz-1)+npz-1]=tgE[*,*,*,0]*tgB[*,*,*,0]+tgE[*,*,*,1]*tgB[*,*,*,1]+tgE[*,*,*,2]*tgB[*,*,*,2]
;   ENDFOR
;  ENDFOR   
; ENDFOR   
    
    
;nnx=20
;nny=20
;nnz=10

;cx=congrid(sx,nnx, /interp,/minus_one)
;cy=congrid(sy,nny, /interp,/minus_one)
;cz=congrid(sz,nnz, /interp,/minus_one)

;cgBx=smooth(cgBx,5,/edge_truncate)

;cgBx=congrid(sgB[*,*,*,0],nnx,nny,nnz,/interp,/minus_one)
;cgBy=congrid(sgB[*,*,*,1],nnx,nny,nnz,/interp,/minus_one)
;cgBz=congrid(sgB[*,*,*,2],nnx,nny,nnz,/interp,/minus_one)
;cgBx=smooth(cgBx,5,/edge_truncate)
;cgBy=smooth(cgBy,5,/edge_truncate)
;cgBz=smooth(cgBz,5,/edge_truncate)


;cgEx=congrid(sgE[*,*,*,0],nnx,nny,nnz,/interp,/minus_one)
;cgEy=congrid(sgE[*,*,*,1],nnx,nny,nnz,/interp,/minus_one)
;cgEz=congrid(sgE[*,*,*,2],nnx,nny,nnz,/interp,/minus_one)
;cgEx=smooth(cgEx,5,/edge_truncate)
;cgEy=smooth(cgEy,5,/edge_truncate)
;cgEz=smooth(cgEz,5,/edge_truncate)

;cEpar=cgEx*cgBx+cgEy*cgBy+cgEz*cgBz
;;cEpar is now a 100x100x50 grid surrounding the area we want.

;tgxr=where((tgx ge xl) and (tgx le xu))
;tgyr=where((tgy ge yl) and (tgy le yu))
;tgzr=where((tgz ge zl) and (tgz le zu));

;tgnx=n_elements(tgxr)
;tgny=n_elements(tgyr)
;tgnz=n_elements(tgzr)

;what about using path_xy in contour?
;mycontour=0
myiso=1
;IF mycontour THEN BEGIN
 ;FOR i=0,tgnz-1 DO BEGIN
 ; contour, Epargrid[tgxr[0]:tgxr[tgnx-1],tgyr[0]:tgyr[tgny-1],tgzr[i]], tgx(tgxr[0]:tgxr[tgnx-1]), tgy(tgyr[0]:tgyr[tgny-1]), level=0, xstyle=1, ystyle=1, path_info=info, path_xy=xy, /path_data_coords
 ; IF n_elements(xy) gt 1 THEN BEGIN
 ;  FOR j=0,n_elements(info)-1 DO BEGIN
 ;   IF (info[j].n gt 4) THEN BEGIN
 ;;    xyz=DBLARR(3,info[j].n-3)
 ;    xyz[0,*]=xy[0,1:info[j].n-3]
 ;    xyz[1,*]=xy[1,1:info[j].n-3]
 ;    xyz[2,*]=tgz[tgzr[i]]
 ;
 ;    xplot3d, xyz[0,*]/lscl, xyz[1,*]/lscl, xyz[2,*]/lscl, $
 ;    xr=[141e6,151e6]/lscl, yr=[146e6,156e6]/lscl, zr=[74e6,90e6]/lscl,xtitle=xt, ytitle=yt, ztitle=zt, /overplot
;    ENDIF
;   ENDFOR
;  ENDIF
; ENDFOR
;ENDIF
IF myiso THEN BEGIN
; xscale=(max(tgx[tgxr])-min(tgx[tgxr]))/DOUBLE(tgnx)/lscl
; yscale=(max(tgy[tgyr])-min(tgy[tgyr]))/DOUBLE(tgny)/lscl
; zscale=(max(tgz[tgzr])-min(tgz[tgzr]))/DOUBLE(tgnz)/lscl

 xscale=(max(x[xr])-min(x[xr]))/DOUBLE(nx)/lscl
 yscale=(max(y[yr])-min(y[yr]))/DOUBLE(ny)/lscl
 zscale=(max(z[zr])-min(z[zr]))/DOUBLE(nz)/lscl

; theModel = Obj_New('IDLgrModel')
; theModel -> Scale, xscale,yscale,zscale	    	    	
; theModel -> Add, mk_iso(Epargrid[tgxr[0]:tgxr[tgnx-1],tgyr[0]:tgyr[tgny-1],tgzr[0]:tgzr[tgnz-1]],0.0,scol=[255,0,0], /low,transparency=0.5,minval=-2e-6)	
; test=obj_new('IDLgrSymbol', theModel)
; if obj_valid(test) THEN test->setproperty, size=1
; if obj_valid(test) THEN test->SetProperty, COLOR = [255, 0, 0]
; xplot3d, [xl,xl]/lscl, [yl,yl]/lscl, [zl,zl]/lscl, symbol=test, /overplot;, $
 
theModel = Obj_New('IDLgrModel')
theModel -> Scale, xscale,yscale,zscale	    	    	
theModel -> Add, mk_iso(E_PARALLEL(xr[0]-1:xr[nx-1]+1,yr[0]-1:yr[ny-1]+1,60:zr[nz-1]+1),1e-3,scol=[255,0,0], /low,transparency=0.5,minval=1e-6)	
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=1
 if obj_valid(test) THEN test->SetProperty, COLOR = [255, 0, 0]
 xplot3d, [xl,xl]/lscl, [yl,yl]/lscl, [z(60),z(60)]/lscl, symbol=test, /overplot;, $
 
 
ENDIF


;;;---this stuff was used to test the calculation of E||----;;
; oOrb3 = OBJ_NEW('orb', COLOR=[0, 0 ,255]) 	;blue orb
; oOrb4 = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 	;red orb
; oOrb5 = OBJ_NEW('orb', COLOR=[128, 128 ,240]) 	;blue orb
; oOrb6 = OBJ_NEW('orb', COLOR=[240, 128,128]) 	;red orb
; rmult=1.0
; xyrat=1.
; frac=0.02
; xylen=1E6
; zlen=1E6
; oOrb3->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
; oOrb4->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
; oSymbol3 = OBJ_NEW('IDLgrSymbol', oOrb3) ;lets visualise the grid
; oSymbol4 = OBJ_NEW('IDLgrSymbol', oOrb4) ;color gridpoints by epar +ve/-ve
; frac=0.01
; oOrb5->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
; oOrb6->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
; oSymbol5 = OBJ_NEW('IDLgrSymbol', oOrb5) ;lets visualise the grid
; oSymbol6 = OBJ_NEW('IDLgrSymbol', oOrb6) ;color gridpoints by epar +ve/-ve


;sgEpar=reform(sgEr[*,*,*,0]*sgBr[*,*,*,0]+sgEr[*,*,*,1]*sgBr[*,*,*,1]+sgEr[*,*,*,2]*sgBr[*,*,*,2])

;temp=dblarr(nx,ny)
;temp[*,ny,0]=y(yr)
;temp[nx,*]=x(xr)

;FOR ix=0,nxr-1 DO BEGIN
; FOR iy=0,nyr-1 DO BEGIN
;  FOR iz=0,nzr-1 DO BEGIN
;   IF sgEpar[ix,iy,iz] gt 0.0D0 THEN BEGIN
;    XPLOT3D, [x[xrr[ix]],x[xrr[ix]]]/lscl, [y[yrr[iy]],y[yrr[iy]]]/lscl, [z[zrr[iz]],z[zrr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol4, THICK=2, /OVERPLOT; +ve Epar=red
;   ENDIF ELSE BEGIN
;    XPLOT3D, [x[xrr[ix]],x[xrr[ix]]]/lscl, [y[yrr[iy]],y[yrr[iy]]]/lscl, [z[zrr[iz]],z[zrr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol3, THICK=2, /OVERPLOT; -ve Epar=blue
;   ENDELSE
;  ENDFOR
; ENDFOR
;ENDFOR
;FOR ix=0,cnx-1 DO BEGIN
; FOR iy=0,cny-1 DO BEGIN
;  FOR iz=0,cnz-1 DO BEGIN
;   IF cEpar[cxr[ix],cyr[iy],czr[iz]] gt 0.0D0 THEN BEGIN
;    XPLOT3D, [cx[cxr[ix]],cx[cxr[ix]]]/lscl, [cy[cyr[iy]],cy[cyr[iy]]]/lscl, [cz[czr[iz]],cz[czr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol6, THICK=2, /OVERPLOT; +ve Epar=red
;   ENDIF ELSE BEGIN
;    XPLOT3D, [cx[cxr[ix]],cx[cxr[ix]]]/lscl, [cy[cyr[iy]],cy[cyr[iy]]]/lscl, [cz[czr[iz]],cz[czr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol5, THICK=2, /OVERPLOT; -ve Epar=blue
;   ENDELSE
;  ENDFOR
; ENDFOR
;ENDFOR
;FOR ix=0,newnx-1 DO BEGIN
; FOR iy=0,newny-1 DO BEGIN
;  FOR iz=0,newnz-1 DO BEGIN
;   IF Epargrid[ix,iy,iz] gt 0.0D0 THEN BEGIN
 ;   XPLOT3D, [tgx[ix],tgx[ix]]/lscl, [tgy[iy],tgy[iy]]/lscl, [tgz[iz],tgz[iz]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol6, THICK=2, /OVERPLOT; +ve Epar=red
 ;  ENDIF ELSE BEGIN
 ;   XPLOT3D, [tgx[ix],tgx[ix]]/lscl, [tgy[iy],tgy[iy]]/lscl, [tgz[iz],tgz[iz]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol5, THICK=2, /OVERPLOT; -ve Epar=blue
 ;  ENDELSE
 ; ENDFOR
 ;ENDFOR
;ENDFOR


END
