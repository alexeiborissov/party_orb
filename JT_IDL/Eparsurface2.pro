;interested in visualising Epar as a surface, but on a subgrid level.
; specifically interested in Eparallel between x=[99,100], y=[99.5,100.5], z=[74.2,75.2]

restore=0
;electric and magnetic field stored? known?
BF='../../../data/VAR378_B.xdr'
EF='../../../data/VAR378_E.xdr'
IF restore THEN BEGIN
print, 'restoring data from Bourdin snapshots'
 RESTORE, filename=BF, /verbose
ENDIF

;xr=where((x ge 99.5e6) and (x le 101e6))
;yr=where((y ge 99.5e6) and (y le 101e6))
;zr=where((z ge 73.2e6) and (z le 76.2e6))
xr=where((x ge 99e6) and (x le 100e6))
yr=where((y ge 99.5e6) and (y le 100.5e6))
zr=where((z ge 73.5e6) and (z le 76.2e6))

xrr=where((x ge 98.5e6) and (x le 100.5e6))
yrr=where((y ge 99.3e6) and (y le 100.8e6))
zrr=where((z ge 73.2e6) and (z le 75.7e6))
nxr=n_elements(xrr)
nyr=n_elements(yrr)
nzr=n_elements(zrr)

print, 'xregion:', xr
print, 'yregion:', yr
print, 'zregion:', zr
nx=n_elements(xr)
ny=n_elements(yr)
nz=n_elements(zr)

sx=x[xr[0]:xr[nx-1]]	; define new x, y, z grids.
sy=y[yr[0]:yr[ny-1]]
sz=z[zr[0]:zr[nz-1]]

IF restore THEN BEGIN

 sgB=B[xr[0]:xr[nx-1],yr[0]:yr[ny-1],zr[0]:zr[nz-1],*]
 sgBr=B[xrr[0]:xrr[nxr-1],yrr[0]:yrr[nyr-1],zrr[0]:zrr[nzr-1],*]
 B=0
 RESTORE, filename=EF, /verbose
 sgE=E(xr[0]:xr[nx-1],yr[0]:yr[ny-1],zr[0]:zr[nz-1],*)
 sgEr=E[xrr[0]:xrr[nxr-1],yrr[0]:yrr[nyr-1],zrr[0]:zrr[nzr-1],*]
 E=0
ENDIF

nnx=20
nny=20
nnz=10

cx=congrid(sx,nnx, /interp,/minus_one)
cy=congrid(sy,nny, /interp,/minus_one)
cz=congrid(sz,nnz, /interp,/minus_one)

;cgBx=smooth(cgBx,5,/edge_truncate)

cgBx=congrid(sgB[*,*,*,0],nnx,nny,nnz,/interp,/minus_one)
cgBy=congrid(sgB[*,*,*,1],nnx,nny,nnz,/interp,/minus_one)
cgBz=congrid(sgB[*,*,*,2],nnx,nny,nnz,/interp,/minus_one)
;cgBx=smooth(cgBx,5,/edge_truncate)
;cgBy=smooth(cgBy,5,/edge_truncate)
;cgBz=smooth(cgBz,5,/edge_truncate)


cgEx=congrid(sgE[*,*,*,0],nnx,nny,nnz,/interp,/minus_one)
cgEy=congrid(sgE[*,*,*,1],nnx,nny,nnz,/interp,/minus_one)
cgEz=congrid(sgE[*,*,*,2],nnx,nny,nnz,/interp,/minus_one)
;cgEx=smooth(cgEx,5,/edge_truncate)
;cgEy=smooth(cgEy,5,/edge_truncate)
;cgEz=smooth(cgEz,5,/edge_truncate)

cEpar=cgEx*cgBx+cgEy*cgBy+cgEz*cgBz
;cEpar is now a 100x100x50 grid surrounding the area we want.

cxr=where((cx ge 99e6) and (cx le 100e6))
cyr=where((cy ge 99.5e6) and (cy le 100.5e6))
czr=where((cz ge 74.2e6) and (cz le 75.2e6))
cnx=n_elements(cxr)
cny=n_elements(cyr)
cnz=n_elements(czr)

;what about using path_xy in contour?

xplot3d, [99e6,99e6]/lscl, [99.5e6,99.5e6]/lscl, [74.2e6,74.2e6]/lscl, $
xr=[99e6,100e6]/lscl, yr=[99.5e6,100.5e6]/lscl, zr=[74.2e6,75.2e6]/lscl,xtitle=xt, ytitle=yt, ztitle=zt

;FOR i=0,cnz-1 DO BEGIN
; contour, cEpar[cxr[0]:cxr[cnx-1],cyr[0]:cyr[cny-1],czr[i]], cx(cxr[0]:cxr[cnx-1]), cy(cyr[0]:cyr[cny-1]), level=0, xstyle=1, ystyle=1, path_info=info, path_xy=xy, /path_data_coords
 ;help, xy
 ;wait, 0.5
; xyz=DBLARR(3,info.n-2)
; xyz[0,*]=xy[0,0:info.n-3]
; xyz[1,*]=xy[1,0:info.n-3]
; xyz[2,*]=cz[czr[i]]
 
 ;xplot3d, xyz[0,*]/lscl, xyz[1,*]/lscl, xyz[2,*]/lscl, $
;xr=[99e6,100e6]/lscl, yr=[99.5e6,100.5e6]/lscl, zr=[74.2e6,75.2e6]/lscl,xtitle=xt, ytitle=yt, ztitle=zt, /overplot
;ENDFOR

nom="Data/"
particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], /op, zscl=lscl, /mm, lscl=lscl, /rel, /symb, /bsymb, frac=0.01



 oOrb3 = OBJ_NEW('orb', COLOR=[0, 0 ,255]) 	;blue orb
 oOrb4 = OBJ_NEW('orb', COLOR=[255, 0 ,0]) 	;red orb
 oOrb5 = OBJ_NEW('orb', COLOR=[128, 128 ,240]) 	;blue orb
 oOrb6 = OBJ_NEW('orb', COLOR=[240, 128,128]) 	;red orb
 rmult=1.0
 xyrat=1.
 frac=0.02
 xylen=1E6
 zlen=1E6
 oOrb3->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 oOrb4->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 oSymbol3 = OBJ_NEW('IDLgrSymbol', oOrb3) ;lets visualise the grid
 oSymbol4 = OBJ_NEW('IDLgrSymbol', oOrb4) ;color gridpoints by epar +ve/-ve
 frac=0.01
 oOrb5->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 oOrb6->Scale, frac*xyrat, frac*xyrat, zlen/xylen*xyrat*frac
 oSymbol5 = OBJ_NEW('IDLgrSymbol', oOrb5) ;lets visualise the grid
 oSymbol6 = OBJ_NEW('IDLgrSymbol', oOrb6) ;color gridpoints by epar +ve/-ve


sgEpar=reform(sgEr[*,*,*,0]*sgBr[*,*,*,0]+sgEr[*,*,*,1]*sgBr[*,*,*,1]+sgEr[*,*,*,2]*sgBr[*,*,*,2])

;temp=dblarr(nx,ny)
;temp[*,ny,0]=y(yr)
;temp[nx,*]=x(xr)

FOR ix=0,nxr-1 DO BEGIN
 FOR iy=0,nyr-1 DO BEGIN
  FOR iz=0,nzr-1 DO BEGIN
   IF sgEpar[ix,iy,iz] gt 0.0D0 THEN BEGIN
    XPLOT3D, [x[xrr[ix]],x[xrr[ix]]]/lscl, [y[yrr[iy]],y[yrr[iy]]]/lscl, [z[zrr[iz]],z[zrr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol4, THICK=2, /OVERPLOT; +ve Epar=red
   ENDIF ELSE BEGIN
    XPLOT3D, [x[xrr[ix]],x[xrr[ix]]]/lscl, [y[yrr[iy]],y[yrr[iy]]]/lscl, [z[zrr[iz]],z[zrr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol3, THICK=2, /OVERPLOT; -ve Epar=blue
   ENDELSE
  ENDFOR
 ENDFOR
ENDFOR
FOR ix=0,cnx-1 DO BEGIN
 FOR iy=0,cny-1 DO BEGIN
  FOR iz=0,cnz-1 DO BEGIN
   IF cEpar[cxr[ix],cyr[iy],czr[iz]] gt 0.0D0 THEN BEGIN
    XPLOT3D, [cx[cxr[ix]],cx[cxr[ix]]]/lscl, [cy[cyr[iy]],cy[cyr[iy]]]/lscl, [cz[czr[iz]],cz[czr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol6, THICK=2, /OVERPLOT; +ve Epar=red
   ENDIF ELSE BEGIN
    XPLOT3D, [cx[cxr[ix]],cx[cxr[ix]]]/lscl, [cy[cyr[iy]],cy[cyr[iy]]]/lscl, [cz[czr[iz]],cz[czr[iz]]]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymbol5, THICK=2, /OVERPLOT; -ve Epar=blue
   ENDELSE
  ENDFOR
 ENDFOR
ENDFOR



;cxmin=where(min(abs(cx-99e6)) eq abs(cx-99e6))
;cxmax=where(min(abs(cx-100e6)) eq abs(cx-100e6))
;cymin=where(min(abs(cy-99e6)) eq abs(cy-99e6))
;cymax=where(min(abs(cy-100e6)) eq abs(cy-100e6))
;czmin=where(min(abs(cz-74.2e6)) eq abs(cz-74.2e6))
;czmax=where(min(abs(cz-75.2e6)) eq abs(cz-75.2e6))

;xby0z0loc=where(min(abs(cEpar[*,cymin,czmin])) eq abs(cEpar[*,cymin,czmin]))
;xby0z0loc=where(min(abs(cEpar[*,cymin,czmin])) eq abs(cEpar[*,cymin,czmin]))


;xplot3d, [cxmin,cxmin]/1e6, [cymin,cymin]/1e6, [czmin,czmax]/1e6, xr=[99,100], yr=[99.5,100.5], zr=[74.2,75.2], linestyle=2 
;xplot3d, [cxmax,cxmax]/1e6, [cymin,cymin]/1e6, [czmin,czmax]/1e6, /oplot, linestyle=2

;
;particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], zscl=1e6, /mm, lscl=1e6, /rel, /symb, /bsymb, frac=0.01, /op
;lscl=1e6

;xscale=(cxr[cnx-1]-cxr[0])/DOUBLE(cnx)/lscl
;yscale=(cyr[cny-1]-cyr[0])/DOUBLE(cny)/lscl
;zscale=(czr[cnz-1]-czr[0])/DOUBLE(cnz)/lscl
;xscale=(max(cx[cxr])-min(cx[cxr]))/DOUBLE(cnx)/lscl
;yscale=(max(cy[cyr])-min(cy[cyr]))/DOUBLE(cny)/lscl
;zscale=(max(cz[czr])-min(cz[czr]))/DOUBLE(cnz)/lscl



xt='x(Mm)'
yt='y(Mm)'
zt='z(Mm)'

;minmax(cEpar)=[-4.12e-6,2.9e-6] 

;theModel = Obj_New('IDLgrModel')
;theModel -> Scale, xscale,yscale,zscale	    	    	
;theModel -> Add, mk_iso(cEpar[cxr[0]:cxr[cnx-1],cyr[0]:cyr[cny-1],czr[0]:czr[cnz-1]],0.0,scol=[0,0,255], /low)	
;;theModel -> Add, mk_iso(cEpar[cxr[0]:cxr[cnx-1],cyr[0]:cyr[cny-1],czr[0]:czr[cnz-1]],0.0,scol=[0,0,255],/low, maxval=3.5e-7);;

;xobjview, theModel;
;;stop
;test=obj_new('IDLgrSymbol', theModel)
;if obj_valid(test) THEN test->setproperty, size=1
;if obj_valid(test) THEN test->SetProperty, COLOR = [0, 0, 255]

;xplot3d, [99e6,99e6]/lscl, [99.5e6,99.5e6]/lscl, [74.2e6,74.2e6]/lscl, symbol=test, $
;xr=[99e6,100e6]/lscl, yr=[99.5e6,100.5e6]/lscl, zr=[74.2e6,75.2e6]/lscl,xtitle=xt, ytitle=yt, ztitle=zt




END
