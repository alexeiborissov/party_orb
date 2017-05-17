;plotting jz contours overlaid with velovect of gas pressure, and magnetic tension and pressure
;@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1
;@xplot3dJT
loadct, 4, /silent
tvlct, r, g, b, /get


;f=11
out=1
od='./figs/newmov3/'
sz=256
sy=128
npts=4096
nom="Data/"

blah='../../laredata/twoloopsunstable/'
savdir='./saves/'

ztop=5
zbottom=-5
xtop=3
xbottom=-1
ytop=1
ybottom=-1

nframes=60
;gpmax=dblarr(nframes)
;lfmax=dblarr(nframes)
;bpmax=dblarr(nframes)
;btmax=dblarr(nframes)
;tfmax=dblarr(nframes)
;
;restore, filename='maxforces.sav', /verbose
;fmax=max([lfmax,gpmax])


FOR f=1,nframes-1 DO BEGIN
;FOR f=1,1 DO BEGIN
narrows=29
xnarrows=narrows
ynarrows=xnarrows
;print, 'getting data..'
fjnam=string(f, format='("modjstore",i3.3,".sav")')
blnam=string(f, format='("blstore",i3.3,".sav")')
;restore, filename=savdir+'kestore.sav', /verbose

ds=getdata(f,/grid, wkdir=blah)

jx_zslicenom=string(f, format='("jx_zslice",i3.3,".eps")')
jy_zslicenom=string(f, format='("jy_zslice",i3.3,".eps")')
;jz_zslicenom=string(f, format='("jz_zslice",i3.3,".eps")')
lorentz_jz_zslice=string(f, format='("lorentz_zslice",i3.3,".eps")')
lorentz_jz_yslice=string(f, format='("lorentz_yslice",i3.3,".eps")')
lorentz_modj_yslice=string(f, format='("jz_bt_bp_gp_yslice",i3.3,".eps")')
lorentz_modj_zslice=string(f, format='("jz_bt_bp_gp_zslice",i3.3,".eps")')

jxtitle=string(ds.time, format='("jx (at z-midplane, t=",f7.2,")")')
jytitle=string(ds.time, format='("jy (at z-midplane, t=",f7.2,")")')
jztitle=string(ds.time, format='("jz (at z-midplane, t=",f7.2,")")')
jztitle2=string(ds.time, format='("jz (at y-midplane, t=",f7.2,")")')
modjtitle=string(ds.time, format='("jz (at y-midplane, t=",f7.2,")")')
modjtitle2=string(ds.time, format='("jz (at z-midplane, t=",f7.2,")")')

grid=ds.grid
delx=grid.x[1]-grid.x[0]
dely=grid.y[1]-grid.y[0]
delz=grid.z[1]-grid.z[0]

x=grid.x
nx = n_elements(x)-1
dx = x(1)-x(0)
y=grid.y
ny = n_elements(y)-1
dy = y(1)-y(0)
z=grid.z
nz = n_elements(z)-1
dz = z(1)-z(0)

data = getdata(f,wkdir=blah,/bx) & bx=data.bx
data = getdata(f,wkdir=blah,/by) & by=data.by
data = getdata(f,wkdir=blah,/bz) & bz=data.bz

ztop=5
zbottom=-5
xtop=3
xbottom=-1
ytop=1
ybottom=-1

uz=where(abs(z-ztop) eq min(abs(z-ztop)))
lz=where(abs(z-zbottom) eq min(abs(z-zbottom)))
ux=where(abs(x-xtop) eq min(abs(x-xtop)))
lx=where(abs(x-xbottom) eq min(abs(x-xbottom)))
uy=where(abs(y-ytop) eq min(abs(y-ytop)))
ly=where(abs(y-ybottom) eq min(abs(y-ybottom)))

sz=where(abs(z-0) eq min(abs(z-0)))-lz
sy=where(abs(y-0) eq min(abs(y-0)))-ly

;;STOP
;; making the box smaller to be xy=[-0.75,0.75], z=[-1,2]
 nx=ux-lx;+1
 ny=uy-ly;+1
 nz=uz-lz;+1
 bgrid = dblarr(nx,ny,nz,3)
 bgrid[*,*,*,0] = bx(lx:ux-1,ly:uy-1,lz:uz-1)
 bgrid[*,*,*,1] = by(lx:ux-1,ly:uy-1,lz:uz-1)
 bgrid[*,*,*,2] = bz(lx:ux-1,ly:uy-1,lz:uz-1)

 print, 'destaggering'
 bgrid=destaggerB(bgrid)

 bx=reform(bgrid[*,*,*,0])
 by=reform(bgrid[*,*,*,1])
 bz=reform(bgrid[*,*,*,2])

;Move grid points so they are at the same locations as B in order to run 
;the null finding code
 xx = x(lx:ux)
 yy = y(ly:uy)
 zz = z(lz:uz)

 data = getdata(f,wkdir=blah,/jx) 
 jx=data.jx[lx:ux,ly:uy,lz:uz]
 data = getdata(f,wkdir=blah,/jy) 
 jy=data.jy[lx:ux,ly:uy,lz:uz]
 data = getdata(f,wkdir=blah,/jz) 
 jz=data.jz[lx:ux,ly:uy,lz:uz]
 jgrid = dblarr(nx+1,ny+1,nz+1,3)
 jgrid[*,*,*,0] = jx
 jgrid[*,*,*,1] = jy
 jgrid[*,*,*,2] = jz
 jgrid=destaggerj(jgrid)
 jz=jgrid[0:nx-1,0:ny-1,0:nz-1,0] 
 jy=jgrid[0:nx-1,0:ny-1,0:nz-1,1]
 jz=jgrid[0:nx-1,0:ny-1,0:nz-1,2]
 
 data = getdata(f,wkdir=blah,/pressure) 
 p=data.pressure[lx:ux-1,ly:uy-1,lz:uz-1]
 p=destaggerP(p)
 modj=sqrt(jx*jx+jy*jy+jz*jz)
 bsq=bgrid[*,*,*,0]*bgrid[*,*,*,0]+bgrid[*,*,*,1]*bgrid[*,*,*,1]+bgrid[*,*,*,2]*bgrid[*,*,*,2]

 BPx=-0.5*(bsq-shift(bsq,1,0,0))/(x[1]-x[0])
 BPy=-0.5*(bsq-shift(bsq,0,1,0))/(y[1]-y[0])
 BPz=-0.5*(bsq-shift(bsq,0,0,1))/(z[1]-z[0])
 BPx[0,*,*]=BPx[1,*,*]
 BPy[*,0,*]=BPy[*,1,*]
 BPz[*,*,0]=BPz[*,*,1]
 
 BTx=bx*(bx-shift(bx,1,0,0))/(x[1]-x[0])+by*(bx-shift(bx,0,1,0))/(y[1]-y[0])+bz*(bx-shift(bx,0,0,1))/(z[1]-z[0])
 BTy=bx*(by-shift(by,1,0,0))/(x[1]-x[0])+by*(by-shift(by,0,1,0))/(y[1]-y[0])+bz*(by-shift(by,0,0,1))/(z[1]-z[0])
 BTz=bx*(bz-shift(bz,1,0,0))/(x[1]-x[0])+by*(bz-shift(bz,0,1,0))/(y[1]-y[0])+bz*(bz-shift(bz,0,0,1))/(z[1]-z[0])
 BTx[0,*,*]=BTx[1,*,*]
 BTy[*,0,*]=BTy[*,1,*]
 BTz[*,*,0]=BTz[*,*,1]
 
 GPx=-0.5*(p-shift(p,1,0,0))/(x[1]-x[0])
 GPx[0,*,*]=GPx[1,*,*]
 GPy=-0.5*(p-shift(p,0,1,0))/(y[1]-y[0])
 GPy[*,0,*]=GPy[*,1,*]
 GPz=-0.5*(p-shift(p,0,0,1))/(z[1]-z[0])
 GPz[*,*,0]=GPz[*,*,1]
  
 TFx=BPx+BTx+GPx
 TFy=BPy+BTy+GPy
 Tfz=BPz+BTz+GPz

 bx=bx[0:nx-1,0:ny-1,0:nz-1]
 by=by[0:nx-1,0:ny-1,0:nz-1]
 bz=bz[0:nx-1,0:ny-1,0:nz-1]
 jx=jx[0:nx-1,0:ny-1,0:nz-1]
 jy=jy[0:nx-1,0:ny-1,0:nz-1]
 jz=jz[0:nx-1,0:ny-1,0:nz-1]
 lfx=(jy*bz-jz*by)
 lfy=(jz*bx-jx*bz)
 lfz=(jx*by-jy*bx)

 ;loadct, 4, /silent
 ;TVLCT, r,g,b, /get
 ;r=reverse(r)
 ;g=reverse(g)
 ;b=reverse(b)
 ;r(0)=0
 ;g(0)=0
 ;b(0)=0
 ;r(255)=255
 ;g(255)=255
 ;b(255)=255
 ;TVLCT, r, g, b
 posc=[0.08,0.88,0.95,0.93]
 posi=[0.12,0.13,0.97,0.78]

 posc2=[0.08,0.91,0.95,0.95]
 posi2=[0.2,0.1,0.95,0.85]

 aspect_ratio=1.5
 myxs=20
 IF NOT(OUT) THEN myxs=myxs*40
 myys=myxs/aspect_ratio


 IF (out) THEN BEGIN
  set_plot,'ps'
  !p.font=0
;  device, /close
  device, filename=od+lorentz_modj_zslice, encapsulated=1, /helvetica, /color, BITS_PER_PIXEL=8
  device, xsize=myxs, ysize=myys
  myth=4
  mycs=2 
  PRINT, 'output to: '+od+lorentz_modj_zslice
  pcol=0
  pbkgcol=255
 ENDIF ELSE BEGIN
  window, 10, ysize=myys, xsize=myxs
  myth=2
  mycs=1
  pcol=255
  pbkgcol=0
 ENDELSE
 !p.background=pbkgcol 
 loadct, 70, /silent
 ;JTblue_red
 TVLCT, r,g,b, /get
 r=reverse(r)
 g=reverse(g)
 b=reverse(b)
 r(0)=0
 g(0)=0
 b(0)=0
 r(255)=255
 g(255)=255
 b(255)=255
 TVLCT, r, g, b

 temp=reform(jz[*,*,sz])
 temp(where(temp gt 6.0))=6.0
 temp(where(temp lt -6.0))=-6.0

 contour,temp, x[lx:ux-1], y[ly:uy-1], /fill, levels=12*findgen(41)/40.-6.0,  xtitle='x', ytitle='y', pos=posi, /iso

 BTxtemp=reform(BTx[*,*,sz])
 BTytemp=reform(BTy[*,*,sz])
 BTx2=congrid(BTxtemp,xnarrows,ynarrows,/CENTER)
 BTy2=congrid(BTytemp,xnarrows,ynarrows,/CENTER)
 LOADCT, 39, /silent
 btmax=max(sqrt(btx2*btx2+bty2*bty2))
 velovect, BTx2, BTy2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(y[ly:uy-1],ynarrows,/CENTER), /overplot, col=190, length=btmax*1.5, thick=myth, pos=posi, /iso
 xyouts, -1, 1.1, 'B Tension',/data, color=190, charsize=mycs
 BPxtemp=reform(BPx[*,*,sz])
 BPytemp=reform(BPy[*,*,sz])
 BPx2=congrid(BPxtemp,xnarrows,ynarrows,/CENTER)
 BPy2=congrid(BPytemp,xnarrows,ynarrows,/CENTER)
 bpmax=max(sqrt(bpx2*bpx2+bpy2*bpy2))
 loadct, 8, /silent
 velovect, BPx2, BPy2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(y[ly:uy-1],ynarrows,/CENTER), /overplot, col=160, length=bpmax*1.5, thick=myth, pos=posi, /iso
 xyouts, 0, 1.1, 'B Pressure',/data, color=160, charsize=mycs
 lfx2=congrid(lfx[*,*,sz],xnarrows,ynarrows,/CENTER)
 lfy2=congrid(lfy[*,*,sz],xnarrows,ynarrows,/CENTER)
 lfmax=max(sqrt(lfx2*lfx2+lfy2*lfy2))
 ;velovect, lfx2, lfy2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(y[ly:uy-1],ynarrows,/CENTER), /overplot, col=pcol, length=lfmax, thick=myth
 ;xyouts, -0.5, 2.1, 'jxB Force',/data, color=pcol, charsize=mycs
 GPxtemp=reform(GPx[*,*,sz])
 GPytemp=reform(GPy[*,*,sz])
 GPx2=congrid(GPxtemp,xnarrows,ynarrows,/CENTER)
 GPy2=congrid(GPytemp,xnarrows,ynarrows,/CENTER)
 LOADCT, 2, /silent
 gpmax=max(sqrt(gpx2*gpx2+gpy2*gpy2))
 velovect, GPx2, GPy2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(y[ly:uy-1],ynarrows,/CENTER), /overplot, col=240, length=gpmax*1.5, thick=myth, pos=posi, /iso
 xyouts, 1, 1.1, 'Gas Pressure',/data, color=240, charsize=mycs
 TFxtemp=reform(TFx[*,*,sz])
 TFytemp=reform(TFy[*,*,sz])
 TFx2=congrid(TFxtemp,xnarrows,ynarrows,/CENTER)
 TFy2=congrid(TFytemp,xnarrows,ynarrows,/CENTER)
 ;LOADCT, 39, /silent
 ;tfmax=max(sqrt(tfx2*tfx2+tfy2*tfy2))
 ;velovect, TFx2, TFy2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(y[ly:uy-1],ynarrows,/CENTER), /overplot, col=240, length=tfmax, thick=myth
 ;xyouts, 2, 2.1, '-gradp+(Btens-Bp)',/data, color=240, charsize=mycharsize
 ;xyouts, -2.5, 10.2, 'Grad P',/data, color=180, charsize=mycharsize
 
 ;JTblue_red
 loadct, 70, /silent
 TVLCT, r,g,b, /get
 r=reverse(r)
 g=reverse(g)
 b=reverse(b)
 r(0)=0
 g(0)=0
 b(0)=0
 r(255)=255
 g(255)=255
 b(255)=255
 TVLCT, r, g, b 
 
 cgColorbar, format='(f5.1)', title=modjtitle2, range=[-6,6],$
 ncolors=251, bottom=3, pos=posc, color=pcol, charthick=myth-1, charsize=1, divisions=4
 
 aspect_ratio=0.5
 myxs=12
 IF NOT(OUT) THEN myxs=myxs*40
 myys=myxs/aspect_ratio
 narrows=22
 loadct, 70, /silent
 ;JTblue_red
 TVLCT, r,g,b, /get
 r=reverse(r)
 g=reverse(g)
 b=reverse(b)
 r(0)=0
 g(0)=0
 b(0)=0
 r(255)=255
 g(255)=255
 b(255)=255
 TVLCT, r, g, b
 IF (out) THEN BEGIN
  device, /close  
  device, filename=od+lorentz_modj_yslice, encapsulated=1, /helvetica, /color, BITS_PER_PIXEL=8
  device, xsize=myxs, ysize=myys
  myth=4  
  mycs=1.2
  PRINT, 'output to: '+od+lorentz_modj_yslice
  pcol=0
  pbkgcol=255
 ENDIF ELSE BEGIN
  window, 9, ysize=myys, xsize=myxs
  myth=2
  mycs=1
  pcol=255
  pbkgcol=0
 ENDELSE
 !p.background=pbkgcol 
 xnarrows=narrows+10
 znarrows=narrows
 ;temp=reform(modj[*,sy,*])
 temp=reform(jz[*,sy,*])
 temp(where(temp gt 6.0))=6.0
 temp(where(temp lt -6.0))=-6.0
 contour,temp, x[lx:ux-1], z[lz:uz-1], /fill, levels=12.*findgen(41)/40.-6.,  xtitle='x', ytitle='z', pos=posi2;, /iso

 ;narrows=narrows-10
 BTxtemp=reform(BTx[*,sy,*])
 BTztemp=reform(BTz[*,sy,*])
 BTx2=congrid(BTxtemp,xnarrows,znarrows,/CENTER)
 BTz2=congrid(BTztemp,xnarrows,znarrows,/CENTER)
 LOADCT, 39, /silent
 btmax=max(sqrt(btx2*btx2+btz2*btz2))
 velovect, BTx2, BTz2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(z[lz:uz-1],znarrows,/CENTER), /overplot, col=190, length=btmax*1.5, thick=myth, pos=posi2
 xyouts, -1, 5.1, 'B Tens',/data, color=190, charsize=mycs
 BPxtemp=reform(BPx[*,sy,*])
 BPztemp=reform(BPz[*,sy,*])
 BPx2=congrid(BPxtemp,xnarrows,znarrows,/CENTER)
 BPz2=congrid(BPztemp,xnarrows,znarrows,/CENTER)
 bpmax=max(sqrt(bpx2*bpx2+bpz2*bpz2))
 loadct, 8, /silent
 velovect, BPx2, BPz2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(z[lz:uz-1],znarrows,/CENTER), /overplot, col=160, length=bpmax*1.5, thick=myth, pos=posi2
 xyouts, 0, 5.1, 'B Pres',/data, color=160, charsize=mycs
 
 lfx2=congrid(reform(lfx[*,sy,*]),xnarrows,znarrows,/CENTER)
 lfz2=congrid(reform(lfz[*,sy,*]),xnarrows,znarrows,/CENTER)
 lfmax=max(sqrt(lfx2*lfx2+lfz2*lfz2))
 ;velovect, lfx2, lfz2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(z[lz:uz-1],znarrows,/CENTER), /overplot, col=pcol, length=lfmax, thick=myth, pos=posi2
 ;xyouts, -0.5, 5.1, 'jxB Force',/data, color=pcol, charsize=mycs
 
 GPxtemp=reform(GPx[*,sy,*])
 GPztemp=reform(GPz[*,sy,*])
 GPx2=congrid(GPxtemp,xnarrows,znarrows,/CENTER)
 GPz2=congrid(GPztemp,xnarrows,znarrows,/CENTER)
 LOADCT, 2, /silent
 gpmax=max(sqrt(gpx2*gpx2+gpz2*gpz2))
 velovect, GPx2, GPz2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(z[lz:uz-1],znarrows,/CENTER), /overplot, col=240, length=gpmax*1.5, thick=myth, pos=posi2
 xyouts, 1, 5.1, 'Gas Pres',/data, color=240, charsize=mycs
 TFxtemp=reform(TFx[*,sy,*])
 TFztemp=reform(TFz[*,sy,*])
 TFx2=congrid(TFxtemp,xnarrows,znarrows,/CENTER)
 TFz2=congrid(TFztemp,xnarrows,znarrows,/CENTER)
 ;LOADCT, 39, /silent
 ;tfmax=max(sqrt(tfx2*tfx2+tfz2*tfz2))
 ;velovect, TFx2, TFz2, congrid(x[lx:ux-1],xnarrows,/CENTER), congrid(z,znarrows,/CENTER), /overplot, col=240, length=tfmax, thick=myth, pos=posi2
 ;xyouts, 2, 2.2, '-gradp+(Btens-Bp)',/data, color=240, charsize=mycharsize
 ;xyouts, -2.5, 10.2, 'Grad P',/data, color=180, charsize=mycharsize
; 
; 
 loadct, 70, /silent
 ;JTblue_red
 TVLCT, r,g,b, /get
 r=reverse(r)
 g=reverse(g)
 b=reverse(b)
 r(0)=0
 g(0)=0
 b(0)=0
 r(255)=255
 g(255)=255
 b(255)=255
 TVLCT, r, g, b 
 cgColorbar,range=[-6,6], format='(f5.1)', title=modjtitle, $
 ncolors=251, bottom=3, pos=posc2, color=pcol, charthick=myth-1, charsize=1, divisions=4

  
 IF (out) THEN BEGIN
  device, /close
  set_plot, 'x'
  !p.font=-1
 ENDIF

ENDFOR
 ;save, gpmax, lfmax, tfmax, bpmax, btmax, filename='maxforces.sav', /compress, /verbose
END
