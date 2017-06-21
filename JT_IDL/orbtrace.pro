;program to trace particles, field lines and current iso-surfaces in Alan's 3D MHD experiments.
;@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1
;@xplot3dJT
loadct, 4, /silent
tvlct, r, g, b, /get

f=02
filename='test.png'
ax=-60
az=45
filewidth=1000
transparent=0.5

npts=16*16
nom="Data/"

arbscl=1e6
lscl=1e6

mythick=4
 steelblue=[188,210,238]
 brightgreen=[0,255,0]
 brightpink=[170,0,220]
 purple=[153,51,255]
 grey=[170,170,170]
 poobrown=[153,76,0]
 olive=[128,128,0]
 salmon=[255,160,122]
 lightsteelblue=[176,196,222]
 lavender=[230,230,250]
 rosybrown=[188,143,143]
 burlywood=[222,184,135]
 tann=[210,180,140]
 newpurple=[128,0,128]
 newpink=[255,153,255]
;col1=lightsteelblue
;col2=tann
col1=newpurple
col2=newpink

IF keyword_set(filename) and keyword_set(transparent) THEN BEGIN
 type=strsplit(filename,".",/extract)
 newfilenames=[type[0]+'_a.'+type[1],type[0]+'_b.'+type[1]]
ENDIF

;print, 'getting data..'
;blah='../../../laredata/oneloopunstable/'
savdir='./saves/'
;fjnam=string(f, format='("modjstore",i3.3,".sav")')
;blnam=string(f, format='("blstore",i3.3,".sav")')
restore, filename=savdir+'kestore.sav', /verbose


;ds=getdata(f,/grid, wkdir=blah);

;grid=ds.grid
;delx=grid.x[1]-grid.x[0]
;dely=grid.y[1]-grid.y[0]
;delz=grid.z[1]-grid.z[0]
;
;x=grid.x
;nx = n_elements(x)-1
;dx = x(1)-x(0)
;y=grid.y
;ny = n_elements(y)-1
;dy = y(1)-y(0)
;z=grid.z
;nz = n_elements(z)-1
;dz = z(1)-z(0)

;data = getdata(f,wkdir=blah,/bx) & bx=data.bx
;data = getdata(f,wkdir=blah,/by) & by=data.by
;data = getdata(f,wkdir=blah,/bz) & bz=data.bz

;ztop=9.9
;zbottom=-9.9
;xtop=3.9
;xbottom=-1.9
;ytop=1.9
;ybottom=-1.9

;uz=where(abs(z-ztop) eq min(abs(z-ztop)))
;lz=where(abs(z-zbottom) eq min(abs(z-zbottom)))
;ux=where(abs(x-xtop) eq min(abs(x-xtop)))
;lx=where(abs(x-xbottom) eq min(abs(x-xbottom)))
;uy=where(abs(y-ytop) eq min(abs(y-ytop)))
;ly=where(abs(y-ybottom) eq min(abs(y-ybottom)))

;;STOP
;; making the box smaller to be xy=[-0.75,0.75], z=[-1,2]
;nx=ux-lx+1
;ny=uy-ly+1
;nz=uz-lz+1
;bgrid = dblarr(nx,ny,nz,3)
;bgrid[*,*,*,0] = bx(lx:ux,ly:uy,lz:uz)
;bgrid[*,*,*,1] = by(lx:ux,ly:uy,lz:uz)
;bgrid[*,*,*,2] = bz(lx:ux,ly:uy,lz:uz);

;print, 'destaggering'
;bgrid=destaggerB(bgrid)

;Move grid points so they are at the same locations as B in order to run 
;the null finding code
;xx = x(lx:ux)
;yy = y(ly:uy)
;zz = z(lz:uz);
nx=100
ny=nx
nz=ny

xx=6*findgen(nx)/99.-2
yy=4*findgen(ny)/99.-2
zz=20*findgen(nz)/99.-10

;undefine, data, bx, by, bz

;PRINT, 'modifying to ignore minimally accelerated particles'
;PRINT, 'now displaying:', n_elements(where(kestore-20.d0 ge 1)), 'particles'

;boost=1.0d0
;startpt=boost*startpt
;endpt=boost*endpt

;cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)
;epts=n_elements(estore)
;STOP
loadct, 4
;stretch, 30,270
tvlct, RR, GG, BB, /GET

;r3=alog10(1.0)
;t3=r3
ref=-2
topref=7
alke=alog10(kestore)
astore=fix((alke)/max(topref-ref)*254)
zerof=where(astore le 0)
astore(zerof)=0
t1=ref
t2=topref

; ch = 0.001
; chmin = 0.0001
; chmax = 0.01
; cepsilon = 1.0d-5
; cmxline = 10000
; csz_stp = size(startpt)


bscl=0.001d0
tscl=20.0d0
q=1.60217653E-19

;fac=arbscl*arbscl*bscl/tscl
;ofac=1.0d0/fac
;ofac10=10.0d0/fac
;ofac1000=1000.0d0/fac

;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t3*ofac,10^t2*ofac], title=mytit, /ylog
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t3*ofac100,10^t2*ofac100], title=mytit, /ylog

;STOP

mytit=' ';'energy gain *50keV [B/10G][L/100km]^2[T/20s]^(-1)'

xt='x (Mm [L/1Mm])'
yt='y (Mm [L/1Mm])'
zt='z (Mm [L/1Mm])'

undefine, data, bx, by, bz

 frac=0.02
 xlen=DOUBLE(max(xx)-min(xx))
 ylen=DOUBLE(max(yy)-min(yy))
 zlen=DOUBLE(max(zz)-min(zz))

 oOrb = OBJ_NEW('orb', COLOR=[170,170,170]) 	    	    	    	    	;initial orb - grey
 oOrb->Scale, frac, ylen/xlen*frac, zlen/xlen*frac
 
 
 oSymb = OBJ_NEW('IDLgrSymbol', oOrb) ;oSymbol for start point
 
; bresult = FILE_TEST(savdir+blnam) 
; IF (bresult eq 0) THEN BEGIN
;  PRINT, 'finding field lines..'
;  linefs=dblarr(3,cmxline,npts)
;  linebs=dblarr(3,cmxline,npts)
;  FOR i=1,npts DO BEGIN
;   linefs[*,*,i-1] = tracefield_3D(startpt[*,i-1]/arbscl,bgrid,xx,yy,zz,ch,chmin,chmax,cepsilon,cmxline)
;   linebs[*,*,i-1] = tracefield_3D(startpt[*,i-1]/arbscl,bgrid,xx,yy,zz,-ch,chmin,chmax,cepsilon,cmxline)
;  ENDFOR
;  save, linefs, linebs, filename=savdir+blnam, /verbose
; ENDIF ELSE BEGIN
;  restore, filename=savdir+blnam, /verbose
; ENDELSE

 
 frac2=frac+astore(0)/256.*0.05
 oOrb2 = OBJ_NEW('orb', COLOR=[rr(astore(0)),gg(astore(0)),bb(astore(0))]) 	;final orb - coloured with ke
 ;oOrb2->Scale, frac, frac, zlen/xlen*frac
 oOrb2->Scale, frac2, ylen/xlen*frac2, zlen/xlen*frac2
 oSymbol2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 for stop point
 ts=[ [startpt[0,0],startpt[0,0]],[startpt[1,0],startpt[1,0]],[startpt[2,0],startpt[2,0]]]
 te=[ [endpt[0,0],endpt[0,0]],[endpt[1,0],endpt[1,0]],[endpt[2,0],endpt[2,0]]]
 
 
iaz=0 
;FOR iaz=0,120 DO BEGIN
  print, 'frame', iaz, '/', 360
  xplot3d, [min(xx),min(xx)]*arbscl/lscl, [min(yy),min(yy)]*arbscl/lscl, [min(zz),min(zz)]*arbscl/lscl, $
  xr=[min(xx),max(xx)]*arbscl/lscl,yr=[min(yy),max(yy)]*arbscl/lscl,zr=[min(zz),max(zz)]*arbscl/lscl, $
  ;ax=ax, az=az, filewidth=filewidth, xtitle=xt, ytitle=yt, ztitle=zt
  ax=ax, az=30+iaz, filewidth=filewidth, xtitle=xt, ytitle=yt, ztitle=zt


;  jresult = FILE_TEST(savdir+fjnam) 
; IF (jresult eq 0) THEN BEGIN
;  data = getdata(f,wkdir=blah,/jx) 
;  jx=data.jx[lx:ux,ly:uy,lz:uz]
;  data = getdata(f,wkdir=blah,/jy) 
;  jy=data.jy[lx:ux,ly:uy,lz:uz]
;  data = getdata(f,wkdir=blah,/jz) 
;  jz=data.jz[lx:ux,ly:uy,lz:uz]
;  modj = sqrt(jx*jx+jy*jy+jz*jz)
;  save, modj, filename=savdir+fjnam
; ENDIF ELSE BEGIN
;  restore, filename=savdir+fjnam, /verbose
; ENDELSE
 
 
 
 
 xscale=(max(xx)-min(xx))/DOUBLE(nx);/lscl
 yscale=(max(yy)-min(yy))/DOUBLE(ny);/lscl
 zscale=(max(zz)-min(zz))/DOUBLE(nz);/lscl
 c1=[255,0,127]
 theModel = Obj_New('IDLgrModel')
 theModel -> Scale, xscale,yscale,zscale	    	    	
 ;theModel -> scale, frac, ylen/xlen*frac2, zlen/xlen*frac2
 
 ;jthresh=5.0
 ;IF (max(modj) gt jthresh) THEN BEGIN
 ; theModel -> Add, mk_iso(modj,jthresh,scol=c1, /low, minval=0.05)	
 ; test=obj_new('IDLgrSymbol', theModel)
 ; if obj_valid(test) THEN test->setproperty, size=1
 ; if obj_valid(test) THEN test->SetProperty, COLOR = c1
 ; xplot3d, [min(xx),min(xx)]*arbscl/lscl, [min(yy),min(yy)]*arbscl/lscl, [min(zz),min(zz)]*arbscl/lscl, symbol=test, /overplot
 ;ENDIF ELSE BEGIN
 ; print, 'value of jthresh not found in lare values of |j|'
 ;ENDELSE
;STOP
 
 FOR i=1,npts DO BEGIN
  frac2=frac+0.01*astore(i-1)/256.
  oOrb2 = OBJ_NEW('orb', COLOR=[rr(astore(i-1)),gg(astore(i-1)),bb(astore(i-1))]) 	;final orb - coloured with ke
  ;oOrb2->Scale, 1.1*frac, 1.1*frac, 1.1*zlen/xlen*frac
  oOrb2->Scale, frac2, ylen/xlen*frac2, zlen/xlen*frac2
  oSymb2 = OBJ_NEW('IDLgrSymbol', oOrb2) ;osymbol2 for stop point
  print, i, npts, format='(i4,"/",i4)'
  ts=[ [startpt[0,i-1],startpt[0,i-1]],[startpt[1,i-1],startpt[1,i-1]],[startpt[2,i-1],startpt[2,i-1]]]
  te=[ [endpt[0,i-1],endpt[0,i-1]],[endpt[1,i-1],endpt[1,i-1]],[endpt[2,i-1],endpt[2,i-1]]]
  ;XPLOT3D, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/lscl, COLOR=[0,0,0], NAME='start', SYMBOL=oSymb, THICK=2, /OVERPLOT
  ;print, ts[*,0]/lscl, ts[*,1]/lscl, ts[*,2]/lscl
  
  ;IF (astore(i-1) gt 128) THEN XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/lscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymbol2, THICK=2, /OVERPLOT
  XPLOT3D, te[*,0]/lscl, te[*,1]/lscl, te[*,2]/lscl, COLOR=[0,0,0], NAME='end', SYMBOL=oSymb2, THICK=2, /OVERPLOT
  ;IF (startpt[2,i-1] eq 200000.0d0) THEN BEGIN
  ;IF (((i-1) mod 32) eq 0) THEN BEGIN
  ; ii = where(linefs[0,*,i-1] gt -8000,nii)
  ; if (ii[0] ne -1) then xplot3d, linefs[0,ii,i-1]*arbscl/lscl, linefs[1,ii,i-1]*arbscl/lscl, linefs[2,ii,i-1]*arbscl/lscl, COLOR=grey, /OVERPLOT
  ; ii = where(linebs[0,*,i-1] gt -8000,nii)
  ; if (ii[0] ne -1) then xplot3d, linebs[0,ii,i-1]*arbscl/lscl, linebs[1,ii,i-1]*arbscl/lscl, linebs[2,ii,i-1]*arbscl/lscl, COLOR=grey, /OVERPLOT
  ;ENDIF
 ENDFOR


END
