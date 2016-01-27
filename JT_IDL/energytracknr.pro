;common somename, tt, FRon
@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 4, /silent
tvlct, r, g, b, /get

snap=41
snapnom=string(snap,format='("l3dsnap",i2,".ninfo")')

;print, "loading lare B-fields"
;ds=getdata(snap,/bx, /by, /bz, wkdir="../julie_dat_files/Data")

;laregrid=ds.grid

;print, "correcting for staggered grid"
;larefieldcorrect, ds.bx,ds.by,ds.bz, ds.grid.npts, newBx,newBy,newBz
;undefine, ds

;bv=dblarr(laregrid.npts[0],laregrid.npts[1],laregrid.npts[2],3)
;bv[*,*,*,0]=newbx
;bv[*,*,*,1]=newby
;bv[*,*,*,2]=newbz
;undefine, newbx, newby, newbz;

;nulls=nullfinder(bv)
;newton_improve, bv, nulls
;print, 'andrews null finder:'
;print, nulls
;print, '--'

;ninfo3d, reform(bv[*,*,*,0]), reform(bv[*,*,*,1]), reform(bv[*,*,*,2]), laregrid.x, laregrid.y, laregrid.z, laregrid.npts, ninfo, ncan, px, py, pz, candidates
;ninfo3d, ds.bx, ds.by, ds.bz, laregrid.x, laregrid.y, laregrid.z, laregrid.npts, ninfo, ncan, px, py, pz, candidates
;print, "ncandidates:", ncan
;print, "coords:", ninfo.coords

;STOP
;save, ncan, ninfo, filename=nom+snapnom, /compress
npts=64
;npts=1280
nom="Data/"
;nom="../relcodev1.0/Data/"
;nom="../teq0/Data/"

print, 'loading: ', nom
ds=getndata(1,/ke)
inike=ds.ke(0)

kestore=dblarr(npts)
for i=1,npts DO BEGIN
 print, i, npts, format='(i4,"/",i4)'
 ds=getndata(i,/ke)  
 maxke=max(ds.ke)
 ;maxketime=ds.time(maxkeloc)
 kestore(i-1)=maxke
 print, maxke
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

;save, kestore, cstore, /compress, /verbose, filename=nom+'../IDL/colorstore_test.sav'
;restore, filename=nom+'../IDL/colorstore_test.sav',/verbose

ref=alog10(1.0e-12)

alke=alog10(kestore-1.99999999d0)
;astore=fix((alke-min(alke))/max(alke-min(alke))*254)
astore=fix((alke-ref)/max(alke-ref)*254)

t1=ref
;t1=alog10(min(kestore)-2.0d0)
t2=alog10(max(kestore)-20.0d0)

;make a new window for the colorbar
;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], title='Particle peak energy gain (eV)', /vertical, /ylog

;make a new window for the colorbar
;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
;title='Particle peak energy gain (eV)', /vertical, /ylog

;particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], zscl=1.0, /Mm, myxyrange=[-1.00,1.00], myzrange=[-1.75,2.75], /rel

;stop
;particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e6, /mm, myxyrange=[-100.0,100.0], myzrange=[-60,60], /rel, ax=-60, az=40, filewidth=800
particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], lscl=1e6, zscl=1e6, /mm, myxyrange=[-100.0,100.0], myzrange=[-60,60], ax=-60, az=40, filewidth=800

for i=2,npts DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e6, /mm, myxyrange=[-100.0,100.0], myzrange=[-60,60]
endfor


;STOP

iniorbplot, 1, op=0, floc=nom, ocol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e6, /km, myxyrange=[-0.4,0.4], myzrange=[-20,20], ax=-77, az=30, filewidth=800
;iniorbplot, 1, op=0, floc=nom, ocol=[0,0,0], zscl=1e6, /km, /rel, myxyrange=[-0.4,0.4], myzrange=[-20,20], ax=-77, az=30, filewidth=800

for i=2,npts DO BEGIN
print, i, npts, format='(i4,"/",i4)'
 iniorbplot, i, /op, floc=nom, ocol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e6, /km
 ;iniorbplot, i, /op, floc=nom, ocol=[0,0,0], zscl=1e6, /km, /rel
endfor

STOP

;plotting the nulls as orbs:
 restore, filename='IDL/'+snapnom, /verbose
 ds=getdata(41,/grid, wkdir="../../julie_dat_files/Data")
 xst=ds.grid.x(0)
 yst=ds.grid.y(0)
 zst=ds.grid.z(0)
 ldx=ds.grid.x(1)-ds.grid.x(0)
 ldy=ds.grid.y(1)-ds.grid.y(0)
 ldz=ds.grid.z(1)-ds.grid.z(0)

 lscl=1e7
 rmult=1.0
 xt='x (Mm)'
 yt='y (Mm)'
 xyrat=1.
 frac=0.3
 xylen=2.0
 zlen=(2.75+1.75)

 oOrb = OBJ_NEW('orb', COLOR=[0, 0 ,255]) 	;blue orb
 oOrb->Scale, frac*xyrat, frac*xyrat, zlen/xylen*frac
 oSymbol = OBJ_NEW('IDLgrSymbol', oOrb) ;make the orb into a symbol

;position array
 n1=[[nulls[0,0]*ldx+xst,nulls[0,0]*ldx+xst],[nulls[0,1]*ldy+yst,nulls[0,1]*ldy+yst],[nulls[0,2]*ldz+zst,nulls[0,2]*ldz+zst]]
 n2=[[nulls[1,0]*ldx+xst,nulls[1,0]*ldx+xst],[nulls[1,1]*ldy+yst,nulls[1,1]*ldy+yst],[nulls[1,2]*ldz+zst,nulls[1,2]*ldz+zst]]

 ;help, te
 ;XPLOT3D, n1[*,0]*10.0, n1[*,1]*10.0, n1[*,2]*10.0, COLOR=[0,0,0], NAME='n1', SYMBOL=oSymbol, THICK=2, /OVERPLOT
 ;XPLOT3D, n2[*,0]*10.0, n2[*,1]*10.0, n2[*,2]*10.0, COLOR=[0,0,0], NAME='n2', SYMBOL=oSymbol, THICK=2, /OVERPLOT

modvxb=0
Epariso=0
modBiso=1
modJiso=0

lsize=10

if modvxb then begin
;isosurface of Epar?
 restore, filename='../../julie_dat_files/modVxB41.sav', /verbose
 theModel = Obj_New('IDLgrModel')
 theModel -> Scale, 2./512, 2./512, 4.5/768	    	; julie's box dimensions
 theModel -> Add, mk_iso(modvxb,0.001,scol=[255,255,255])	; note threshold choice as 2nd arg of mk_iso
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=lsize
 xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.75*lsize,-1.75*lsize], symbol=test, /overplot
endif

; stop 
 if epariso then begin
;isosurface of Epar?
 restore, filename='../../julie_dat_files/Epar41.sav', /verbose
 theModel = Obj_New('IDLgrModel')
 theModel -> Scale, 2./512, 2./512, 4.5/768	    	; julie's box dimensions
 theModel -> Add, mk_iso(Epar,0.00005,scol=[255,255,255])	; note threshold choice as 2nd arg of mk_iso
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=lsize
 xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.75*lsize,-1.75*lsize], symbol=test, /overplot
 endif
 
 
 ;STOP
if modbiso then begin
;isosurface of b?
 restore, filename='../../julie_dat_files/modB41.sav', /verbose
 theModel = Obj_New('IDLgrModel')
 theModel -> Scale, 2./512, 2./512, 4.5/768	    	; julie's box dimensions
 theModel -> Add, mk_iso(modB,0.2,scol=[255,255,255], /low)	; note threshold choice as 2nd arg of mk_iso
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=lsize
 if obj_valid(test) THEN test->SetProperty, COLOR = [0, 240, 0]
 xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.75*lsize,-1.75*lsize], symbol=test, /overplot
endif
 
 
if modjiso then begin 
;isosurface of current? 
 ;restore, filename='../../julie_dat_files/Jtest41.sav', /verbose
 ;modJ=sqrt(J[*,*,*,0]*J[*,*,*,0]+J[*,*,*,1]*J[*,*,*,1]+J[*,*,*,2]*J[*,*,*,2])
 ;save, modJ, filename='../julie_dat_files/modJ41.sav', /compress, /verbose
 restore, filename='../julie_dat_files/modJ41.sav', /verbose 
 theModel = Obj_New('IDLgrModel')
 theModel -> Scale, 2./512, 2./512, 4.5/768	    	; julie's box dimensions
 theModel -> Add, mk_iso(modJ,0.0001,scol=[0,240,0], maxval=0.001)	; note threshold choice as 2nd arg of mk_iso
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=lsize
 xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-1.75*lsize,-1.75*lsize], symbol=test, /overplot
endif   
   
STOP



iniorbplot, 1, op=0, floc=nom, ocol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e6, /me, /rel, myxyrange=[-40.,40.], myzrange=[-20,20]

for i=2,npts DO BEGIN
print, i, npts, format='(i4,"/",i4)'
 iniorbplot, i, /op, floc=nom, ocol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e6, /me, /rel, myxyrange=[-40.,40.], myzrange=[-20,20]
endfor

        
END
