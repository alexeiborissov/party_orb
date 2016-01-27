@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1

loadct, 2, /silent

f=41

print, 'getting data..'
;ds=getdata(41,/vx,/vy,/vz,/bx,/by,/bz)
ds=getdata(f,/grid, wkdir="../../julie_dam_files/Data")

grid=ds.grid
delx=grid.x[1]-grid.x[0]
dely=grid.y[1]-grid.y[0]
delz=grid.z[1]-grid.z[0]

jnom="../../julie_dam_files/"+string(f, format='("newJ", I2,".sav")')

restore, jnom, /verbose
modj=sqrt(j[*,*,*,0]*j[*,*,*,0]+j[*,*,*,1]*j[*,*,*,1]+j[*,*,*,2]*j[*,*,*,2])

window, 0
contour, modj[*,*,368], grid.x, grid.y, /fill, nlev=20, /iso, xr=[-0.6,0.6], yr=[-0.6,0.6], title='|j(z=0)|', ytitle='y', xtitle='x'
loadct, 4, /silent
tvlct, r, g, b, /get

npts=21*1*1;+1
nom="Data/"
;nom="../relcodev1.0/Data/"
;nom="../teq0/Data/"

print, 'loading: ', nom
ds=getrdata(1,/ke)
inike=ds.ke(0)

kestore=dblarr(npts)
for i=1,npts DO BEGIN
 print, i, npts, format='(i4,"/",i4)'
 ds=getrdata(i,/ke)  
 maxke=max(ds.ke)
 ;maxketime=ds.time(maxkeloc)
 kestore(i-1)=maxke
; if ds.z[0] eq 0.5e6 THEN oplot, [ds.x[0], ds.x[0]]*1e-6, [ds.y[0],ds.y[0]]*1e-6, thick=2, psym=4, symsize=3
 if ds.z[0] eq 0.5e4 THEN oplot, [ds.x[0], ds.x[0]]*1e-4, [ds.y[0],ds.y[0]]*1e-4, thick=2, psym=4, symsize=2
endfor
cstore=fix((kestore-min(kestore))/max(kestore-min(kestore))*254)

;save, kestore, cstore, /compress, /verbose, filename=nom+'../IDL/colorstore_test.sav'
;restore, filename=nom+'../IDL/colorstore_test.sav',/verbose

ref=alog10(1.0e-12)

alke=alog10(kestore-2.0d0)
;astore=fix((alke-min(alke))/max(alke-min(alke))*254)
astore=fix((alke-ref)/max(alke-ref)*254)

t1=ref
;t1=alog10(min(kestore)-2.0d0)
t2=alog10(max(kestore)-2.0d0)

;make a new window for the colorbar
;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], title='Particle peak energy gain (eV)', /vertical, /ylog

;make a new window for the colorbar
cgwindow
cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], $ 
title='Particle peak energy gain (eV)', /vertical, /ylog

;particletrack, 1, floc=nom, xyzt, lcol=[0,0,0], zscl=1.0, /Mm, myxyrange=[-1.00,1.00], myzrange=[-1.75,2.75], /rel


particletrack, 1, floc=nom, xyzt, lcol=[r(astore(0)),g(astore(0)),b(astore(0))], zscl=1e3, $
;/km, myxyrange=[-10.0,10.0], myzrange=[-17.5,27.5], /rel, /symb
/km, myxyrange=[-10.0,10.0], myzrange=[-5.0,15.0], /rel, /symb


for i=2,npts DO BEGIN
 particletrack, i, /op, floc=nom, xyzt, lcol=[r(astore(i-1)),g(astore(i-1)),b(astore(i-1))], zscl=1e3, /km, /rel, /symb,myxyrange=[-10.0,10.0], myzrange=[-5.0,15.0]
endfor

stop
;if modjiso then begin 

 ds=getdata(f,/vx, /vy, /vz, wkdir="../../julie_dam_files/Data")
 modv=sqrt(ds.vx*ds.vx+ds.vy*ds.vy+ds.vz*ds.vz)
 modv[0:128,*,*]=0.0d0
 modv[*,0:128,*]=0.0d0
 modv[*,384:512,*]=0.0d0
 modv[384:512,*,*]=0.0d0
 
 lsize=10
 theModel = Obj_New('IDLgrModel')
 ;theModel -> Scale, 2./512, 2./512, 4.5/768	    	; julie's box dimensions
 theModel -> Scale, 2./512, 2./512, 2.0/768
 theModel -> Add, mk_iso(modv,0.0005,scol=[0,240,0], maxval=0.001)	; note threshold choice as 2nd arg of mk_iso
 test=obj_new('IDLgrSymbol', theModel)
 if obj_valid(test) THEN test->setproperty, size=lsize
 xplot3d, [-1*lsize,-1*lsize], [-1*lsize,-1*lsize], [-0.5*lsize,-0.5*lsize], symbol=test, /overplot
;endif




END
