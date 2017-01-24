;PRO FLINE
; FLINE
;  numerically integrates field lines according to field specified in derivs
;  then displays
;  REQUIRES - ODEINT, derivs, RKQC
@xplot3dJT
@ODEINT
;@ps_fonts
common somename, epsilon

tempdir2='../../protons/epsilon
tempdir='../epsilon'
newoutdir='./allepsilon'

 npts=11*11*5
 dims=[1.,1.,1.]
 mysymsize=[1.0,1.0,1.0]*0.6*dims  

;;----- define variables ----------------
eps=0.0000001
h1=0.0001     	    ; step size?
string1='derivs2'

nfolders=21
bigkestore=dblarr(2,nfolders)
avgkestore=dblarr(2,nfolders)

FOR e=0, nfolders-1 DO BEGIN

es=string(e,format='(i2.2)')
print, 'sorting epsilon '+es

result = FILE_TEST(tempdir+es+'/kestore.sav') 
IF (result eq 0) THEN BEGIN
 nom=tempdir+es+'/Data/'
 print, 'loading: ', nom
 ds=getrdata(1,/ke, /gyro, wkdir=nom+'DataR/')
 inike=ds.ke(0)
 estore=0
 kestore=dblarr(npts)
 maxgyro=dblarr(npts)
 tstop=dblarr(npts)
 startpt=dblarr(3,npts)
 endpt=dblarr(3,npts)
 ns=n_elements(ds.t)
 kestore[0]=inike
 maxg=max(ds.gyror)
 maxgyro[0]=maxg
 startpt[*,0]=[ds.x[0],ds.y[0],ds.z[0]]
 endpt[*,0]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
 tstop[0]=ds.t[ns-1]
 for i=1,npts DO BEGIN
  tflag=0
  print, i, npts, format='(i4,"/",i4)'
  ds=getrdata(i,/ke,/gyro, wkdir=nom+'DataR/')
  ns=n_elements(ds.t)
  startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
  endpt[*,i-1]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
  tstop[i-1]=ds.t[ns-1]
  maxke=max(ds.ke)
  kestore[i-1]=maxke
  print, startpt[*,i-1], 'max ke=', maxke
  maxg=max(ds.gyror)
  maxgyro[i-1]=maxg
 endfor
 save, npts, kestore, startpt, endpt, tstop, maxgyro, filename=tempdir+es+'/kestore.sav'
ENDIF ELSE BEGIN
restore, filename=tempdir+es+'/kestore.sav', /verbose
ENDELSE
bigkestore(0,e)=max(kestore)
avgkestore(0,e)=mean(kestore)

result2 = FILE_TEST(tempdir2+es+'/kestore.sav') 
IF (result2 eq 0) THEN BEGIN
 nom=tempdir2+es+'/Data/'
 print, 'loading: ', nom
 ds=getrdata(1,/ke, /gyro, wkdir=nom+'DataR/')
 inike=ds.ke(0)
 estore=0
 kestore=dblarr(npts)
 maxgyro=dblarr(npts)
 tstop=dblarr(npts)
 startpt=dblarr(3,npts)
 endpt=dblarr(3,npts)
 ns=n_elements(ds.t)
 kestore[0]=inike
 maxg=max(ds.gyror)
 maxgyro[0]=maxg
 startpt[*,0]=[ds.x[0],ds.y[0],ds.z[0]]
 endpt[*,0]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
 tstop[0]=ds.t[ns-1]
 for i=1,npts DO BEGIN
  tflag=0
  print, i, npts, format='(i4,"/",i4)'
  ds=getrdata(i,/ke,/gyro, wkdir=nom+'DataR/')
  ns=n_elements(ds.t)
  startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
  endpt[*,i-1]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
  tstop[i-1]=ds.t[ns-1]
  maxke=max(ds.ke)
  kestore[i-1]=maxke
  print, startpt[*,i-1], 'max ke=', maxke
  maxg=max(ds.gyror)
  maxgyro[i-1]=maxg
 endfor
 save, npts, kestore, startpt, endpt, tstop, maxgyro, filename=tempdir2+es+'/kestore.sav'
ENDIF ELSE BEGIN
restore, filename=tempdir2+es+'/kestore.sav', /verbose
ENDELSE
bigkestore(1,e)=max(kestore)
avgkestore(1,e)=mean(kestore)
ENDFOR

!p.background=255
!p.font = 1
loadct, 39
psa, file = 'kemax.eps', xs = 16, ys = 9, /inches, /encaps, /color, bits_per_pixel=24
;window, 0
!p.font = 1
plot, findgen(nfolders), bigkestore[0,*]/1e6, xtitle="epsilon(t)", ytitle="kinetic energy gain [MeV]", psym=7, col=0, charsize=3, symsize=4, thick=5, yr=[0.05,5], /ylog
oplot, findgen(nfolders), bigkestore[1,*]/1e6, psym=6, linestyle=2, col=0, symsize=4, thick=5
oplot, findgen(nfolders), avgkestore[0,*]/1e6, psym=7, linestyle=2, col=240, symsize=4, thick=5
oplot, findgen(nfolders), avgkestore[1,*]/1e6, psym=6, linestyle=2, col=240, symsize=4, thick=5
legend, ["electron","proton"], psym=[7,6],col=[0,0], /bottom, /right, textcolors=0, charsize=3, symsize=[4,4], thick=5, box=0
legend, ["max","average"], linestyle=[0,0],col=[0,240], /top, /left, textcolors=0, charsize=3, symsize=[4,4], thick=5, box=0
pse
!p.background=0
;ENDFOR
t1=3
t2=6.3
loadct, 4
cgwindow
cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t1,10^t2], /vertical, /ylog 
;title='Particle peak energy gain (eV)', /vertical, /ylog


end
