;program to only plot the start and end points of particle traces at different times, TOGETHER WITH FAN PLANES, NULLS AND SPINES.
;@ninfo3d.pro
eps=0.0000001
h1=0.01
string1='derivs'
FRon=1
;@xplot3dJT
loadct, 4, /silent
tvlct, r, g, b, /get

f=1059
filename='test.png'
ax=-76
az=15
filewidth=1000
transparent=0.5
out=0

od='./figs/'

;npts=2800
;npts=100
 npts=16*16*16

;FOR e=10,10 DO BEGIN
;tempdir='../../electrons_maxwellian/epsilon'
tempdir='.'
;es=string(e,format='(i2.2)')
;print, 'sorting epsilon '+es

result = FILE_TEST(tempdir+'/saves/kestore.sav') 
IF (result eq 0) THEN BEGIN
 print, 'loading electron data'
 ds=getrdata(1,/ke,/gyro, wkdir=tempdir+'/Data/DataR/')
 ns=n_elements(ds.t)
 ;loadct, 2
 estore=0
 kestore=dblarr(npts)
 inike=dblarr(npts)
 finke=dblarr(npts)
 tstop=dblarr(npts)
 startpt=dblarr(3,npts)
 endpt=dblarr(3,npts)
 mrg=dblarr(npts)
 mrp=dblarr(npts)
 ns=n_elements(ds.t)
 kestore[0]=ds.ke(0)
 inike[0]=ds.ke(0)
 finke[0]=ds.ke[ns-1]
 startpt[*,0]=[ds.x[0],ds.y[0],ds.z[0]]
 endpt[*,0]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
 tstop[0]=ds.t[ns-1]
 mrg[0]=max(ds.gyror)
 mrp[0]=max(ds.gyrop)
 for i=1,npts-1 DO BEGIN
  tflag=0
  print, i, npts, format='(i4,"/",i4)'
  ds=getrdata(i,/ke,/gyro, wkdir=tempdir+'/Data/DataR/')
  ns=n_elements(ds.t)
  inike[i]=ds.ke[0]
  finke[i]=ds.ke[ns-1]
  startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
  endpt[*,i-1]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
  tstop[i-1]=ds.t[ns-1]
  maxke=max(ds.ke)
  kestore[i-1]=maxke
  mrg[i-1]=max(ds.gyror)
  mrp[i-1]=max(ds.gyrop)
  print, startpt[*,i-1], 'max ke=', maxke
 endfor
 save, npts, kestore, startpt, endpt, tstop, inike, finke, filename=tempdir+'/saves/kestore.sav'
ENDIF ELSE BEGIN
;STOP
restore, filename=tempdir+'/saves/kestore.sav', /verbose
ENDELSE

;data = data(0:i-1) ; trim the data at the end.
;STOP
;!p.multi=[0,1,2]
;plot, findgen(Npts), inike
loadct, 3
;stop
t=1e6
bsz = 0.2
minEN=1

;myEN = sqrt(data)
myEN=inike;*6.242e18  ;(convert to eV)
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
print, 'min/max initial electron energy (log eV):',min(l_EN),max(l_EN)
nmyEN = n_elements(myEN)
;h = histogram(l_EN,binsize=bsz,min=lf_sz_min, locations=xbin)
h = histogram(l_EN,binsize=bsz, locations=xbin)
n_h = n_elements(h)

myEN=finke;*6.242e18  ;(convert to eV)
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
print, 'min/max final electron energy (log eV):',min(l_EN),max(l_EN)
nmyEN = n_elements(myEN)
;h = histogram(l_EN,binsize=bsz,min=lf_sz_min, locations=xbin)
h2 = histogram(l_EN,binsize=bsz, locations=xbin2)
n_h = n_elements(h)


tempx=10^xbin
tempy=h
ntempx = N_ELEMENTS(tempx)
histo_tempx1 = CONGRID(tempx, ntempx*2)
histo_tempy1 = CONGRID(tempy, ntempx*2, /CENTER)
tempx=10^xbin2
tempy=h2
ntempx = N_ELEMENTS(tempx)
histo_tempx2 = CONGRID(tempx, ntempx*2)
histo_tempy2 = CONGRID(tempy, ntempx*2, /CENTER)

 aspect_ratio=1.8
 myxs=20
 IF NOT(OUT) THEN myxs=myxs*50
 myys=myxs/aspect_ratio
loadct, 39

IF (out) THEN BEGIN	; set up plot area
  set_plot,'ps'
  !p.font=0
  device, filename=od+'pespectra.eps', encapsulated=1, /helvetica, /color, BITS_PER_PIXEL=8
  device, xsize=myxs, ysize=myys
  mycharsize=1
  myth=4
  PRINT, 'output to: '+od+'pespectra.eps'
 ENDIF ELSE BEGIN
  window, 5, ysize=myys, xsize=myxs
  mycharsize=1
  myth=2
 ENDELSE 


plot, histo_tempx1, histo_tempy1, xtitle='energy (eV)', ytitle='N', /xlog, xr=[1,1e6], /ylog, yr=[1e0,1e4]
oplot, histo_tempx2, histo_tempy2, col=240

;oplot, 10^xbin, 0.0001*(10^xbin)^1.5, thick=4



oplot, [t/11604.505,t/11604.505], [1,1e8], linestyle=1
xyouts, 10, 2, "KE (T=1E6K)", /DATA, charsize=0.9

;legend, ['ini electron KE','ini proton KE','final electron KE','final proton KE'], linestyle=[0,2,0,2], color=[0,0,240,240], /bottom, /right, charsize=0.9


 IF (out) THEN BEGIN
  device, /close
  set_plot, 'x'
  !p.font=-1
 ENDIF

END
