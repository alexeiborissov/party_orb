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

;npts=2800
npts=1280

 nom="Data/"
 print, 'loading: ', nom
 ds=getrdata(1,/ke,/gyro)
 inike=ds.ke(0)
 ;loadct, 2
 estore=0
 kestore=dblarr(npts)
 tstop=dblarr(npts)
 startpt=dblarr(3,npts)
 endpt=dblarr(3,npts)
 mrg=dblarr(npts)
 mrp=dblarr(npts)
 ns=n_elements(ds.t)
 kestore[0]=inike
 startpt[*,0]=[ds.x[0],ds.y[0],ds.z[0]]
 endpt[*,0]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
 tstop[0]=ds.t[ns-1]
 mrg[0]=max(ds.gyror)
 mrp[0]=max(ds.gyrop)
 for i=1,npts DO BEGIN
  tflag=0
  print, i, npts, format='(i4,"/",i4)'
  ds=getrdata(i,/ke,/gyro)
  ns=n_elements(ds.t)
  startpt[*,i-1]=[ds.x(0),ds.y(0),ds.z(0)]
  endpt[*,i-1]=[ds.x[ns-1],ds.y[ns-1],ds.z[ns-1]]
  tstop[i-1]=ds.t[ns-1]
  maxke=max(ds.ke)
  kestore[i-1]=maxke
  mrg[i-1]=max(ds.gyror)
  mrp[i-1]=max(ds.gyrop)
  print, startpt[*,i-1], 'max ke=', maxke
 endfor

END
