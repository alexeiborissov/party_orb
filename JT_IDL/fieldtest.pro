;program to test the value of B found at a certain position to the value used by Fortran
print, 'getting lare data..'
blah='../data/VAR378_B.xdr'
nx=1024
ny=1024
nz=256

restore, filename=blah, /verbose
print, 'getting particle data..'
ds=getrdata(1,/fields)

inipos=[ds.x(0),ds.y(0),ds.z(0)]


print, min(abs(ds.x(0)-x),minloc)
print, 'particle begins at', inipos, format='(a,"[",e9.2,",",e9.2,",",e9.2,"]")'

tempx=min(abs(ds.x(0)-x),minlocx)
if ((ds.x(0) ge x(minlocx)) and (ds.x(0) le x(minlocx+1))) THEN xloc=minlocx
if ((ds.x(0) ge x(minlocx-1)) and (ds.x(0) le x(minlocx))) THEN xloc=minlocx-1
tempy=min(abs(ds.y(0)-y),minlocy)
if ((ds.y(0) ge y(minlocy)) and (ds.y(0) le y(minlocy+1))) THEN yloc=minlocy
if ((ds.y(0) ge y(minlocy-1)) and (ds.y(0) le y(minlocy))) THEN yloc=minlocy-1
tempz=min(abs(ds.z(0)-z),minlocz)
if ((ds.z(0) ge z(minlocz)) and (ds.z(0) le z(minlocz+1))) THEN zloc=minlocz
if ((ds.z(0) ge z(minlocz-1)) and (ds.z(0) le z(minlocz))) THEN zloc=minlocz-1

delx=x[1:nx-1]-x[0:nx-2]
dely=y[1:ny-1]-y[0:ny-2]
delz=z[1:nz-1]-z[0:nz-2]

coarseloc=[xloc,yloc,zloc]

dx=(inipos[0]-x[xloc])/delx[xloc]
dy=(inipos[1]-y[yloc])/dely[yloc]
dz=(inipos[2]-z[zloc])/delz[zloc]

dr=[dx,dy,dz]

interptempbx=vdcalc3d(B[*,*,*,0], xloc,yloc,zloc,dx,dy,dz)
interptempby=vdcalc3d(B[*,*,*,1], xloc,yloc,zloc,dx,dy,dz)
interptempbz=vdcalc3d(B[*,*,*,2], xloc,yloc,zloc,dx,dy,dz)

print, 'IDL says B should be:', interptempbx[0],interptempbz[0],interptempbz[0]
print, 'Fortran B is:', ds.B[0,0],ds.B[1,0],ds.B[2,0]
END
;xloc=0
;yloc=0
;zloc=0

;print, 'cell no.', xloc, yloc, zloc
;print, 'bx(0,0,0)=',b[xloc,yloc,zloc,0], format='(a,e12.5)'
;print, 'bx(1,0,0)=',b[xloc+1,yloc,zloc,0], format='(a,e12.5)'
;print, 'bx(0,1,0)=',b[xloc,yloc+1,zloc,0], format='(a,e12.5)'
;print, 'bx(0,0,1)=',b[xloc,yloc,zloc+1,0], format='(a,e12.5)'
;print, 'bx(1,1,0)=',b[xloc+1,yloc+1,zloc,0], format='(a,e12.5)'
;print, 'bx(0,1,1)=',b[xloc,yloc+1,zloc+1,0], format='(a,e12.5)'
;print, 'bx(1,0,1)=',b[xloc+1,yloc,zloc+1,0], format='(a,e12.5)'
;print, 'bx(1,1,1)=',b[xloc+1,yloc+1,zloc+1,0], format='(a,e12.5)'
;print, '----------'
;print, 'by(0,0,0)=',b[xloc,yloc,zloc,1], format='(a,e12.5)'
;print, 'by(1,0,0)=',b[xloc+1,yloc,zloc,1], format='(a,e12.5)'
;print, 'by(0,1,0)=',b[xloc,yloc+1,zloc,1], format='(a,e12.5)'
;print, 'by(0,0,1)=',b[xloc,yloc,zloc+1,1], format='(a,e12.5)'
;print, 'by(1,1,0)=',b[xloc+1,yloc+1,zloc,1], format='(a,e12.5)'
;print, 'by(0,1,1)=',b[xloc,yloc+1,zloc+1,1], format='(a,e12.5)'
;print, 'by(1,0,1)=',b[xloc+1,yloc,zloc+1,1], format='(a,e12.5)'
;print, 'by(1,1,1)=',b[xloc+1,yloc+1,zloc+1,1], format='(a,e12.5)'
;print, '----------'
;print, 'bz(0,0,0)=',b[xloc,yloc,zloc,2], format='(a,e12.5)'
;print, 'bz(1,0,0)=',b[xloc+1,yloc,zloc,2], format='(a,e12.5)'
;print, 'bz(0,1,0)=',b[xloc,yloc+1,zloc,2], format='(a,e12.5)'
;print, 'bz(0,0,1)=',b[xloc,yloc,zloc+1,2], format='(a,e12.5)'
;print, 'bz(1,1,0)=',b[xloc+1,yloc+1,zloc,2], format='(a,e12.5)'
;print, 'bz(0,1,1)=',b[xloc,yloc+1,zloc+1,2], format='(a,e12.5)'
;print, 'bz(1,0,1)=',b[xloc+1,yloc,zloc+1,2], format='(a,e12.5)'
;print, 'bz(1,1,1)=',b[xloc+1,yloc+1,zloc+1,2], format='(a,e12.5)'
;STOP
