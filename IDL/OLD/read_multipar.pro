PRO read_multipar,filename,xstart=xstart,ystart=ystart,zstart=zstart,tstart=tstart,ekinstart=ekinstart,alphastart=alphastart,eperpstart=eperpstart,eparstart=eparstart,vtotstart=vtotstart,vperpstart=vperpstart,vparstart=vparstart,Exstart=Exstart,Eystart=Eystart,Ezstart=Ezstart,Bxstart=Bxstart,Bystart=Bystart,Bzstart=Bzstart,BStart=BStart,mu=mu,xend=xend,yend=yend,zend=zend,tend=tend,ekinend=ekinend,alphaend=alphaend,eperpend=eperpend,eparend=eparend,vtotend=vtotend,vperpend=vperpend,vparend=vparend,Exend=Exend,Eyend=Eyend,Ezend=Ezend,Bxend=Bxend,Byend=Byend,Bzend=Bzend,Bend=Bend,unlost=unlost,lost=lost,nounlost=nounlost,countlostparticles=countlostparticles


;changed to use double precision

; !this inputs: at start: x,y,z, vpar
; !              at end: mu, tt, x, y, z, vpar



n=file_lines(filename)   ;22869
A=dblarr(37,n)
openr,lun,filename,/get_lun
readf,lun,A
free_lun,lun
xstart=dblarr(n,/nozero)
xstart=transpose(A(0,*))
ystart=dblarr(n,/nozero)
ystart=transpose(A(1,*))
zstart=dblarr(n,/nozero)
zstart=transpose(A(2,*))

tstart=dblarr(n,/nozero)
tstart=transpose(A(3,*))

ekinstart=dblarr(n,/nozero)
ekinstart=transpose(A(4,*))

alphastart=dblarr(n,/nozero)
alphastart=transpose(A(5,*))

eperpstart=dblarr(n,/nozero)
eperpstart=transpose(A(6,*))

eparstart=dblarr(n,/nozero)
eparstart=transpose(A(7,*))

vtotstart=dblarr(n,/nozero)
vtotstart=transpose(A(8,*))

vperpstart=dblarr(n,/nozero)
vperpstart=transpose(A(9,*))

;----------------------------------

vparstart=dblarr(n,/nozero)
vparstart=transpose(A(10,*))

Exstart=dblarr(n,/nozero)
Exstart=transpose(A(11,*))
Eystart=dblarr(n,/nozero)
Eystart=transpose(A(12,*))
Ezstart=dblarr(n,/nozero)
Ezstart=transpose(A(13,*))

Bxstart=dblarr(n,/nozero)
Bxstart=transpose(A(14,*))
Bystart=dblarr(n,/nozero)
Bystart=transpose(A(15,*))
Bzstart=dblarr(n,/nozero)
Bzstart=transpose(A(16,*))
Bstart=dblarr(n,/nozero)
Bstart=transpose(A(17,*))

mu=dblarr(n,/nozero)
mu=transpose(A(18,*))

xend=dblarr(n,/nozero)
xend=transpose(A(19,*))
yend=dblarr(n,/nozero)
yend=transpose(A(20,*))
zend=dblarr(n,/nozero)
zend=transpose(A(21,*))

tend=dblarr(n,/nozero)
tend=transpose(A(22,*))

ekinend=dblarr(n,/nozero)
ekinend=transpose(A(23,*))

alphaend=dblarr(n,/nozero)
alphaend=transpose(A(24,*))

eperpend=dblarr(n,/nozero)
eperpend=transpose(A(25,*))

eparend=dblarr(n,/nozero)
eparend=transpose(A(26,*))

vtotend=dblarr(n,/nozero)
vtotend=transpose(A(27,*))

vperpend=dblarr(n,/nozero)
vperpend=transpose(A(28,*))

vparend=dblarr(n,/nozero)
vparend=transpose(A(29,*))

Exend=dblarr(n,/nozero)
Exend=transpose(A(30,*))
Eyend=dblarr(n,/nozero)
Eyend=transpose(A(31,*))
Ezend=dblarr(n,/nozero)
Ezend=transpose(A(32,*))

Bxend=dblarr(n,/nozero)
Bxend=transpose(A(33,*))
Byend=dblarr(n,/nozero)
Byend=transpose(A(34,*))
Bzend=dblarr(n,/nozero)
Bzend=transpose(A(35,*))
Bend=dblarr(n,/nozero)
Bend=transpose(A(36,*))

;id's of particles still in the trap (not lost to photosphere)
;so can plot, say,
;plot,xend[unlost],yend[unlost],psym=3
unlost=where(tend eq max(tend),nounlost)
lost=where(tend lt max(tend), countlostparticles)



end
