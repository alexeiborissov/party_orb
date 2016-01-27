FUNCTION loaddata_mks_n,filename


;note units are MKS except energy is in eV.
;filename='multipar_fmt_t10.nohead'

read_multipar,filename,xstart=xstart,ystart=ystart,zstart=zstart,tstart=tstart,ekinstart=ekinstart,alphastart=alphastart,eperpstart=eperpstart,eparstart=eparstart,vtotstart=vtotstart,vperpstart=vperpstart,vparstart=vparstart,Exstart=Exstart,Eystart=Eystart,Ezstart=Ezstart,Bxstart=Bxstart,Bystart=Bystart,Bzstart=Bzstart,BStart=BStart,mu=mu,xend=xend,yend=yend,zend=zend,tend=tend,ekinend=ekinend,alphaend=alphaend,eperpend=eperpend,eparend=eparend,vtotend=vtotend,vperpend=vperpend,vparend=vparend,Exend=Exend,Eyend=Eyend,Ezend=Ezend,Bxend=Bxend,Byend=Byend,Bzend=Bzend,Bend=Bend,unlost=unlost,lost=lost,nounlost=nounlost,countlostparticles=countlostparticles

ratio_totalenergy=ekinend/ekinstart
ratio_parenergy=eparend/eparstart
ratio_perpenergy=eperpend/eperpstart

ratio_b=bend/bstart

;normalisation coeffs.
L=10.0d6
T0=100.0d0
Vscl=L/T0
Q = - 1.6022d-19
M = 9.1095d-31 
EkinScl=-1/Q*M*Vscl*Vscl
B0=0.01d0
E0 = Vscl * B0

xstart=TEMPORARY(xstart)*L
ystart=TEMPORARY(ystart)*L
zstart=TEMPORARY(zstart)*L
tstart=TEMPORARY(tstart)*T0
ekinstart=TEMPORARY(ekinstart)*EkinScl
alphastart=TEMPORARY(alphastart)*!RaDeg
eperpstart=TEMPORARY(eperpstart)*EkinScl
eparstart=TEMPORARY(eparstart)*EkinScl
vtotstart=TEMPORARY(vtotstart)*Vscl
vperpstart=TEMPORARY(vperpstart)*Vscl
vparstart=TEMPORARY(vparstart)*Vscl
Exstart=TEMPORARY(Exstart)*E0
Eystart=TEMPORARY(Eystart)*E0
Ezstart=TEMPORARY(Ezstart)*E0
Bxstart=TEMPORARY(Bxstart)*B0
Bystart=TEMPORARY(Bystart)*B0
Bzstart=TEMPORARY(Bzstart)*B0
BStart=TEMPORARY(BStart)*B0
mu=TEMPORARY(mu)*M*Vscl^2/B0
xend=TEMPORARY(xend)*L
yend=TEMPORARY(yend)*L
zend=TEMPORARY(zend)*L
tend=TEMPORARY(tend)*T0
ekinend=TEMPORARY(ekinend)*EkinScl
alphaend=TEMPORARY(alphaend)*!RaDeg
eperpend=TEMPORARY(eperpend)*EkinScl
eparend=TEMPORARY(eparend)*EkinScl
vtotend=TEMPORARY(vtotend)*Vscl
vperpend=TEMPORARY(vperpend)*Vscl
vparend=TEMPORARY(vparend)*Vscl
Exend=TEMPORARY(Exend)*E0
Eyend=TEMPORARY(Eyend)*E0
Ezend=TEMPORARY(Ezend)*E0
Bxend=TEMPORARY(Bxend)*B0
Byend=TEMPORARY(Byend)*B0
Bzend=TEMPORARY(Bzend)*B0
Bend=TEMPORARY(Bend)*B0

lost=where(tend lt max(tend))
unlost=where(tend eq max(tend))

da={xstart:xstart,ystart:ystart,zstart:zstart,tstart:tstart,ekinstart:ekinstart,alphastart:alphastart, eperpstart:eperpstart,$
eparstart:eparstart, vtotstart:vtotstart, vperpstart:vperpstart, vparstart:vparstart, Exstart:Exstart, Eystart:Eystart, Ezstart:Ezstart,$
Bxstart:Bxstart, Bystart:Bystart, Bzstart:Bzstart,BStart:BStart, mu:mu, xend:xend, yend:yend, zend:zend, tend:tend, ekinend:ekinend, alphaend:alphaend,$
eperpend:eperpend, eparend:eparend, vtotend:vtotend, vparend:vparend, vperpend:vperpend, Exend:Exend, Eyend:Eyend, Ezend:Ezend,$
Bxend:Bxend, Byend:Byend, Bzend:Bzend, Bend:Bend, lost:lost, unlost:unlost}

return,da

end
