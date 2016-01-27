;works with temporary output files, using gettempdata
loadct, 39

items = ['+ve charge','-ve charge']
;       sym = [1,2,3]
;       legend,items,psym=sym
;attempting to diagnose particle code problems using temporary output.
FOR i=1,48 do begin

e=gettempdata(i,wkdir="Data/DataR/", sub='e')
p=gettempdata(i,wkdir="Data/DataR/", sub='p')


PRINT, 'max(e.t)=',max(e.t), ', max(p.t)=',max(p.t)


ENDFOR


STOP

pmodb=sqrt(p.b[0,*]*p.b[0,*]+p.b[1,*]*p.b[1,*]+p.b[2,*]*p.b[2,*])
emodb=sqrt(e.b[0,*]*e.b[0,*]+e.b[1,*]*e.b[1,*]+e.b[2,*]*e.b[2,*])
pmode=sqrt(p.e[0,*]*p.e[0,*]+p.e[1,*]*p.e[1,*]+p.e[2,*]*p.e[2,*])
emode=sqrt(e.e[0,*]*e.e[0,*]+e.e[1,*]*e.e[1,*]+e.e[2,*]*e.e[2,*])
pmodrdot=sqrt(p.drdt[0,*]*p.drdt[0,*]+p.drdt[1,*]*p.drdt[1,*]+p.drdt[2,*]*p.drdt[2,*])
emodrdot=sqrt(e.drdt[0,*]*e.drdt[0,*]+e.drdt[1,*]*e.drdt[1,*]+e.drdt[2,*]*e.drdt[2,*])
pmodbxderivs=sqrt(p.dbdx[0,*]*p.dbdx[0,*]+p.dbdy[0,*]*p.dbdy[0,*]+p.dbdz[0,*]*p.dbdz[0,*])
pmodbyderivs=sqrt(p.dbdx[1,*]*p.dbdx[1,*]+p.dbdy[1,*]*p.dbdy[1,*]+p.dbdz[1,*]*p.dbdz[1,*])
pmodbzderivs=sqrt(p.dbdx[2,*]*p.dbdx[2,*]+p.dbdy[2,*]*p.dbdy[2,*]+p.dbdz[2,*]*p.dbdz[2,*])
emodbxderivs=sqrt(e.dbdx[0,*]*e.dbdx[0,*]+e.dbdy[0,*]*e.dbdy[0,*]+e.dbdz[0,*]*e.dbdz[0,*])
emodbyderivs=sqrt(e.dbdx[1,*]*e.dbdx[1,*]+e.dbdy[1,*]*e.dbdy[1,*]+e.dbdz[1,*]*e.dbdz[1,*])
emodbzderivs=sqrt(e.dbdx[2,*]*e.dbdx[2,*]+e.dbdy[2,*]*e.dbdy[2,*]+e.dbdz[2,*]*e.dbdz[2,*])
pmodexderivs=sqrt(p.dedx[0,*]*p.dedx[0,*]+p.dedy[0,*]*p.dedy[0,*]+p.dedz[0,*]*p.dedz[0,*])
pmodeyderivs=sqrt(p.dedx[1,*]*p.dedx[1,*]+p.dedy[1,*]*p.dedy[1,*]+p.dedz[1,*]*p.dedz[1,*])
pmodezderivs=sqrt(p.dedx[2,*]*p.dedx[2,*]+p.dedy[2,*]*p.dedy[2,*]+p.dedz[2,*]*p.dedz[2,*])
emodexderivs=sqrt(e.dedx[0,*]*e.dedx[0,*]+e.dedy[0,*]*e.dedy[0,*]+e.dedz[0,*]*e.dedz[0,*])
emodeyderivs=sqrt(e.dedx[1,*]*e.dedx[1,*]+e.dedy[1,*]*e.dedy[1,*]+e.dedz[1,*]*e.dedz[1,*])
emodezderivs=sqrt(e.dedx[2,*]*e.dedx[2,*]+e.dedy[2,*]*e.dedy[2,*]+e.dedz[2,*]*e.dedz[2,*])


BEfields=0
DIFFEQNS=0
Bfieldderivs=0
moreBfields=0
moreBfieldderivs=0
Efieldderivs=1
EparandT=0
R=0
Uparandgamma=0
dmodbds=0
rg=0

myxr=[0,200]

cs=2
i=0
IF BEfields THEN BEGIN
;first, lets plot changes in |B| and |E|
 window, i, ysize=650
 !p.multi=[0,1,2]
 myyr=[-5e-4,5e-4]
 mytit=string(pmodb[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emodb[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodb-pmodb[0], xmargin=[14,8], xr=myxr, psym=-1, yr=myyr, title='delta |B|', ytitle='|B|-|B|(t=0)', xtitle='NSTP'
 oplot, e.nstp,emodb-emodb[0], col=240, psym=-4
 legend, items, psym=[-1,-4], col=[255,240]
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 mytit=string(pmode[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emode[0],format='("(e ini val=",E10.2,")")')
 myyr=[-6e-5,6e-5]
 plot, p.nstp, pmode-pmode[0], xmargin=[14,8], xr=myxr, psym=-1, yr=myyr, title='delta |E|', ytitle='|E|-|E|(t=0)', xtitle='NSTP'
 oplot, e.nstp,emode-emode[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 
!p.multi=[0]
i=i+1
ENDIF

IF DIFFEQNS THEN BEGIN
 ; now changes in differential equation variables, U, gamma and R
 window, i, ysize=1000
 !p.multi=[0,1,3]
 myyr=[-500,200]
 mytit=string(p.dudt[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dudt[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dudt-p.dudt[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dudt', xtitle='NSTP', ytitle='dudt-dudt(t=0)', charsize=cs
 oplot, e.nstp,e.dudt-e.dudt[0], col=240, psym=-4
 legend, items, psym=[-1,-4], col=[255,240]
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[0,0.04]
 mytit=string(p.dgammadt[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dgammadt[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dgammadt-p.dgammadt[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dgammadt', xtitle='NSTP', ytitle='dgammadt-dgammadt(t=0)', charsize=cs
 oplot, e.nstp,e.dgammadt-e.dgammadt[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[0,0.05]
 mytit=string(pmodrdot[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emodrdot[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodrdot-pmodrdot[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta |dRdt|', xtitle='NSTP', ytitle='|dRdt|-|dRdt|(t=0)', charsize=cs
 oplot, e.nstp,emodrdot-emodrdot[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 !p.multi=0
 i=i+1
ENDIF

IF Bfieldderivs THEN BEGIN
; changes in magnetic field derivatives?
 window, i, ysize=1000
 !p.multi=[0,1,3]
 myyr=[-0.1,0.1]
 mytit=string(pmodbxderivs[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emodbxderivs[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodbxderivs-pmodbxderivs[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbxds', xtitle='NSTP', ytitle='dbxds-dbxds(t=0)', charsize=cs
 oplot, e.nstp,emodbxderivs-emodbxderivs[0], col=240, psym=-4
 legend, items, psym=[-1,-4], col=[255,240]
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.1,0.1]
 mytit=string(pmodbyderivs[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emodbyderivs[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodbyderivs-pmodbyderivs[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbyds', xtitle='NSTP', ytitle='dbyds-dbyds(t=0)', charsize=cs
 oplot, e.nstp,emodbyderivs-emodbyderivs[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.1,0.1]
 mytit=string(pmodbzderivs[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(pmodbzderivs[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodbzderivs-pmodbzderivs[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbzds', xtitle='NSTP', ytitle='dbzds-dbzds(t=0)', charsize=cs
 oplot, e.nstp,pmodbzderivs-pmodbzderivs[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 !p.multi=0
 i=i+1
ENDIF

IF moreBfields THEN BEGIN
;first, lets plot changes in |B| and |E|
 window, i, ysize=900, xsize=400
 !p.multi=[0,1,3]
 myyr=[-0.02,0.02]
 mytit=string(p.b[0,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.b[0,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.b[0,*]-p.b[0,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta bx', xtitle='NSTP', ytitle='bx-bx(t=0)', charsize=cs
 oplot, e.nstp,e.b[0,*]-e.b[0,0], col=240, psym=-4
 legend, items, psym=[-1,-4], col=[255,240]
  xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.02,0.02]
 mytit=string(p.b[1,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.b[1,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.b[1,*]-p.b[1,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta by', xtitle='NSTP', ytitle='by-by(t=0)', charsize=cs
 oplot, e.nstp,e.b[1,*]-e.b[1,0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.02,0.02]
 mytit=string(p.b[0,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.b[0,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.b[2,*]-p.b[2,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta bz', xtitle='NSTP', ytitle='bz-bz(t=0)', charsize=cs
 oplot, e.nstp,e.b[2,*]-e.b[2,0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
!p.multi=[0]
i=i+1
ENDIF

IF moreBfieldderivs THEN BEGIN
; changes in magnetic field derivatives?
 window, i, ysize=900, xsize=1100
 !p.multi=[0,3,3]
 myyr=[-0.1,0.1]
 mytit=string(p.dbdx[0,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdx[0,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdx[0,*]-p.dbdx[0,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbxdx', xtitle='NSTP', ytitle='dbxdx-dbxdx(t=0)', charsize=cs
 oplot, e.nstp,e.dbdx[0,*]-e.dbdx[0,0], col=240, psym=-4
 legend, items, psym=[-1,-4], col=[255,240]
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.1,0.1]
 mytit=string(p.dbdx[1,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdx[1,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdx[1,*]-p.dbdx[1,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbydx', xtitle='NSTP', ytitle='dbydx-dbydx(t=0)', charsize=cs
 oplot, e.nstp,e.dbdx[1,*]-e.dbdx[1,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.1,0.1]
 mytit=string(p.dbdx[2,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdx[2,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdx[2,*]-p.dbdx[2,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbzdx', xtitle='NSTP', ytitle='dbzdx-dbzdx(t=0)', charsize=cs
 oplot, e.nstp,e.dbdx[2,*]-e.dbdx[2,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.1,0.1]
 mytit=string(p.dbdy[0,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdy[0,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdy[0,*]-p.dbdy[0,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbxdy', xtitle='NSTP', ytitle='dbxdy-dbxdy(t=0)', charsize=cs
 oplot, e.nstp,e.dbdy[0,*]-e.dbdy[0,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 ;legend, items, psym=[-1,-4], col=[255,240]
 myyr=[-0.1,0.1]
 mytit=string(p.dbdy[1,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdy[1,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdy[1,*]-p.dbdy[1,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbydy', xtitle='NSTP', ytitle='dbydy-dbydy(t=0)', charsize=cs
 oplot, e.nstp,e.dbdy[1,*]-e.dbdy[1,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.1,0.1]
 mytit=string(p.dbdy[2,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdy[2,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdy[2,*]-p.dbdy[2,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbzdy', xtitle='NSTP', ytitle='dbzdy-dbzdy(t=0)', charsize=cs
 oplot, e.nstp,e.dbdy[2,*]-e.dbdy[2,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.02,0.02]
 mytit=string(p.dbdz[0,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdz[0,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdz[0,*]-p.dbdz[0,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbxdz', xtitle='NSTP', ytitle='dbxdz-dbxdz(t=0)', charsize=cs
 oplot, e.nstp,e.dbdz[0,*]-e.dbdz[0,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 ;legend, items, psym=[-1,-4], col=[255,240]
 myyr=[-0.02,0.02]
 mytit=string(p.dbdz[1,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdz[1,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdz[1,*]-p.dbdz[1,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbydz', xtitle='NSTP', ytitle='dbydz-dbydz(t=0)', charsize=cs
 oplot, e.nstp,e.dbdz[1,*]-e.dbdz[1,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.02,0.02]
 mytit=string(p.dbdz[2,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbdz[2,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbdz[2,*]-p.dbdz[2,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dbzdz', xtitle='NSTP', ytitle='dbzdz-dbzdz(t=0)', charsize=cs
 oplot, e.nstp,e.dbdz[2,*]-e.dbdz[2,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 !p.multi=0
 i=i+1
ENDIF


IF Efieldderivs THEN BEGIN
; changes in electric field derivatives?
 window, i, ysize=1000
 !p.multi=[0,1,3]
  myyr=[-1e-3,1e-3]
 mytit=string(pmodexderivs[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emodexderivs[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodexderivs-pmodexderivs[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dexds', xtitle='NSTP', ytitle='dexds-dexds(t=0)', charsize=cs
 oplot, e.nstp,emodexderivs-emodexderivs[0], col=240, psym=-4
  xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 legend, items, psym=[-1,-4], col=[0,240]
  mytit=string(pmodeyderivs[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emodeyderivs[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodeyderivs-pmodeyderivs[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta deyds', xtitle='NSTP', ytitle='deyds-deyds(t=0)', charsize=cs
 oplot, e.nstp,emodeyderivs-emodeyderivs[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 mytit=string(pmodezderivs[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(emodezderivs[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, pmodezderivs-pmodezderivs[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta dezds', xtitle='NSTP', ytitle='dezds-dezds(t=0)', charsize=cs
 oplot, e.nstp, emodezderivs-emodezderivs[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 !p.multi=0
i=i+1
ENDIF

IF EparandT THEN BEGIN
 WINDOW, i, ysize=600
 !p.multi=[0,1,2]
 myyr=[-3e-8,5e-8]
 mytit=string(p.epar[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.epar[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.epar-p.epar[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta Epar', xtitle='NSTP', ytitle='t-t(t=0)', charsize=cs
 oplot, e.nstp,e.epar-e.epar[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 legend, items, psym=[-1,-4], col=[255,240]
 myyr= [0,2e-4]
 mytit=string(p.t[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.t[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.t-p.t[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta t', xtitle='NSTP', ytitle='t-t(t=0)', charsize=cs
 oplot, e.nstp,e.t-e.t[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 !p.multi=0
 i=i+1
ENDIF

IF R THEN BEGIN
 ; change in position Rx/Ry/Rz
 window, i, ysize=1000, xsize=600
 !p.multi=[0,1,3]
 myyr=[-0.0005,0.001]
 mytit=string(p.R[0,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.R[0,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.R[0,*]-p.R[0,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta Rx', xtitle='NSTP', ytitle='Rx-Rx(t=0)', charsize=cs
 oplot, e.nstp, e.R[0,*]-e.R[0,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 legend, items, psym=[-1,-4], col=[255,240]
 myyr=[-0.0005,0.001]
 mytit=string(p.R[1,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.R[1,0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.R[1,*]-p.R[1,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta Ry', xtitle='NSTP', ytitle='Ry-Ry(t=0)', charsize=cs
 oplot, e.nstp, e.R[1,*]-e.R[1,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 myyr=[-0.005,0.01]
 mytit=string(p.R[2,0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.R[2,0],format='("(e ini val=",E10.2,")")') 
 plot, p.nstp, p.R[2,*]-p.R[2,0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta Rz', xtitle='NSTP', ytitle='Rz-Rz(t=0)', charsize=cs
 oplot, e.nstp, e.R[2,*]-e.R[2,0], col=240, psym=-4
 xyouts, 2./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 2./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 !p.multi=0
 i=i+1
ENDIF


IF Uparandgamma THEN BEGIN
 WINDOW, i, ysize=600
 !p.multi=[0,1,2]
 myyr=[-0.04,0.04]
 mytit=string(p.u[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.u[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.u-p.u[0], xmargin=[14,8], xr=myxr, psym=-2, yr=[-0.04,0.04], $
 title='delta U', xtitle='NSTP', ytitle='u-u(t=0)', charsize=cs
 oplot, e.nstp,e.u-e.u[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 legend, items, psym=[-1,-4], col=[255,240], /bottom, /right
 myyr=[0,1e-3]
 mytit=string(p.gamma[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.gamma[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.gamma-p.gamma[0], xmargin=[14,8], xr=myxr, psym=-2, yr=[0,1e-3], $
 title='delta gamma', xtitle='NSTP', ytitle='gamma-gamma(t=0)', charsize=cs
 oplot, e.nstp,e.gamma-e.gamma[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 !p.multi=0
 i=i+1
ENDIF

IF dmodbds THEN BEGIN
cs=1
!p.multi=0
 WINDOW, i, ysize=400
 myyr=[-0.05,0.05]
 mytit=string(p.dbds[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.dbds[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.dbds-p.dbds[0], xmargin=[14,8], xr=myxr, psym=-2, yr=[-0.05,0.05], $
 title='delta DBDS', xtitle='NSTP', ytitle='DBDS-DBDS(t=0)', charsize=cs
 oplot, e.nstp,e.dbds-e.dbds[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 legend, items, psym=[-1,-4], col=[255,240], /bottom, /right
 !p.multi=0
 i=i+1
ENDIF

IF rg THEN BEGIN
cs=1
!p.multi=0
 WINDOW, i, ysize=400
 myyr=[0,5e-6]
 mytit=string(p.rg[0],format='("(p ini val=",E10.2,")")')
 mytit2=string(e.rg[0],format='("(e ini val=",E10.2,")")')
 plot, p.nstp, p.rg-p.rg[0], xmargin=[14,8], xr=myxr, psym=-2, yr=myyr, $
 title='delta rg', xtitle='NSTP', ytitle='rg-rg(t=0)', charsize=cs
 oplot, e.nstp,e.rg-e.rg[0], col=240, psym=-4
 xyouts, 3./4.*myxr[1], myyr[0]+0.9*(myyr[1]-myyr[0]), mytit
 xyouts, 3./4.*myxr[1], myyr[0]+0.84*(myyr[1]-myyr[0]), mytit2, col=240
 legend, items, psym=[-1,-4], col=[255,240], /bottom, /right
 !p.multi=0
 i=i+1
ENDIF

END
