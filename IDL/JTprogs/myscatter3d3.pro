;myscatter3d - 
; requires pvels
; sets up a plot axis, and picks which points are plotted by pvels.


@ODEINT
@mk_vector
heap_gc
common somename, tt, FRon

eps=0.0000001
h1=0.01
string1='derivs'
len=1e6
FRon=1
tt=0.0

thisPalette = obj_New('IDLgrPalette')
thisPalette->loadct, 34

oOrb = obj_new('orb', COLOR=[0, 255 ,0]) 
oOrb->Scale, .02*len, .02*len, .125*len
oSymbol = obj_new('IDLgrSymbol', oOrb) ;oSymbol is green orb for start point


;npts=16*16
nom="../Data/"
;nom="../../Efield/Data/"
;particletrack, 1, n_c=32, tscl=0.0, op=0, floc=nom, lscl=len,xyzt

xplot3d, dblarr(11),dblarr(11),findgen(11)-5, xrange=[-1,1]*len, yrange=[-1,1]*len, zrange=[-6,6]*len, $
color=[100,149,237], TITLE=string(tt,format='("ExB, t=",f5.2)')

mvals=0
count=0
jjj=[128,256,512,1024]
;jjj=[4,12,100,130,140,160,190,230,250]
mstore=dblarr(4,n_elements(jjj))
;FIRST LOOP FIGURES OUT SCALED VALUES
FOR jj=0,n_elements(jjj)-1 DO BEGIN
 pvels, jjj[jj], wkdir=nom, mvals=mvals, lscl=1e6, bscl=0.01, tscl=100.0, /silent, stride=30
 mstore(*,jj)=mvals
ENDFOR
 newmax=[max(mstore[0,*]),max(mstore[1,*]),max(mstore[2,*]),max(mstore[3,*])]

; SECOND LOOP PLOTS - NOTE THE LACK OF /SILENT SWITCH
FOR jj=0,n_elements(jjj)-1 DO BEGIN
 pvels, jjj[jj], wkdir=nom, mvals=newmax, lscl=1e6, bscl=0.01, tscl=100.0, stride=30
ENDFOR


        
END
