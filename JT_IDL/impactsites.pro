;program to only plot the start and end points of particle traces at different times, TOGETHER WITH FAN PLANES, NULLS AND SPINES.

loadct, 4, /silent
tvlct, r, g, b, /get

out=1

od='./figs/'

npts=16*16*16
t=1e6
bsz = 0.3
minEN=1


;FOR e=10,10 DO BEGIN
;tempdir='../../electrons_maxwellian/epsilon'
dir1='../t000_noetabkg/'
dir2='../t090_noetabkg/'
dir3='../t115_noetabkg/'
dir4='../t190_noetabkg/'
dir5='../t235_noetabkg/'
dir6='../t270_noetabkg/'

restore, filename=dir2+'/saves/kestore.sav', /verbose

lowerendz=where(endpt[2,*] lt -9e6)
upperendz=where(endpt[2,*] gt 9e6)

;STOP
posi=[0.1,0.55,0.93,0.95]
posi2=[0.1,0.07,0.93,0.47]

;---PLOTTING---;
aspect_ratio=0.75
myxs=20
IF NOT(OUT) THEN myxs=myxs*50
myys=myxs/aspect_ratio
;loadct, 73
loadct, 4, /silent

TVLCT, r,g,b, /get
r(0)=0
g(0)=0
b(0)=0
r(255)=255
g(255)=255
b(255)=255
TVLCT, r, g, b
 
ref=-2
topref=7
alke=alog10(kestore)
astore=fix((alke)/max(topref-ref)*254)
zerof=where(astore le 0)
astore(zerof)=0
t1=ref
t2=topref

;cgwindow
;cgColorbar, Divisions=4, Minor=10, Format='(e8.1)', Ticklen=-0.25, /window,range=[10^t3*ofac,10^t2*ofac], title=mytit, /ylog

 
IF (out) THEN BEGIN	; set up plot area
  set_plot,'ps'
  !p.font=0
  device, filename=od+'isites_t90.eps', encapsulated=1, /helvetica, /color, BITS_PER_PIXEL=8
  device, xsize=myxs, ysize=myys
  mycharsize=1
  myth=4
  PRINT, 'output to: '+od+'isites_t90.eps'
 ENDIF ELSE BEGIN
  window, 5, ysize=myys, xsize=myxs
  mycharsize=1
  myth=2
 ENDELSE 

;labels=['t=0','t=90','t=115','t=190','t=235','t=270']
;jtcols=[50,90,130,170,210,240]

 plot, [-2,4], [-2,2], /nodata, pos=posi, title='upper boundary impacts', xtitle='x (Mm)', ytitle='y (Mm)', thick=myth, charsize=mycharsize, charthick=myth-1, /iso
 FOR i=1,n_elements(upperendz)-1 DO BEGIN
  oplot, [endpt[0,upperendz[i]],endpt[0,upperendz[i]]]*1e-6, [endpt[1,upperendz[i]],endpt[1,upperendz[i]]]*1e-6, psym=1, symsize=2, col=astore[upperendz[i]], thick=myth
 ENDFOR

 plot, [-2,4], [-2,2], /nodata, pos=posi2, /noerase, title='lower boundary impacts', xtitle='x (Mm)', ytitle='y (Mm)', thick=myth, charsize=mycharsize, charthick=myth-1, /iso
 FOR i=1,n_elements(lowerendz)-1 DO BEGIN
  oplot, [endpt[0,lowerendz[i]],endpt[0,lowerendz[i]]]*1e-6, [endpt[1,lowerendz[i]],endpt[1,lowerendz[i]]]*1e-6, psym=1, symsize=2, col=astore[lowerendz[i]], thick=myth
 ENDFOR

;legend, labels, linestyle=[0,0,0,0,0,0], color=jtcols, /top, /right, charsize=0.9, thick=myth, charthick=myth-1


 IF (out) THEN BEGIN
  device, /close
  set_plot, 'x'
  !p.font=-1
 ENDIF

END
