;program to only plot the start and end points of particle traces at different times, TOGETHER WITH FAN PLANES, NULLS AND SPINES.

loadct, 4, /silent
tvlct, r, g, b, /get

out=0

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

restore, filename=dir6+'/saves/kestore.sav', /verbose

PRINT, max(kestore)

STOP
restore, filename=dir1+'/saves/kestore.sav', /verbose
myEN=inike
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
nmyEN = n_elements(myEN)
hi1 = histogram(l_EN,binsize=bsz, locations=xbin_i1)
n_hi1 = n_elements(hi1)

dEi1 = fltarr((size(xbin_i1))(1))
for j = 1,(size(xbin_i1))(1)-1 do begin
	dEi1(j) = 10^(xbin_i1(j)) - 10^(xbin_i1(j-1))
endfor
dEi1(0) = 1.0
tempx=xbin_i1
tempy=alog10(hi1/dEi1)
ntempx = N_ELEMENTS(tempx)
histo_xi1 = CONGRID(tempx, ntempx*2)
histo_yi1 = CONGRID(tempy, ntempx*2, /CENTER)
myEN=finke
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
nmyEN = n_elements(myEN)
hf1 = histogram(l_EN,binsize=bsz, locations=xbin_f1)
n_hf1 = n_elements(hf1)
dEf1 = fltarr((size(xbin_f1))(1))
for j = 1,(size(xbin_f1))(1)-1 do begin
	dEf1(j) = 10^(xbin_f1(j)) - 10^(xbin_f1(j-1))
endfor
dEf1(0) = 1.0
tempx=xbin_f1
tempy=alog10(hf1/dEf1)
ntempx = N_ELEMENTS(tempx)
histo_xf1 = CONGRID(tempx, ntempx*2)
histo_yf1 = CONGRID(tempy, ntempx*2, /CENTER)
;---;
restore, filename=dir2+'/saves/kestore.sav', /verbose
myEN=finke
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
nmyEN = n_elements(myEN)
hf2 = histogram(l_EN,binsize=bsz, locations=xbin_f2)
n_hf2 = n_elements(hf2)
dEf2 = fltarr((size(xbin_f2))(1))
for j = 1,(size(xbin_f2))(1)-1 do begin
	dEf2(j) = 10^(xbin_f2(j)) - 10^(xbin_f2(j-1))
endfor
dEf2(0) = 1.0
tempx=xbin_f2
tempy=alog10(hf2/dEf2)
ntempx = N_ELEMENTS(tempx)
histo_xf2 = CONGRID(tempx, ntempx*2)
histo_yf2 = CONGRID(tempy, ntempx*2, /CENTER)
;---;
restore, filename=dir3+'/saves/kestore.sav', /verbose
myEN=finke
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
nmyEN = n_elements(myEN)
hf3 = histogram(l_EN,binsize=bsz, locations=xbin_f3)
n_hf3 = n_elements(hf3)
dEf3 = fltarr((size(xbin_f3))(1))
for j = 1,(size(xbin_f3))(1)-1 do begin
	dEf3(j) = 10^(xbin_f3(j)) - 10^(xbin_f3(j-1))
endfor
dEf3(0) = 1.0
tempx=xbin_f3
tempy=alog10(hf3/dEf3)
ntempx = N_ELEMENTS(tempx)
histo_xf3 = CONGRID(tempx, ntempx*2)
histo_yf3 = CONGRID(tempy, ntempx*2, /CENTER)
;---;
restore, filename=dir4+'/saves/kestore.sav', /verbose
myEN=finke
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
nmyEN = n_elements(myEN)
hf4 = histogram(l_EN,binsize=bsz, locations=xbin_f4)
n_hf4 = n_elements(hf2)
dEf4 = fltarr((size(xbin_f4))(1))
for j = 1,(size(xbin_f4))(1)-1 do begin
	dEf4(j) = 10^(xbin_f4(j)) - 10^(xbin_f4(j-1))
endfor
dEf4(0) = 1.0
tempx=xbin_f4
tempy=alog10(hf4/dEf4)
ntempx = N_ELEMENTS(tempx)
histo_xf4 = CONGRID(tempx, ntempx*2)
histo_yf4 = CONGRID(tempy, ntempx*2, /CENTER)
;---;
restore, filename=dir5+'/saves/kestore.sav', /verbose
myEN=finke
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
nmyEN = n_elements(myEN)
hf5 = histogram(l_EN,binsize=bsz, locations=xbin_f5)
n_hf5 = n_elements(hf5)
dEf5 = fltarr((size(xbin_f5))(1))
for j = 1,(size(xbin_f5))(1)-1 do begin
	dEf5(j) = 10^(xbin_f5(j)) - 10^(xbin_f5(j-1))
endfor
dEf5(0) = 1.0
tempx=xbin_f5
tempy=alog10(hf5/dEf5)
ntempx = N_ELEMENTS(tempx)
histo_xf5 = CONGRID(tempx, ntempx*2)
histo_yf5 = CONGRID(tempy, ntempx*2, /CENTER)

;---;
restore, filename=dir6+'/saves/kestore.sav', /verbose
myEN=finke
myEN=myEN(sort(myEN))
ii = where(myEN gt minEN)
myEN = myEN(ii)
l_EN = alog10(myEN) 
nmyEN = n_elements(myEN)
hf6 = histogram(l_EN,binsize=bsz, locations=xbin_f6)
n_hf6= n_elements(hf6)
dEf6 = fltarr((size(xbin_f6))(1))
for j = 1,(size(xbin_f6))(1)-1 do begin
	dEf6(j) = 10^(xbin_f6(j)) - 10^(xbin_f6(j-1))
endfor
dEf6(0) = 1.0
tempx=xbin_f6
tempy=alog10(hf6/dEf6)
ntempx = N_ELEMENTS(tempx)
histo_xf6 = CONGRID(tempx, ntempx*2)
histo_yf6 = CONGRID(tempy, ntempx*2, /CENTER)
;---;

;---PLOTTING---;
aspect_ratio=1.3
myxs=20
IF NOT(OUT) THEN myxs=myxs*50
myys=myxs/aspect_ratio
loadct, 73

TVLCT, r,g,b, /get
r(0)=0
g(0)=0
b(0)=0
r(255)=255
g(255)=255
b(255)=255
TVLCT, r, g, b
 
IF (out) THEN BEGIN	; set up plot area
  set_plot,'ps'
  !p.font=0
  device, filename=od+'dN_tdep.eps', encapsulated=1, /helvetica, /color, BITS_PER_PIXEL=8
  device, xsize=myxs, ysize=myys
  mycharsize=1
  myth=4
  PRINT, 'output to: '+od+'dN_tdep.eps'
 ENDIF ELSE BEGIN
  window, 5, ysize=myys, xsize=myxs
  mycharsize=1
  myth=2
 ENDELSE 

labels=['t=0','t=90','t=115','t=190','t=235','t=270']
jtcols=[50,90,130,170,210,240]

plot, histo_xi1, histo_yi1, xtitle='energy 10^ (eV)', ytitle='dN', xr=[-1,7], yr=[-8,2], linestyle=1, thick=myth, xthick=mythick, ythick=mythick, charthick=myth-1
oplot, histo_xf1, histo_yf1, col=jtcols[0], thick=myth
oplot, histo_xf2, histo_yf2, col=jtcols[1], thick=myth
oplot, histo_xf3, histo_yf3, col=jtcols[2], thick=myth
oplot, histo_xf4, histo_yf4, col=jtcols[3], thick=myth
oplot, histo_xf5, histo_yf5, col=jtcols[4], thick=myth
oplot, histo_xf6, histo_yf6, col=jtcols[5], thick=myth


oplot, [alog10(t/11604.505),alog10(t/11604.505)], [-10,10], linestyle=5
xyouts, 10, 2, "KE (T=1E6K)", /DATA, charsize=0.9


legend, labels, linestyle=[0,0,0,0,0,0], color=jtcols, /top, /right, charsize=0.9, thick=myth, charthick=myth-1


 IF (out) THEN BEGIN
  device, /close
  set_plot, 'x'
  !p.font=-1
 ENDIF

END
