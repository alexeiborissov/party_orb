
.r read_RV
.r localmax
xmaxs=localmax(x)
oplot,x[xmaxs],y[xmaxs],psym=2,color=getcolor('red',2)
vparint=fltarr(n_elements(t))
for i=1,n_elements(t)-1 do vparint[i]=(t[i]-t[i-1])*(vpar[i]^2+vpar[i-1]^2)/2d0
xmaxs=[0,xmaxs]
j=0.0d0
for i=1,n_elements(xmaxs)-1 do begin $
  j=[j,total(vparint[xmaxs[i-1]:xmaxs[i]])]
xmaxs=xmaxs[1:n_elements(xmaxs)-1]
j=j[1:n_elements(xmaxs)-1]
plot,t[xmaxs],j
print,max(j),min(j)

print,max(j)-min(j)

plot,t[xmaxs],j,xtitle='time - bounce points only',ytitle='J'
write_png,"j.png",tvrd()
