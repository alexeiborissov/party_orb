
.r read_RV
.r localmax
;make array containing steps which are turning points
xmaxs=localmax(x)
oplot,x[xmaxs],y[xmaxs],psym=2,color=getcolor('red',2)
vparint=fltarr(n_elements(t))
;create array with vpar contribution from each timestep
for i=1,n_elements(t)-1 do vparint[i]=(t[i]-t[i-1])*(vpar[i]^2+vpar[i-1]^2)/2d0
;don't add in 0, as it's not a proper bounce point!
;add in 0 to turning points array
;xmaxs=[0,xmaxs]
j=0.0d0

;add up all vpar contributions for a round trip, counting endpoints only half
for i=1,n_elements(xmaxs)-1 do begin $
  j=[j,total(vparint[xmaxs[i-1]:xmaxs[i]])-0.5*vparint[xmaxs[i-1]]+0.5*vparint[xmaxs[i]+1]]
;vparint[i] refers to the average vpar between t[i-1] and t[i]
;hence the need to subtract the first vparint and add on one at the end

;remove last element from xmaxs (just need starting time of 'round trips')
xmaxs=xmaxs[0:n_elements(xmaxs)-2]


;remove first element from j (that was put in just to keep it double etc)
help,j
j=j[1:n_elements(xmaxs)-1]
help,j,xmaxs

plot,t[xmaxs],j
print,max(j),min(j)

print,max(j)-min(j)

plot,t[xmaxs],j,xtitle='time - bounce points only',ytitle='J (needs multiplied by mass and possibly tscl)'
write_png,"j.png",tvrd()


print,'variation of J, in percent of average value'
print,(max(j)-min(j))/total(j)*n_elements(j)*100

window
plot,x,y,xtitle="x / 10Mm",ytitle="y / 10 Mm",title="Points used as start/end of integration"
oplot,x[xmaxs],y[xmaxs],color=getcolor('red',2),psym=2

write_png,"j_points.png",tvrd(/true)