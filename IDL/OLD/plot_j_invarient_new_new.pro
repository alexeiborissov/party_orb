;modified 13/1/10: Keith Grady
;plots J (calculated by integrating vpar over time) and trip lengths 

.r read_RV
.r localmax
;make array containing steps which are turning points
xmaxs=localmax(x)
oplot,x[xmaxs],y[xmaxs],psym=2,color=getcolor('red',2)
vparint=fltarr(n_elements(t))
length=fltarr(n_elements(t))

;create array with vpar contribution from each timestep
for i=1L,n_elements(t)-1 do vparint[i]=(t[i]-t[i-1])*(vpar[i]^2+vpar[i-1]^2)/2d0
for i=1L,n_elements(t)-1 do length[i]=sqrt((x[i]-x[i-1])^2+(y[i]-y[i-1])^2+(z[i]-z[i-1])^2)
;don't add in 0, as it's not a proper bounce point!
;add in 0 to turning points array
;xmaxs=[0,xmaxs]
j=0.0d0
triplengths=0.0d0
maxe=0.0d0
mine=0.0d0
maxe1=0.0d0
mine1=0.0d0
maxe2=0.0d0
mine2=0.0d0
maxe3=0.0d0
mine3=0.0d0

;add up all vpar contributions for a round trip, counting endpoints only half
for i=1,n_elements(xmaxs)-1 do begin $
  j=[j,total(vparint[xmaxs[i-1]:xmaxs[i]])-0.5*vparint[xmaxs[i-1]]+0.5*vparint[xmaxs[i]]] & $
  triplengths=[triplengths,total(length[xmaxs[i-1]:xmaxs[i]])] & $
  maxe=[maxe,max(et[xmaxs[i-1]:xmaxs[i]])] & $
  mine=[mine,min(et[xmaxs[i-1]:xmaxs[i]])] & $
  maxe1=[maxe1,max(e1[xmaxs[i-1]:xmaxs[i]])] & $
  mine1=[mine1,min(e1[xmaxs[i-1]:xmaxs[i]])] & $
  maxe2=[maxe2,max(e2[xmaxs[i-1]:xmaxs[i]])] & $
  mine2=[mine2,min(e2[xmaxs[i-1]:xmaxs[i]])] & $
  maxe3=[maxe3,max(e3[xmaxs[i-1]:xmaxs[i]])] & $
  mine3=[mine3,min(e3[xmaxs[i-1]:xmaxs[i]])] & $
  

;vparint[i] refers to the average vpar between t[i-1] and t[i]
;hence the need to subtract the first vparint and add on one at the end
help,xmaxs
;remove last element from xmaxs (just need starting time of 'round trips')
xmaxs=xmaxs[0:n_elements(xmaxs)-2]


;remove first element from j (that was put in just to keep it double etc)
help,j
j=j[1:n_elements(j)-1]
triplengths=triplengths[1:n_elements(triplengths)-1]
maxe=maxe[1:n_elements(maxe)-1]
mine=mine[1:n_elements(mine)-1]
maxe1=maxe1[1:n_elements(maxe1)-1]
mine1=mine1[1:n_elements(mine1)-1]
maxe2=maxe2[1:n_elements(maxe2)-1]
mine2=mine2[1:n_elements(mine2)-1]
maxe3=maxe3[1:n_elements(maxe3)-1]
mine3=mine3[1:n_elements(mine3)-1]
help,j,xmaxs

plot,t[xmaxs],j
print,max(j),min(j)

print,max(j)-min(j)

print,'variation of J, in percent of average value'
j_var=(max(j)-min(j))/total(j)*n_elements(j)*100
print,j_var

plot,t[xmaxs],j,xtitle='time - right side (x>0) bounce points only',ytitle='J (needs multiplied by mass and possibly tscl)',title='variation between min/max and average j: '+string(j_var,format='(f10.5)')+'%'
write_png,"j.png",tvrd()

plot,t[xmaxs],triplengths,xtitle='start time - right side (x>0) bounce points only',ytitle='length of round trip'
write_png,"length.png",tvrd()



loadct,39
plot,x,y,xtitle="x / 10Mm",ytitle="y / 10 Mm",title="Points used as start/end of integration"
oplot,x[xmaxs],y[xmaxs],color=110,psym=2

write_png,"j_points.png",tvrd(/true)


;$gimp j.png length.png j_points.png       
