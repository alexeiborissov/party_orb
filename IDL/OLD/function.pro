a=0.9
b=0.9
cc=0.4
esp=1.
Lv=1.
L=1.
np=100
ymin=0.1
ymax=5
y=fltarr(np)
dy0dy=fltarr(np)
y=ymin+findgen(np)*(ymax-ymin)/(np-1)
for T=1.05,20,0.1 do begin
dY0dY=1./2.*(1.+tanh((y-Lv/L)*a))/(1.+y/(cc*T)^esp)+ $ 
1./2.*(cc*T)^esp*alog(1.+y/(cc*T)^esp)*  $
(1.-tanh((y-Lv/L)*a)^2)*a-1./2.* $
(1.-tanh((y-Lv/L)*b)^2)*b*y+1./2.-1./2.*tanh((y-Lv/L)*b)
plot,y,dy0dy
wait,0.2
endfor
end
