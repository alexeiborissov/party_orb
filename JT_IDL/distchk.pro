
;nom="./Data/"
nom="../quicksurveyP/Data/"
ds=getrdata(10, wkdir=nom+'DataR/')

nt=n_elements(ds.t)
r=dblarr(nt-1)

for i=1,nt-1 do begin

r(i-1)=sqrt((ds.x(i)-ds.x(i-1))*(ds.x(i)-ds.x(i-1))+(ds.y(i)-ds.y(i-1))*(ds.y(i)-ds.y(i-1))+(ds.z(i)-ds.z(i-1))*(ds.z(i)-ds.z(i-1)))

endfor



END
