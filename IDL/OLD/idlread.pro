

openr,10,'data.dat'


t = dblarr(10000)
r1 = dblarr(10000)
r2 = dblarr(10000)
r3 = dblarr(10000)
energy = dblarr(10000)

i=0

while not eof(10) do begin

readf,10,aa,bb,cc,dd,ee

t(i) = aa
r1(i) = bb
r2(i) = cc
r3(i) = dd
energy(i) = ee

i=i+1

endwhile

close,10

end






