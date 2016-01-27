
nom="../quicksurveyP/Data/"

gstore=dblarr(2,50)

FOR i=1,50 DO BEGIN

 ds=getrdata(i,/gyro, /fields, /vpar, /ke,/rel)
 gstore[0,i-1]=max(ds.gyror)
 
 dsp=getrdata(i,/gyro, /fields, /vpar, /ke,/rel, wkdir=nom+'DataR/')
 gstore[1,i-1]=max(dsp.gyror)
ENDFOR




END
