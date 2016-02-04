;quick test program to read unformatted data from files which should contain
; particle code output dumps of the lare2d/3d data.

nx=128
ny=128
nz=1
nframes=2

partybx=dblarr(nx,ny,nz,nframes)
partyby=dblarr(nx,ny,nz,nframes)
partybz=dblarr(nx,ny,nz,nframes) 
partyvx=dblarr(nx,ny,nz,nframes)
partyvy=dblarr(nx,ny,nz,nframes)
partyvz=dblarr(nx,ny,nz,nframes) 

;path might need to change depending on rel or non-rel version of the code

openr, lun, 'Data/DataR/bx2idl.dat', /GET_LUN, /F77_UNFORMATTED
READU, lun, partybx
free_lun, lun
openr, lun, 'Data/DataR/by2idl.dat', /GET_LUN, /F77_UNFORMATTED
READU, lun, partyby
free_lun, lun
openr, lun, 'Data/DataR/bz2idl.dat', /GET_LUN, /F77_UNFORMATTED
READU, lun, partybz
free_lun, lun
openr, lun, 'Data/DataR/vx2idl.dat', /GET_LUN, /F77_UNFORMATTED
READU, lun, partyvx
free_lun, lun
openr, lun, 'Data/DataR/vy2idl.dat', /GET_LUN, /F77_UNFORMATTED
READU, lun, partyvy
free_lun, lun
openr, lun, 'Data/DataR/vz2idl.dat', /GET_LUN, /F77_UNFORMATTED
READU, lun, partyvz
free_lun, lun

window, 0
shade_surf, partybx


END
