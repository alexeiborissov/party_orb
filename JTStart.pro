COMMON wkdirs, wkdir_global, tt, FRon, len, an
;@Start
@SDF/Start
;@ IDL/StartCFD.pro
;.r IDL/getdata.pro  	    ; reads in data
;.r IDL/getndata.pro  	    ; reads in nonrelativistic data
.r IDL/getrdata.pro  	    ; reads in relativistic data
.r IDL/getsdata.pro 	    ; reads in summary datafiles
.r IDL/iniorbplot.pro	    ; plot initial array of positions
.r IDL/legend.pro   	    ; creates quick legend
.r IDL/mk_vector.pro	    ; makes pointy 3d arrows
.r IDL/ODEINT.pro   	    ; ODE integrator
.r IDL/orb__define.pro	    ; used to define 3d orbs
.r IDL/particletrack.pro    ; for 1 particle highlights track in 3d
;.r IDL/pvels.pro    	    ; creates 3d vector plot
;.r IDL/remove.pro   	    ; removes elements from array
;.r IDL/regionzoom.pro	    ; like particletrack but zoomable
.r IDL/symbol_obj.pro	    ; ?
;.r IDL/quickbifur.pro	    ; a program to (coarsely) look for zeros in 1d
.r IDL/getenergy.pro
;.r IDL/StartLARE.pro
;.R
;.r ODEINT
.r fieldwrapping
;.r JTblue_red

q0 = 1.602176565d-19 ; elementary charge [C]
m0 = 9.10938291d-31  ; electron mass [kg]
v0 = 2.99792458d8    ; speed of light [m/s]
kb = 1.3806488d-23   ; Boltzmann's constant [J/K]
mu0 = 4.0d-7 * !dpi  ; Magnetic constant [N/A^2]
epsilon0 = 8.8541878176203899d-12 ; Electric constant [F/m]
h_planck = 6.62606957d-34 ; Planck constant [J s]
Q0 = 1.60217646d-19 ; proton charge [C]
M0 = 9.10938188d-31 ; electron mass [kg]
kb = 1.3806503d-23  ; Boltzmann's constant [J/K]
wkdir_global="Data"

.r trilinear_3D tracefield_3D 
.r destaggerB
.r xplot3dJT



