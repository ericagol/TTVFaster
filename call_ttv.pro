pro call_ttv,jmax0
; This routine gives an example of a call of compute_ttv.pro
; which computes the first-order eccentricity TTVs, from
; Agol & Deck (2015).
; It uses parameters appropriate for the outer two planets
; of Kepler-62e/f.
common order,jmax
; Set jmax:
jmax=jmax0
; Read in parameters from file:
data=dblarr(10) & openr,1,'kepler62ef_planets.txt' & readf,1,data & close,1
; Compute 40 transits of inner and 20 of outer planet:
nt1=40 & nt2=20
; Set up uniform ephemeris to compute TTVs:
time1=data[2]+dindgen(nt1)*data[1]
time2=data[7]+dindgen(nt2)*data[6]
; Rearrange planet parameters for ingest to compute_ttv.pro:
param=data
param[1]=data[2]
param[2]=data[1]
param[6]=data[7]
param[7]=data[6]
; Call compute_ttv.pro:
compute_ttv,[time1,time2],param,model
; Model contains the transit times.
; Plot the results:
!p.multi=[0,1,2]
; Subtract the transit times to only plot the TTVs:
plot,time1,model[0:nt1-1]-time1,xr=[-100,5500],xs=1
plot,time2,model[nt1:*]-time2,xr=[-100,5500],xs=1
char=get_kbrd(1)
return
end
