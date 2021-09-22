; IDL code for plotting the C example run of TTVFaster
pro plot_times
readcol,"output_test_ic2",planet,epoch,time,format='d,d,d'
p0 = where(planet eq 0)
Eph = linfit(epoch[p0],time[p0])
ep0 = epoch[p0]
t0 = time[p0]
ttv0 = t0-ep0*Eph[1]-Eph[0]

p1 = where(planet eq 1)
Eph = linfit(epoch[p1],time[p1])
ep1 = epoch[p1]
t1 = time[p1]
ttv1 = t1-ep1*Eph[1]-Eph[0]
ttv0 *=24.0*60.0;
ttv1 *=24.0*60.0

plot,t0,ttv0,psym=2,xtitle = 'Time in days', ytitle = 'TTV in minutes'
oplot,t1,ttv1,psym=2,color=cgcolor("red")
stop

end
