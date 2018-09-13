spawn, 'ls workspace/*initial', ifiles
ni = N_ELEMENTS(ifles)
cube = dblarr(64,100,32,1)
for i = 0, ni -1 DO BEGIN
	read_4d,dat,ifiles[i], lims = lims
	rstart = lims.rstart
	rend = lims.rend
	tstart = lims.tstart
	tend = lims.tend
	cube[*,istart:iend,tstart:tend,0] = dat
endfor

end
