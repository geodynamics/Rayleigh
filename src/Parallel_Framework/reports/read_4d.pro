PRO read_4d, file,data, lims = lims
	OPENU, 13,file, /f77_unformatted
	n1 = 0L
	n2 = 0L
	n3 = 0L
	n4 = 0L
	readu,13, n1
	readu,13, n2
	readu,13, n3
	readu,13, n4
	data = dblarr(n1,n2,n3,n4)
	readu,13, data
	readu,13, n1
	readu,13, n2
	readu,13, n3
	readu,13, n4
	lims = {rstart:n1-1,rend:n2-1,tstart:n3-1,tend:n4-1}	
	close, 13

END
