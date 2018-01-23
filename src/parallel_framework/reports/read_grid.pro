pro read_grid, r,costheta
	n = 0L
	OPENU, 13, 'grid', /f77_unformatted
	READU,13,n
	r = dblarr(n)
	READU,13,r
	READU,13,n
	costheta = dblarr(n)
	readu,13,costheta
	close, 13

end
