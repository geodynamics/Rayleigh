PRO READ_TIME, file, res, swap = swap
	if (keyword_set(swap)) THEN BEGIN
		openu,13,file, /swap_if_little_endian
	ENDIF ELSE BEGIN
		openu,13,file
	ENDELSE
	ncol = 0L
	nrow = 0L
	ntimers = 0L
	nr = 0L
	lmax = 0L
	niter = 0L
	READU,13, ncol
	READU,13, nrow
	READU,13, ntimers
	READU,13, nr
	READU,13, lmax
	READU,13,niter
	np = nrow*ncol
	col_rank = LONARR(np)
	row_rank = LONARR(np)
	alltimes = DBLARR(ntimers,np)
	READU,13,col_rank
	READU,13,row_rank
	READU,13,alltimes

	close, 13

	res = {alltimes:alltimes, col_rank:col_rank, row_rank:row_rank, $
		np:np, nrow:nrow,ncol:ncol,lmax:lmax,nr:nr,niter:niter}

End

