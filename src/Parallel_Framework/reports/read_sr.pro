PRO read_sr, file,res
	OPENU, 13,file, /f77_unformatted
	np = 0L
	readu,13, np
	sc = lonarr(np)
	sd = lonarr(np)
	rc = lonarr(np)
	rd = lonarr(np)
	readu,13, sc
	readu,13, sd
	readu,13, rc
	readu,13, rd

	close, 13

	res = {np:np, scount:sc,sdisp:sd,rcount:rc,rdisp:rd}
END
