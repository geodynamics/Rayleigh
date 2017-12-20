PRO READ_SHELLAVG, file, res
	endian_tag = 0L
	version = 0L
	nrec = 0L
	nr = 0L
	nq = 0L
	CLOSE, 13
	OPENR, 13, file
	READU,13, endian_tag
	If (endian_tag ne 314) THEN BEGIN
		endian1 = endian_tag
		CLOSE,13
		OPENR,13,file, /swap_if_little_endian
		READU,13,endian_tag
		IF (endian_tag ne 314) THEN BEGIN
			print, 'Unable to discern endianess of file!'
			print, "Expected integer value 314 in first 4 bytes"
			print, "Found : ", endian1
			print, "And   : ", endian_tag
			STOP
		ENDIF

	Endif
	; Read the data dimensions
	READU,13, version
	READU,13, nrec
	READU,13, nr
	READU,13, nq


	; Read the header arrays
	qvals      = LONARR(nq)
	radius     = DBLARR(nr)

	READU,13,qvals
	READU,13,radius


	; Read the individual records
	IF (version eq 1) THEN BEGIN
		vals = DBLARR(nr,nq,nrec)
		nmom = 1
	ENDIF
	IF (version ge 2) THEN BEGIN
		vals = DBLARR(nr,4,nq,nrec)
		nmom = 4
	ENDIF
	tmp  = DBLARR(nr)
	time = DBLARR(nrec)
	iter = LONARR(nrec)
	it = 0L
	tm = 0.0d0
	If (version eq 1) THEN BEGIN
	FOR i = 0, nrec-1 DO BEGIN
		FOR k = 0, nq -1 DO BEGIN
			READU,13,tmp
			vals[*,k,i] = tmp
		ENDFOR
		READU,13,tm
		READU,13,it
		time[i] = tm
		iter[i] = it
	ENDFOR
	ENDIF ELSE BEGIN
	FOR i = 0, nrec-1 DO BEGIN
		FOR k = 0, nq -1 DO BEGIN
			For m = 0, nmom-1 DO BEGIN
			READU,13,tmp
			vals[*,m,k,i] = tmp
			ENDFOR
		ENDFOR
		READU,13,tm
		READU,13,it
		time[i] = tm
		iter[i] = it
	ENDFOR

	ENDELSE

	CLOSE, 13

	; Build a lookup table for the quantity codes
	qmax = 901L
	lut = LONARR(qmax+1)
	lut[*] = qmax*2
	FOR i = 0, nq -1 DO BEGIN
		lut[qvals[i]] =  i
	ENDFOR
	
	res = {vals:vals, radius:radius, qvals:qvals, time:time, $
			iter:iter, version:version, lut:lut}

END
