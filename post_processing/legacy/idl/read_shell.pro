PRO READ_SHELL, file, res
	endian_tag = 0L
	version = 0L
	nrec = 0L
	ntheta = 0L
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
	READU,13, ntheta
	READU,13, nr
	READU,13, nq
	nphi       = 2*ntheta


	; Read the header arrays
	qvals      = LONARR(nq)
	radius     = DBLARR(nr)
	shell_inds = LONARR(nr)
	costheta   = DBLARR(ntheta)

	READU,13,qvals
	READU,13,radius
	READU,13,shell_inds
	READU,13,costheta


	; Read the individual records
	vals = DBLARR(nphi,ntheta,nr,nq,nrec)
	tmp  = DBLARR(nphi,ntheta)
	time = DBLARR(nrec)
	iter = LONARR(nrec)
	it = 0L
	tm = 0.0d0
	FOR i = 0, nrec-1 DO BEGIN
		FOR k = 0, nq -1 DO BEGIN
			FOR j = 0, nr -1 DO BEGIN
				READU,13,tmp
				vals[*,*,j,k,i] = tmp
			ENDFOR
		ENDFOR
		READU,13,tm
		READU,13,it
		time[i] = tm
		iter[i] = it
	ENDFOR
	CLOSE, 13

	; Build a lookup table for the quantity codes
        qmax = 901L
        lut = LONARR(qmax+1)
        lut[*] = qmax*2
        FOR i =	0, nq -1 DO BEGIN
                lut[qvals[i]] =	 i
        ENDFOR

	
	res = {vals:vals, radius:radius, costheta:costheta,qvals:qvals, shell_inds:shell_inds, time:time, $
			iter:iter, version:version, lut:lut}

END
