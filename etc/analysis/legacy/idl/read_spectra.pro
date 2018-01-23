PRO READ_SPECTRA, file, res
	endian_tag = 0L
	version = 0L
	nrec = 0L
	lmax = 0L
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
	READU,13, lmax
	READU,13, nr
	READU,13, nq


	; Read the header arrays
	qvals      = LONARR(nq)
	radius     = DBLARR(nr)
	shell_inds = LONARR(nr)

	READU,13,qvals
	READU,13,radius
	READU,13,shell_inds


	; Read the individual records
	vals = DCOMPLEXARR(lmax+1,lmax+1,nr,nq,nrec)
	tmp  = DBLARR(lmax+1,lmax+1)
	time = DBLARR(nrec)
	iter = LONARR(nrec)
	it = 0L
	tm = 0.0d0
	FOR i = 0, nrec-1 DO BEGIN
		FOR im = 0, 1 DO BEGIN
		FOR k = 0, nq -1 DO BEGIN
			FOR j = 0, nr -1 DO BEGIN
				READU,13,tmp
				if (im eq 0) THEN BEGIN
					if (j eq 1) and (k eq 0) then begin
						;if i eq 0 then begin
						;for mm = 0, 127 do print, mm, tmp[2,mm]
						;endif
					endif
					vals[*,*,j,k,i] = dcomplex(tmp,0)
				endif else begin
					vals[*,*,j,k,i] = vals[*,*,j,k,i]+dcomplex(0,tmp)
				endelse
			ENDFOR
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

	pow = Total(abs(vals)^2,2)
	ke = reform(pow[*,*,lut[1],*])+reform(pow[*,*,lut[2],*])+reform(pow[*,*,lut[3],*])
	res = {vals:vals, radius:radius, qvals:qvals, shell_inds:shell_inds, time:time, $
			iter:iter, version:version, lut:lut, pow:pow,ke:ke}

END
