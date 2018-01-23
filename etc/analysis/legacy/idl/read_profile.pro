pro read_profile, file, res
	CLOSE, 13

	openr, 13, file
	endian_tag = 0L
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



	nr = 0L
	n2 = 0L
	readu,13,nr
	readu,13,n2
	nq = n2-1
	vals = dblarr(nr,nq)


	radius = dblarr(nr)
	readu,13, radius
	readu,13, vals

	close, 13
	res = {radius:radius, vals:vals}

end
