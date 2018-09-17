PRO read_4d, file,data
	OPENU, 13,file, /f77_unformatted
	n1 = 0L
	n2 = 0L
	n3 = 0L
	n4 = 0L
	readu,13, n1
	readu,13, n2
	readu,13, n3
	readu,13, n4
	data = dcomplexarr(n1,n2,n3,n4)
	readu,13, data
	close, 13

END
