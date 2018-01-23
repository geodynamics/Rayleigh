Module TestSuite
	Use Test_SHT
	Use Test_Cheby
	Use ProblemSize
	Implicit None
	Logical :: test_mode = .false.
	Logical :: test_legendre = .false.
	Logical :: test_chebyshev = .false.
	Namelist /Test_Namelist/ test_mode,ntest_legendre, test_legendre, test_chebyshev
Contains
	Subroutine Test_Lib()
		If (my_rank .eq. 0) Then
			write(6,*)'Initiating library function tests.'
		Endif
		If (test_legendre) Call Test_Spherical_Transforms()
		If (test_chebyshev) Call Test_Chebyshev_Transforms()
	End Subroutine Test_Lib


End Module TestSuite
