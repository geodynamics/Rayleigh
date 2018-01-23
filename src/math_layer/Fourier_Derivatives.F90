Module Fourier_Derivatives

Interface d_by_dphi
	Module Procedure d_by_dphi3D
End Interface

Contains
	!/////////////////////////////////////////////
	!  d_by_dm3D 
	!  Computes derivative of arrin along first dimension
	!  Arrin is assumed to be in spectral space in dimension 1 
	!  Storage of m values is assumed to be in FFTW's r2c in place format
	!	Note that this has no spatial scale factor (d_by_dphi vs d_by_dx)	
	Subroutine d_by_dphi3D(arrin,arrout)
		Implicit None
		Integer :: i, j, k, m, ashape(1:3),ni,nj,nk
		Real*8, Intent(InOut) :: arrin(:,1:,1:)
		Real*8, Intent(InOut) :: arrout(:,1:,1:)
		ashape = shape(arrin)
		ni = ashape(2)
		nj = ashape(3)
		nk = ashape(1)
		Do j = 1, nj
			Do i = 1, ni
				Do k = 1, nk,2
					m = (k-1)/2
					arrout(k,i,j) = -m*arrin(k+1,i,j)
					arrout(k+1,i,j) = m*arrin(k,i,j)
				Enddo
			Enddo
		Enddo
	End Subroutine d_by_dphi3D



End Module Fourier_Derivatives
