Module Legendre_Polynomials
	! NOTE - need to convert everything here except for the last step to quad precision eventually

	! This module contains all code necessary to generate and store the Legendre polynomials,
	! their colocation points, and their associated integration weights.
	! The polynomials computed are actually the renormalized associated legendre polynomials -
	!  - meaning that they carry the spherical harmonic normalization.
	Real*16, Allocatable :: coloc(:), gl_weights(:)
	Integer :: n_theta 
	Integer :: l_max, n_m
	Integer :: m_mod = 1	! Only calculate p_lms for every m_mod'th m
	Integer, Allocatable :: m_values(:),n_l(:),n_l_even(:),n_l_odd(:)
	Logical :: parity = .true.
	Real*16 ::	PiQuad  = 3.1415926535897932384626433832795028841972q+0
	Type, Public :: even_odd_sep
		Real*8, Allocatable :: even(:)
		Real*8, Allocatable :: odd(:)
	End Type even_odd_sep
	Type, Public :: even_odd_sepi
		Integer, Allocatable :: even(:)
		Integer, Allocatable :: odd(:)
	End Type even_odd_sepi
	Type, Public :: p_lm_array 
		Real*8, Allocatable :: data(:,:)
	End Type p_lm_array

	Type, Public :: p_lm_array_quad
		Real*16, Allocatable :: data(:,:)
	End Type p_lm_array_quad

	Type(p_lm_array_quad), Allocatable :: p_lmq(:) 
    Type(p_lm_array), Allocatable :: p_lm(:), ip_lm(:)
	Type(p_lm_array), Allocatable :: p_lm_odd(:), p_lm_even(:)
	Type(p_lm_array), Allocatable :: ip_lm_odd(:), ip_lm_even(:) ! i means 'integration weights included'
	Type(even_odd_sep), Allocatable :: lvals(:)
    Type(even_odd_sepi), Allocatable :: lvalsi(:)
Contains

Subroutine Finalize_Legendre()
	Implicit None
	DeAllocate(coloc, gl_weights)
	DeAllocate(m_values,n_l)
	Call DeAllocate_Plms(depar = .true.)
End Subroutine Finalize_Legendre

Subroutine DeAllocate_Plms(depar)
	Implicit None
	Integer :: i
	Logical, Optional, Intent(In) :: depar
	If (allocated(p_lm)) Then
		Do i = 1, n_m
			If (allocated(p_lm(i)%data)) Then
				DeAllocate(p_lm(i)%data)
			Endif
		Enddo
		DeAllocate(p_lm)
	Endif
	If (present(depar)) Then
		If (depar .and. parity) Then
			Write(6,*)'DeAllocating Parity Arrays'
			Call DeAllocate_Parity_Plms()
		Endif
	Endif
End Subroutine DeAllocate_Plms

Subroutine DeAllocate_Parity_Plms()
	Implicit None
	Integer :: i,m
	If (allocated(p_lm_odd)) Then
		Do i = 1, n_m
			If (allocated(p_lm_odd(i)%data)) Then
				DeAllocate(p_lm_odd(i)%data)
			Endif
		Enddo
		DeAllocate(p_lm_odd)
	Endif
	If (allocated(p_lm_even)) Then
		Do i = 1, n_m
			If (allocated(p_lm_even(i)%data)) Then
				DeAllocate(p_lm_even(i)%data)
			Endif
		Enddo
		DeAllocate(p_lm_even)
	Endif
	If (allocated(n_l_even)) DeAllocate(n_l_even)
	If (allocated(n_l_odd)) DeAllocate(n_l_odd)
	If (allocated(lvals)) Then
		Do m = 1, n_m
			If (allocated(lvals(m)%even)) DeAllocate(lvals(m)%even)
			If (allocated(lvals(m)%odd)) DeAllocate(lvals(m)%odd)
		Enddo
		DeAllocate(lvals)
	Endif
	If (allocated(lvalsi)) Then
		Do m = 1, n_m
			If (allocated(lvalsi(m)%even)) DeAllocate(lvalsi(m)%even)
			If (allocated(lvalsi(m)%odd)) DeAllocate(lvalsi(m)%odd)
		Enddo
		DeAllocate(lvalsi)
	Endif

End Subroutine DeAllocate_Parity_Plms

Subroutine Initialize_Legendre(nt,lmax,mval,parity_in)
	Implicit None
	Real*16 :: coloc_min,coloc_max
	Logical, Intent(In) :: parity_in
	Integer, Intent(in) :: nt,lmax, mval(:)
	

	parity = parity_in
	! Set up the grid
	coloc_min = -1
	coloc_max = 1
	n_theta = nt
	l_max = lmax
	Allocate(coloc(1:n_theta))
	Allocate(gl_weights(1:n_theta))
	n_m = size(mval)
	Allocate(m_values(1:n_m))
	m_values(:) = mval(:)
	Call Find_Colocation(coloc_min, coloc_max,coloc,gl_weights,n_theta)
	Call Compute_Plms()
End Subroutine Initialize_Legendre


Subroutine Compute_Plms()
	! Subroutine Compute_Plms(m_values,n_theta, l_max)
	! We feed this a list of m_values (presumably distributed across processors)
	! And also l_max.  This is sufficient to initialize the legendre polynomials
	Implicit None
	Real*16 ::  x,tmp,factorial_ratio,amp, renorm
	Integer :: i, m, l, mv, ntmax

	n_m = size(m_values)
	
	Allocate(n_l(1:n_m))
	! Now calculate the p_lms
	! y_lm(theta) is stored as ylm(m)%(theta,el)	
	! Forget about de-aliasing right now

	Allocate(p_lm(1:n_m))
    Allocate(p_lmq(1:n_m))
	Allocate(ip_lm(1:n_m))

    ntmax = n_theta
    If (parity) Then
	    Allocate(p_lm_odd(1:n_m))
	    Allocate(p_lm_even(1:n_m))
	    Allocate(ip_lm_odd(1:n_m))
	    Allocate(ip_lm_even(1:n_m))
	    Allocate(n_l_even(1:n_m))
	    Allocate(n_l_odd(1:n_m))
	    Allocate(lvals(1:n_m))
        Allocate(lvalsi(1:n_m))
        ntmax = n_theta/2
    Endif


    ! Compute P_lm(theta) for all l's at each m
    ! Calculation is done in quad precision.  Storage is done in double.
    ! One m at a time to save memory

	Do m = 1, n_m
		n_l(m) = l_max-m_values(m)+1
		Allocate(p_lmq(m)%data(1:ntmax,m_values(m):l_max))




	! First, fill in the l = m pieces (closed form expression)
	! and the l = m+1 pieces
	factorial_ratio = 1.0q0

		mv = m_values(m)
		Call compute_factorial_ratio(mv,factorial_ratio)
		amp = ((mv+0.5q0)/(2.0q0*PiQuad))**0.5q0		
		amp = amp*factorial_ratio
		Do i = 1, ntmax
			x = coloc(i)
			tmp = 1.0q0-x*x
			
			If (mod(mv,2) .eq. 1) Then
				!odd m
				p_lmq(m)%data(i,mv) = -amp*tmp**(mv/2+0.5q0)
			Else
				!even m
				p_lmq(m)%data(i,mv) = amp*tmp**(mv/2)
			Endif
			If (mv .lt. l_max) then
				p_lmq(m)%data(i,mv+1) = p_lmq(m)%data(i,mv)*x*(2.0q0*mv+3)**0.5q0
			Endif
		Enddo


    	!General recursion for l > m+1
		mv = m_values(m)		
		Do l = mv+2, l_max		
			Do i = 1, ntmax
				x = coloc(i)
				amp = (l-1)**2-mv*mv
				amp = amp/ (4.0q0*(l-1)**2-1.0q0)
				amp = amp**0.5q0
				tmp = p_lmq(m)%data(i,l-1)*x-amp*p_lmq(m)%data(i,l-2)
				amp = (4.0q0*l*l-1.0q0)/(l*l-mv*mv)
				p_lmq(m)%data(i,l) = tmp*amp**0.5q0
			Enddo
		Enddo

	    If (parity) Then
		    Call parity_resort(m)
	    Else
            !Fill in double precision arrays
            ! Add normalization for integration
            Allocate(ip_lm(m)%data(1:ntmax,m_values(m):l_max))
            Allocate(p_lm(m)%data(m_values(m):l_max,1:ntmax))
		    mv = m_values(m)
		    Do l = mv, l_max
			    Do i = 1, ntmax
				    
				    p_lm(m)%data(l,i)  = p_lmq(m)%data(i,l)
                    renorm = 2.0q0*PiQuad*gl_weights(i)
                    tmp = p_lmq(m)%data(i,l)*renorm
				    ip_lm(m)%data(i,l) = tmp
			    Enddo
		    Enddo
        Endif
        DeAllocate(p_lmq(m)%data)

    Enddo 
    ! DeAllocate data structures that will no longer be used.
    ! "data" attributes have already been deallocated.
    DeAllocate(p_lmq)
    If (parity) Then
        DeAllocate(p_lm)
        DeAllocate(ip_lm)
    Endif    

End Subroutine Compute_Plms

Subroutine Parity_Resort(m)
	Implicit None
    Integer, Intent(In) :: m
	Integer :: l, indeven, indodd,partest, i
    Real*16 :: renorm, tmp
	! Resort the p_lms into even and odd arrays


		n_l_even(m) = 0
		n_l_odd(m) = 0
		Do l = m_values(m), l_max
			partest = l-m_values(m)
			If (Mod(partest,2) .eq. 1) Then
				n_l_odd(m) = n_l_odd(m)+1
			Else
				n_l_even(m) = n_l_even(m)+1
			Endif
		Enddo


		If (n_l_even(m) .gt. 0) Then
			Allocate(ip_lm_even(m)%data(1:n_theta/2,1:n_l_even(m)))
			Allocate(p_lm_even(m)%data(1:n_l_even(m),1:n_theta/2))
			Allocate(lvals(m)%even(1:n_l_even(m)))
			Allocate(lvalsi(m)%even(1:n_l_even(m)))
		Endif
		If (n_l_odd(m) .gt. 0) Then
			Allocate(ip_lm_odd(m)%data(1:n_theta/2,1:n_l_odd(m)))
			Allocate(p_lm_odd(m)%data(1:n_l_odd(m),1:n_theta/2))
			Allocate(lvals(m)%odd(1:n_l_odd(m)))
			Allocate(lvalsi(m)%odd(1:n_l_odd(m)))
		Endif
		indeven = 1
		indodd = 1
		Do l = m_values(m), l_max
			partest = l-m_values(m)
			If (Mod(partest,2) .eq. 1) Then
				 lvals(m)%odd(indodd) = l
				lvalsi(m)%odd(indodd) = l
                Do i = 1, n_theta/2
                    renorm = 2.0q0*PiQuad*gl_weights(i)
                    tmp = p_lmq(m)%data(i,l)*renorm
				    ip_lm_odd(m)%data(i,indodd) = tmp
				    p_lm_odd(m)%data(indodd,i) = p_lmq(m)%data(i,l)
                Enddo
				indodd = indodd +1
                
			Else
				 lvals(m)%even(indeven) = l
				lvalsi(m)%even(indeven) = l
                Do i = 1, n_theta/2
                    renorm = 2.0q0*PiQuad*gl_weights(i)
                    tmp = p_lmq(m)%data(i,l)*renorm
    				ip_lm_even(m)%data(i,indeven) = tmp
    				p_lm_even(m)%data(indeven,i) =   p_lmq(m)%data(i,l)
                Enddo
				indeven = indeven+1
			Endif
		Enddo



End Subroutine Parity_Resort
Subroutine Find_Colocation(x1,x2,abscissas, weights, order_n)
	Implicit None
	! Compute the Gauss-Legendre colocation points and discretized weights
	! appropriate for the interval x1 < x < x2
	! Based heavily on Numerical Recipes Volume 2
	! Variables have been renamed for clarity
	! Legendre polynomial calculation is accomplished in a separate subroutine
	!	to help with readability.
	Real*16, Intent(Out) :: abscissas(1:), weights(1:)
	Real*16, Intent(In) :: x1,x2
	Integer, Intent(In) :: order_n
	Real*16 :: pn, ith_root,deriv_pn, new_guess
	Real*16 :: midpoint, scaling
	Integer :: i, n_roots
	Logical :: converged
	Real*16  :: eps, del
	midpoint = 0.5q0*(x1+x2)
	scaling  = 0.5q0*(x2-x1)
	n_roots  = (order_n+1)/2	! Roots are symmetric - exploit symmetry
	
	eps = 3.0Q-14


	Do i = 1, n_roots
		ith_root = cos(PiQuad*(i-0.25q0)/(order_n+0.5q0))
		converged = .false.
		Do While (.not. converged)			
			Call nth_legendre(ith_root,order_n,pn,deriv_pn)

			new_guess = ith_root - pn/deriv_pn
			del = abs(ith_root-new_guess)
			ith_root = new_guess
			if (del .le. eps) Then
				converged = .true.	
			Endif
			abscissas(i) = midpoint-scaling*ith_root
			abscissas(order_n+1-i) = midpoint+scaling*ith_root
			weights(i) = 2.0q0*scaling/((1.0q0-ith_root*ith_root)*deriv_pn*deriv_pn)
			weights(order_n+1-i) = weights(i)
		Enddo 
	Enddo
End Subroutine Find_Colocation

Subroutine nth_legendre(x,n,pn,deriv_pn)
	! Calculates the value of the nth legendre polynomial at location x
	! Returns x, p_n(x) and d_by_dx (p_n(x))
	Implicit None
	Real*16, Intent(Out) :: pn, deriv_pn
	Real*16, Intent(In) :: x
	Integer, Intent(In) :: n
	Integer :: j
	Real*16 :: pn_minus1, pn_minus2
	pn = 1.0q0	!p0
	pn_minus1 = 0.0q0
	! Use recursion relation to build p_order_n(x)
	Do j = 1, n
		pn_minus2 = pn_minus1
		pn_minus1 = pn
		pn = ( (2.0q0*j-1.0q0)*x*pn_minus1 - (j-1.0q0)*pn_minus2 )/j			
	Enddo
	deriv_pn = n*(x*pn-pn_minus1)/(x*x-1.0q0)
End Subroutine nth_legendre


Subroutine factorial(n,f)
	Integer :: n
	Real(16), Intent(Out) :: f
	Integer :: i
	f = 1.0q0
	Do i = 1, n
		f = f*i
	Enddo	
End Subroutine factorial

Subroutine compute_factorial_ratio(m,ratio)
		! build sqrt( (2m-1)!!(2m-1)!!/(2m)!! ) stably
		Integer, Intent(In) :: m
		Real*16, Intent(Out) :: ratio
		Integer :: i
		ratio = 1.0q0
		Do	i = 1, m
			ratio = ratio*((i-0.5q0)/i)**0.5q0  !ratio = ratio*(2m-1)/(2m)
		Enddo
End Subroutine compute_factorial_ratio
End Module Legendre_Polynomials
