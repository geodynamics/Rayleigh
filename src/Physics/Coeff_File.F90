    Type :: Coefficient_Set
        Integer :: version
        Integer :: nx
        Integer :: nconst
        Integer :: nfunc
        Integer, Allocatable :: cset(:)
        Integer, Allocatable :: fset(:)
        Real*8, Allocatable  :: x(:)
        Real*8, Allocatable  :: functions(:)
        Real*8, Allocatable  :: constants(:)

        Contains
        Procedure :: Init => Initialize_Coefficient_Set
        Procedure :: Input => Read_Coefficient_Set
    End Type Coefficient_Set

    Subroutine Initialize_Coefficient_Set(self,xin,nconst,nfunc)
        Implicit None
        Class(Coefficient_Set) :: self
        Real*8, Intent(In) :: xin(1:)
        Integer, Intent(In) :: nconst, nfunc
        self%nx = size(xin)
        Allocate(self%x(1:self%nx))
        self%x(:) = xin(:)

        Allocate(self%cset(1:nconst))
        Allocate(self%constants(1:nconst))
        self%nconst = nconst
        
        Allocate(self%fset(1:nfunc))
        Allocate(self%functions(1:self%nx, 1:nfunc))
        self%nfunc = nfunc
        
        self%cset = 0
        self%fset = 0
        self%functions = 0.0d0
        self%constants = 0.0d0
    End Subroutine Initialize_Coefficient_Set

    !Note:  Should make the read routine parallel
    Subroutine Read_Coefficient_Set(self, filename)
        Implicit None
        Character*120, Intent(In) :: filename
        Character*120 :: ref_file
        Integer :: pi_integer,nx_old
        Integer :: i, k, j

        Real*8, Allocatable :: ref_arr_old(:,:), rtmp(:), rtmp2(:)
        Real*8, Allocatable :: old_radius(:)

        ref_file = Trim(my_path)//filename

        Open(unit=15,file=ref_file,form='unformatted', status='old',access='stream')
        !Verify Endianness
        Read(15)pi_integer
        If (pi_integer .eq. 314) Then

            ! Read in constants and their 'set' flags
            Read(15) self%eqversion
            Read(15) self%cset(1:self%nconst)
            Read(15) self%fset(1:self%nfunc)
            Read(15) self%constants(1:self%nconst)
            
            ! Read the file's old grid
            Read(15) nx_old
            Allocate(f_x_old(1:nx_old,1:self%nfunc)) 
            Allocate(x_old(1:nx_old))

            Read(15)(x_old(i),i=1,nx_old)
            Do k = 1, self%nfunc
                Read(15)(f_x_old(i,k) , i=1 , nx_old)
            Enddo

            !Check to see if radius is tabulated in ascending or descending order.
            !If it is found to be in ascending order, reverse the radius and the 
            !input array of functions
            If (x_old(1) .lt. x_old(nx_old)) Then

                If (my_rank .eq. 0) Write(6,*)'Reversing Radial Indices in Equation Coefficients File!'

                Allocate(rtmp(1:nr_ref))

                rtmp = x_old
                Do i = 1, nx_old
                    x_old(i) = rtmp(nx_old-i+1)
                Enddo

                Do k = 1, self%nfunc
                    rtmp(:) = f_x_old(:,k)
                    Do i = 1, nx_old
                        f_x_old(i,k) = rtmp(nr_ref-i+1)
                    Enddo
                Enddo

                DeAllocate(rtmp)

            Endif

            Close(15)


            If (nr_ref .ne. n_r) Then
                !Interpolate onto the current radial grid if necessary
                !Note that the underlying assumption here is that same # of grid points
                ! means same grid - come back to this later for generality
                Allocate(   rtmp2(1:n_r))
                Allocate( rtmp(1:nr_ref))

                Do k = 1, n_ra_functions

                    rtmp(:) = ref_arr_old(:,k)
                    rtmp2(:) = 0.0d0

                    Call Spline_Interpolate(rtmp, old_radius, rtmp2, radius)

                    ra_functions(1:n_r,k) = rtmp2
                Enddo

                DeAllocate(rtmp,rtmp2)
            Else

                ! Bit redundant here, but may want to do filtering on ref_arr array
                ra_functions(1:n_r,1:n_ra_functions) = &
                       ref_arr_old(1:n_r,1:n_ra_functions)

                If (my_rank .eq. 0) Then
                    call stdout%print("WARNING:  nr = nr_old.  Assuming grids are the same.")
                Endif
            Endif
            DeAllocate(ref_arr_old,old_radius)
            
            ! Finally, if the logarithmic derivatives of rho, T, nu, kappa, and eta were
            ! not specified, then we compute them here.
            ! only calculate the log derivative if the function was set, otherwise there
            ! are divide by zero issues
            If ((fset(8) .eq. 0) .and. (fset(1) .eq. 1)) Then
                Call log_deriv(ra_functions(:,1), ra_functions(:,8)) ! dlnrho
            Endif
            If ((fset(9) .eq. 0) .and. (fset(8) .eq. 1)) Then
                Call log_deriv(ra_functions(:,8), ra_functions(:,9), no_log=.true.) !d2lnrho
            Endif
            If ((fset(10) .eq. 0) .and. (fset(4) .eq. 1)) Then
                Call log_deriv(ra_functions(:,4), ra_functions(:,10)) !dlnT
            Endif
            If ((fset(11) .eq. 0) .and. (fset(3) .eq. 1)) Then
                Call log_deriv(ra_functions(:,3), ra_functions(:,11)) !dlnnu
            Endif
            If ((fset(12) .eq. 0) .and. (fset(5) .eq. 1)) Then
                Call log_deriv(ra_functions(:,5), ra_functions(:,12)) !dlnkappa
            Endif
            If ((fset(13) .eq. 0) .and. (fset(7) .eq. 1)) Then
                Call log_deriv(ra_functions(:,7), ra_functions(:,13)) !dlneta
            Endif
        Else
            Write(6,*)'Error.  This file appears to be corrupt (check Endian convention).'
            Write(6,*)'Pi integer: ', pi_integer
        Endif

    End Subroutine Read_Custom_Reference_File
