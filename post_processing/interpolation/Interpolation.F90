Module Interpolation
    !Use IFPORT !If Intel
    Use OMP_LIB
    Implicit None
    Logical :: to_cartesian = .True.
    Logical :: to_spherical =.False.
    Integer, Parameter :: TrueDouble = Selected_real_kind(14,99)
    Integer, Public :: iteration, iteration_step, initial_iteration, final_iteration 
    Integer, Public :: nr, nth, nphi, onr
    Integer :: nfloat  ! DEAL WITH THIS!
    Integer :: nthrd=0
    Integer, Public :: ncube=0
    Real*8 :: dpi = 3.1415926535897932384626433832795028841972d0
    Real*8 :: rmin, rmax
    Real*8, Allocatable, Dimension(:) :: newx, newy, newz !Old grid
    Real*8, Allocatable, Dimension(:) :: oldtheta, oldr !Old grid
    Real*8, Allocatable, Dimension(:) :: phi, theta, radius !New grid
    Real*8, Allocatable, Dimension(:) :: sinphi(:), cosphi(:), sintheta(:), costheta(:) ! Old grid trig
    Real*8, Allocatable, Dimension(:,:,:) :: olddata, data, newdata, mag_data
    Real*8, Allocatable, Dimension(:,:,:) :: rdata, tdata, pdata
    Character(1024) :: input_file = 'None', output_file = 'None', grid_file = 'None'
    Character(1024) :: input_rfile = 'None', input_tfile = 'None', input_pfile = 'None'
    Character(1024) :: output_xfile = 'None', output_yfile = 'None', output_zfile = 'None'
    Character(1024) :: mag_file='None'
    Character(24) :: quantity
    Logical :: double_precision_output = .false. ! output
    Logical :: single_precision_input  = .false.
    Logical :: sphere_mean = .false., phi_mean = .false.
    Logical :: verbose = .false., vector_mode=.false., output_mag=.false.
    Real*8  :: rmax_zero = -1.0d0, rmin_zero = -1.0d0

Contains
  
    Subroutine Read_Grid(gfile1,gfile2)
        Implicit None
        Character(1024), Intent(In) :: gfile1, gfile2
        Integer*4 :: grid_unit, ir, it, etag, version
        Real(TrueDouble) ::  dth, dr, dphi
        Character(1024) :: igfile
        grid_unit = 11

        igfile = gfile1
        If (gfile2 .ne. 'None') igfile = gfile2

        If (verbose) Then
            Write(6,*)'Reading grid file: ', Trim(igfile)
            Write(6,*)""
        Endif

        Open(grid_unit, file=Trim(igfile), form='unformatted', access='stream',status='old')
        Read(grid_unit) etag
        If (gfile2 .eq. 'None') Read(grid_unit) version
        Read(grid_unit) onr
        Read(grid_unit) nth
        Read(grid_unit) nphi
        Allocate(oldr(onr), oldtheta(nth), phi(nphi), theta(nth), radius(nr))
        Allocate(sintheta(nth),costheta(nth),sinphi(nphi),cosphi(nphi))
        Read(grid_unit) oldr
        Read(grid_unit) oldtheta

        Close(grid_unit)

        rmin = Minval(oldr)
        rmax = Maxval(oldr)

        If (verbose) Then
            Write(6,*) 'Grid read successfully.'
            Write(6,*)''
            Write(6,*) 'Input dimensions:' 
            Write(6,*) '     nr = ', onr
            Write(6,*) ' ntheta = ', nth
            Write(6,*) '   nphi = ', nphi
            Write(6,*) '   rmin = ', rmin
            Write(6,*) '   rmax = ', rmax
            Write(6,*) ''
        Endif

        dphi = 2d0*dpi/Dble(nphi)
        dr = (rmax-rmin)/Dble(nr-1)

        Do it=0,nphi-1
            phi(it+1) = dphi*Dble(it) - dpi
            sinphi(it+1)=sin(phi(it+1))
            cosphi(it+1)=cos(phi(it+1))
        EndDo

        Do it=0,nth-1
            theta(it+1) = oldtheta(it+1) 
            sintheta(it+1)=sin(theta(it+1))
            costheta(it+1)=cos(theta(it+1))
        EndDo

        Do ir=0,nr-1
            radius(ir+1) = rmax - dr*Dble(ir)
        EndDo

        If (to_cartesian) Then
            Allocate(newx(ncube),newy(ncube),newz(ncube))
            dr = 2d0*rmax/Dble(ncube-1)
            Do it=0,ncube-1
                newx(it+1) = dr*Dble(it)-rmax
            EndDo
            newy = newx
            newz = newx
        Else
            ! Leaving this here as a reminder for spherical interpolation
            Allocate(olddata(nphi,nth,onr),data(nphi,nth,nr))
            olddata=0.0d0
            data=0.0d0
        EndIf

    End Subroutine Read_Grid

    Subroutine Interpolate()
        Implicit None

        If (verbose) Then
            Write(6,*)''
            Write(6,*)'Interpolating onto N^3 Cartesian grid with N=',ncube
            Write(6,*)''
        Endif
        
        If (vector_mode) Then

            Call Read_Grid(input_rfile, grid_file)
            Call Allocate_Data()

            Call Read_Data(rdata,input_rfile)
            Call Read_Data(tdata,input_tfile)
            Call Read_Data(pdata,input_pfile)

            Call Compute_vx(olddata,rdata,tdata,pdata)
            Call Interp3D(olddata,newdata)
            Call Write_Data(newdata,output_xfile)     

            If (output_mag) mag_data=newdata**2

            Call Compute_vy(olddata,rdata,tdata,pdata)
            Call Interp3D(olddata,newdata)
            Call Write_Data(newdata,output_yfile)

            If (output_mag) mag_data = mag_data+newdata**2

            Call Compute_vz(olddata,rdata,tdata)
            Call Interp3D(olddata,newdata)
            Call Write_Data(newdata,output_zfile)

            If (output_mag) Then
                mag_data = mag_data+newdata**2
                newdata = sqrt(mag_data)
                Call Write_Data(newdata,mag_file)
            Endif
            
        Else

            Call Read_Grid(input_file, grid_file)
            Call Allocate_Data()

            Call Read_Data(olddata,input_file)
            Call Interp3D(olddata,newdata)
            Call Write_Data(newdata,output_file)                

        Endif

    End Subroutine Interpolate

    Subroutine Allocate_Data()
        Implicit None

        Allocate(olddata(nphi,nth,onr))
        Allocate(newdata(ncube,ncube,ncube))
        olddata = 0.0d0
        newdata = 0.0d0

        If (vector_mode) Then

            Allocate(  rdata(nphi,nth,onr))
            Allocate(  tdata(nphi,nth,onr))
            Allocate(  pdata(nphi,nth,onr))

            rdata = 0.0d0
            tdata = 0.0d0
            pdata = 0.0d0
            if (output_mag) then
                Allocate(mag_data(ncube,ncube,ncube))
                mag_data = 0.0d0
            Endif
        Endif

    End Subroutine Allocate_Data

    Subroutine Write_Data(odata,ofile)
        Implicit None
        Real*8, Intent(In) :: odata(:,:,:)
        Character(1024), Intent(In) :: ofile
        Real*4, Allocatable :: buff(:,:,:)
        Integer :: interp_unit=12
        
        If (verbose) Then
            Write(6,*)' '
            Write(6,*)'Interpolation complete.'
            Write(6,*)'Output data min: ', Minval(odata)
            Write(6,*)'Output data max: ', Maxval(odata)
            Write(6,*)''
            Write(6,*)'Writing output to: ', Trim(ofile)
            Write(6,*)''
        Endif

        If (double_precision_output) Then
            Open(interp_unit, file=Trim(ofile), form='unformatted', access='stream', status='unknown')
            Write(interp_unit) odata
            Close(interp_unit)
        Else
            Allocate(buff(ncube,ncube,ncube))
            buff = Real(odata)
            Open(interp_unit, file=Trim(ofile), form='unformatted', access='stream', status='unknown')
            Write(interp_unit) buff
            Close(interp_unit)
            Deallocate(buff)
        EndIf

    End Subroutine Write_Data

    Subroutine Compute_Vx(odata,irdata,itdata,ipdata)
        Implicit None
        Real*8, Intent(Out) :: odata(:,:,:)
        Real*8, Intent(In), Dimension(:,:,:) :: irdata, itdata,  ipdata
        Integer :: i, j, p
        Do j = 1, onr
        Do i = 1, nth
        Do p = 1, nphi
            odata(p,i,j) = (irdata(p,i,j)*sintheta(i) + itdata(p,i,j)*costheta(i))*cosphi(p)
            odata(p,i,j) = odata(p,i,j)-ipdata(p,i,j)*sinphi(p)
        Enddo
        Enddo
        Enddo
    End Subroutine Compute_Vx

    Subroutine Compute_Vy(odata,irdata,itdata,ipdata)
        Implicit None
        Real*8, Intent(Out) :: odata(:,:,:)
        Real*8, Intent(In), Dimension(:,:,:) :: irdata, itdata,  ipdata
        Integer :: i, j, p
        Do j = 1, onr
        Do i = 1, nth
        Do p = 1, nphi
            odata(p,i,j) = (irdata(p,i,j)*sintheta(i) + itdata(p,i,j)*costheta(i))*sinphi(p)
            odata(p,i,j) = odata(p,i,j)+ipdata(p,i,j)*cosphi(p)
        Enddo
        Enddo
        Enddo
    End Subroutine Compute_Vy

    Subroutine Compute_Vz(odata,irdata,itdata)
        Implicit None
        Real*8, Intent(Out) :: odata(:,:,:)
        Real*8, Intent(In), Dimension(:,:,:) :: irdata, itdata
        Integer :: i, j
        Do j = 1, onr
        Do i = 1, nth
            odata(:,i,j) = irdata(:,i,j)*costheta(i) - itdata(:,i,j)*sintheta(i)
        Enddo
        Enddo
    End Subroutine Compute_Vz

    Subroutine Read_Data(indata,infile)
        Implicit None
        Real*8, Intent(InOut) :: indata(:,:,:)
        Character(1024), Intent(In) :: infile
        Integer :: interp_unit, ir, ip, chunk, i,j
        Integer*4 :: version, etag
        Real*8, Allocatable, Dimension(:,:) :: otmp, ntmp, phi_avg
        Real*8, Allocatable, Dimension(:) :: tiweights(:)
        Real*4, Allocatable, Dimension(:,:,:) :: buff
        Real*8 :: ovnphi, sphere_avg

        interp_unit=12

        If (verbose) Then
            Write(6,*)""
            Write(6,*)'Reading data from: ', Trim(infile)
            Write(6,*)""
        Endif

        Open(interp_unit, file=Trim(infile), form='unformatted', access='stream', status='unknown')
        If (grid_file .eq. 'None') Then
            ! New spherical_3D will contain grid info in header
            Read(interp_unit) etag
            Read(interp_unit) version
            Read(interp_unit) onr
            Read(interp_unit) nth
            Read(interp_unit) nphi
            Read(interp_unit) oldr
            Read(interp_unit) oldtheta
        Endif

        If (single_precision_input) Then
            Allocate(buff(nphi,nth,onr))
            Read(interp_unit) buff
            indata = Dble(buff)
            Deallocate(buff)
        Else
            Read(interp_unit) indata
        EndIf
        Close(interp_unit)

        If (verbose) Then
            Write(6,*) 'Data read successfully.'
            Write(6,*)''
            Write(6,*) 'Input data min: ', minval(indata)
            Write(6,*) 'Input data max: ', maxval(indata)
        Endif

        !If desired, zero out data outside of [rmin,rmax]
        If (rmax_zero .ge. 0.0d0) Then
            If (verbose) Then
                Write(6,*)''
                Write(6,*)'Settin data to 0 for r greater than: ', rmax_zero
                Write(6,*)''
            Endif
            Do i = 1, onr
                If (oldr(i) .gt. rmax_zero) Then
                    olddata(:,:,i) = 0.0d0
                Endif
            Enddo
        Endif

        If (rmin_zero .ge. 0.0d0) Then
            If (verbose) Then
                Write(6,*)''
                Write(6,*)'Setting data to 0 for r less than: ', rmin_zero
                Write(6,*)''
            Endif
            Do i = 1, onr
                If (oldr(i) .lt. rmin_zero) Then
                    olddata(:,:,i) = 0.0d0
                Endif
            Enddo
        Endif

        !If desired, remove means
        If (sphere_mean .or. phi_mean)  Then
            Allocate(phi_avg(nth,onr))
            ovnphi=1.0d0/DBLE(nphi)
            Do j = 1, onr
                Do i = 1, nth
                    phi_avg(i,j) = Sum(olddata(:,i,j))*ovnphi
                Enddo
            Enddo

            If (phi_mean .and. (.not. sphere_mean)) Then
                If (verbose) Then
                    Write(6,*)''
                    Write(6,*)'Removing phi mean.'
                    Write(6,*)''
                Endif
                Do j = 1, onr
                    Do i = 1, nth
                        olddata(:,i,j) = olddata(:,i,j)-phi_avg(i,j)
                    Enddo
                Enddo
            Endif

            If (sphere_mean) Then
                If (verbose) Then
                    Write(6,*)''
                    Write(6,*)'Removing spherical mean.'
                    Write(6,*)''
                Endif
                Allocate(tiweights(1:nth))
                ! multiply by 0.25 rather than summing weights to normalize
                tiweights(1) = sin(oldtheta(1))*(oldtheta(1)-oldtheta(2))*0.25d0
                tiweights(nth) = sin(oldtheta(nth))*(oldtheta(nth-1)-oldtheta(nth))*0.25d0
                Do j = 2, nth-1
                    tiweights(j) = sin(oldtheta(j))*(oldtheta(j-1)-oldtheta(j+1))*0.25d0
                Enddo
                Do i = 1, onr
                    sphere_avg = SUM(phi_avg(:,i)*tiweights(:))
                    olddata(:,:,i) = olddata(:,:,i)-sphere_avg
                Enddo
                DeAllocate(tiweights)
            Endif
            DeAllocate(phi_avg)
        Endif

    End Subroutine Read_Data

    Subroutine Finalize_Interp()
        Deallocate(olddata, newdata, oldtheta, oldr, phi, theta, radius)
        If (to_cartesian) Then
            Deallocate(newx, newy, newz)
        EndIf
        If (vector_mode) Then
            DeAllocate(rdata,tdata,pdata)
        Endif
    End Subroutine Finalize_Interp


  Subroutine Interp3d(indata, outdata)
    Real*8, Intent(In) :: indata(1:,1:,1:)
    Real*8, Intent(Out) :: outdata(1:,1:,1:)
    Real*8, Dimension(:), Allocatable:: ox, oy, oz ! Mapped to Old grid
    Real*8, Dimension(:), Allocatable :: nx, ny, nz !Mapped to New grid
    Integer :: onx, ony, onz !Size of ox and oy
    Integer :: nnx, nny, nnz !Size of nx and ny
    Integer :: ix(4), iy(4), iz(4)  !Integer storage array for locations of 64 nearest-neighbors on the old grid
    Real*8 :: xc, yc, zc, rc, tc, pc !Current location on the new grid
    Real*8 :: mu, mu2, mu3 !Normalized current grid coordinate
    Real*8 :: akm1, ak, akp1, akp2 !Storage for old data
    Real*8 :: a0, a1, a2, a3 !Interpolant coefficients
    Real*8 :: ry(4,4), rz(4) !Array temporary after x interpolation
    Integer :: i1, i2, i3, ii, ij, ik, chunk, nelem !loop variables
    Integer :: ttest

    !phi,oldtheta,oldr,newx,newy,newz

    nnx = Size(newx)
    nny = Size(newy)
    nnz = Size(newz)
    
    onx = Size(phi)
    ony = Size(oldtheta)
    onz = Size(oldr)

    Allocate(ox(1:onx))
    Allocate(oy(1:ony))
    Allocate(oz(1:onz))
    ox(:) = phi(:)
    oy(:) = oldtheta(:)
    oz(:) = oldr(:)

    Allocate(nx(1:nnx))
    Allocate(ny(1:nny))
    Allocate(nz(1:nnz))
    nx(:) = newx(:)
    ny(:) = newy(:)
    nz(:) = newz(:)

    chunk = 2
    nelem = 0
    !$OMP PARALLEL DO SHARED(indata,outdata,nnx,nny,nnz,onx,ony,onz,nx,ny,nz,ox,oy,oz,chunk,rmax,rmin) &
    !$OMP PRIVATE(i1,i2,i3,ii,ij,ik,zc,yc,xc,pc,tc,rc,ix,iy,iz,mu,mu2,mu3,akm1,ak,akp1,akp2,a0,a1,a2,a3,ry,rz,ttest) &
    !$OMP SCHEDULE(DYNAMIC,chunk)
    Do i3=1,nnz
       zc = nz(i3)
       Do i2=1,nny
          yc = ny(i2)
          Do i1=1,nnx
             xc = nx(i1)

             !Convert cart to sphere
             rc = sqrt(xc*xc+yc*yc+zc*zc)
             If ((rc .gt. rmax).or.(rc .lt. rmin)) Then
                outdata(i1,i2,i3) = 0d0
                nelem=nelem+1
             Else
                tc = acos(zc/rc)
                pc = atan2(yc,xc)

                !For tricubic convolution interpolation, need nearest 64 points from old grid
                !  in each direction. Find them in O(N) time.
                ii = 1
                Do While ((ii .lt. onz) .and. (rc .lt. oz(ii))) 
                   ii = ii + 1
                EndDo
                iz(2) = Min(Max(ii-1,1),onz)
                iz(3) = Min(iz(2)+1,onz)
                iz(4) = Min(iz(2)+2,onz)
                iz(1) = Max(iz(2)-1,1)

                ii = 1
                Do While ((ii .lt. ony) .and. (tc .lt. oy(ii)))  ! Note -- had to change this due to reversed theta grid (was tc .gt. oy(ii)) Nick F
                   ii = ii + 1
                EndDo
                iy(2) = Min(Max(ii-1,1),ony)
                iy(3) = Min(iy(2)+1,ony)
                iy(4) = Min(iy(2)+2,ony)
                iy(1) = Max(iy(2)-1,1)

                ii = 1
                Do While ((ii .lt. onx) .and. (pc .gt. ox(ii))) 
                   ii = ii + 1
                EndDo
                ix(2) = Min(Max(ii-1,1),onx)
                ix(3) = Min(ix(2)+1,onx)
                ix(4) = Min(ix(2)+2,onx)
                ix(1) = Max(ix(2)-1,1)
                
                mu = (pc-ox(ix(2)))/(ox(ix(3))-ox(ix(2)))
                mu2 = mu*mu
                mu3 = mu*mu2
                Do ik=1,4
                   Do ij=1,4
                      akm1 = indata(ix(1),iy(ij),iz(ik))
                      ak   = indata(ix(2),iy(ij),iz(ik))
                      akp1 = indata(ix(3),iy(ij),iz(ik))
                      akp2 = indata(ix(4),iy(ij),iz(ik))
                      
                      a3 = ak 
                      a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
                      a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
                      a2 = -0.5d0*akm1 + 0.5d0*akp1
                      ry(ij,ik) = a0*mu3+a1*mu2+a2*mu+a3           
                   EndDo
                EndDo
                
                mu = (tc-oy(iy(2)))/(oy(iy(3))-oy(iy(2)))
                mu2 = mu*mu
                mu3 = mu*mu2             
                Do ik=1,4             
                   akm1 = ry(1,ik)
                   ak   = ry(2,ik)
                   akp1 = ry(3,ik)
                   akp2 = ry(4,ik)
                   
                   a3 = ak 
                   a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
                   a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
                   a2 = -0.5d0*akm1 + 0.5d0*akp1
                   rz(ik) = a0*mu3+a1*mu2+a2*mu+a3
                EndDo
                
                mu = (rc-oz(iz(2)))/(oz(iz(3))-oz(iz(2)))
                mu2 = mu*mu
                mu3 = mu*mu2

                akm1 = rz(1)
                ak   = rz(2)
                akp1 = rz(3)
                akp2 = rz(4)
             
                a3 = ak 
                a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
                a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
                a2 = -0.5d0*akm1 + 0.5d0*akp1
                outdata(i1,i2,i3) = a0*mu3+a1*mu2+a2*mu+a3
             EndIf
          EndDo
       EndDo
    EndDo
    !$OMP END PARALLEL DO 

    DeAllocate(ox,oy,oz)
    DeAllocate(nx,ny,nz)
    !Print*, 'Percent zeros:', dble(nelem)/dble(nnx*nny*nnz), '%'
  End Subroutine Interp3d


!////////////// Old Routines
    Subroutine Interpolate_Old()
        Implicit None
        Integer :: interp_unit, ir, ip, chunk, i,j
        Integer*4 :: reclen, etag
        Real*8, Allocatable, Dimension(:,:) :: otmp, ntmp, phi_avg
        Real*8, Allocatable, Dimension(:) :: tiweights(:)
        Real*4, Allocatable, Dimension(:,:,:) :: buff
        Real*8 :: ovnphi, sphere_avg

        interp_unit=12

        If (verbose) Then
            Write(6,*)""
            Write(6,*)'Reading input file: ', Trim(input_file)
            Write(6,*)""
        Endif
        reclen = nfloat*onr*nth*nphi

        If (single_precision_input) Then
            Allocate(buff(nphi,nth,onr))
            Open(interp_unit, file=Trim(input_file), form='unformatted', access='stream', status='unknown')
            Read(interp_unit) etag
            If (etag .eq. 314) Then
                !Call Read_Grid_New(interp_unit)
            Else
                ! a little silly, but fseek works differently between Intel and GNU
                Close(interp_unit)
                !Call Read_Grid_Legacy()
                Open(interp_unit, file=Trim(input_file), form='unformatted', access='stream', status='unknown')
            Endif
            Read(interp_unit) buff
            Close(interp_unit)
            olddata = Dble(buff)
            Deallocate(buff)
        Else
            Open(interp_unit, file=Trim(input_file), form='unformatted', access='stream', status='unknown')
            Read(interp_unit) etag
            If (etag .eq. 314) Then
                !Call Read_Grid_New(interp_unit)
            Else
                ! a little silly, but fseek works differently between Intel and GNU
                Close(interp_unit)
                !Call Read_Grid_Legacy()
                Open(interp_unit, file=Trim(input_file), form='unformatted', access='stream', status='unknown')
            Endif
            Read(interp_unit) olddata
            Close(interp_unit)
        EndIf

        If (verbose) Then
            Write(6,*) 'Data read successfully.'
            Write(6,*)''
            Write(6,*) 'Input dimensions:' 
            Write(6,*) '     nr = ', onr
            Write(6,*) ' ntheta = ', nth
            Write(6,*) '   nphi = ', nphi
            Write(6,*) '   rmin = ', rmin
            Write(6,*) '   rmax = ', rmax
            Write(6,*) ''
            Write(6,*) 'Input data min: ', minval(olddata)
            Write(6,*) 'Input data max: ', maxval(olddata)
        Endif

    !If desired, zero out data outside of [rmin,rmax]
    If (rmax_zero .ge. 0.0d0) Then
        If (verbose) Then
            Write(6,*)''
            Write(6,*)'Settin data to 0 for r greater than: ', rmax_zero
            Write(6,*)''
        Endif
        Do i = 1, onr
            If (oldr(i) .gt. rmax_zero) Then
                olddata(:,:,i) = 0.0d0
            Endif
        Enddo
    Endif

    If (rmin_zero .ge. 0.0d0) Then
        If (verbose) Then
            Write(6,*)''
            Write(6,*)'Setting data to 0 for r less than: ', rmin_zero
            Write(6,*)''
        Endif
        Do i = 1, onr
            If (oldr(i) .lt. rmin_zero) Then
                olddata(:,:,i) = 0.0d0
            Endif
        Enddo
    Endif

    !If desired, remove means
    If (sphere_mean .or. phi_mean)  Then
        Allocate(phi_avg(nth,onr))
        ovnphi=1.0d0/DBLE(nphi)
        Do j = 1, onr
            Do i = 1, nth
                phi_avg(i,j) = Sum(olddata(:,i,j))*ovnphi
            Enddo
        Enddo

        If (phi_mean .and. (.not. sphere_mean)) Then
            If (verbose) Then
                Write(6,*)''
                Write(6,*)'Removing phi mean.'
                Write(6,*)''
            Endif
            Do j = 1, onr
                Do i = 1, nth
                    olddata(:,i,j) = olddata(:,i,j)-phi_avg(i,j)
                Enddo
            Enddo
        Endif

        If (sphere_mean) Then
            If (verbose) Then
                Write(6,*)''
                Write(6,*)'Removing spherical mean.'
                Write(6,*)''
            Endif
            Allocate(tiweights(1:nth))
            ! multiply by 0.25 rather than summing weights to normalize
            tiweights(1) = sin(oldtheta(1))*(oldtheta(1)-oldtheta(2))*0.25d0
            tiweights(nth) = sin(oldtheta(nth))*(oldtheta(nth-1)-oldtheta(nth))*0.25d0
            Do j = 2, nth-1
                tiweights(j) = sin(oldtheta(j))*(oldtheta(j-1)-oldtheta(j+1))*0.25d0
            Enddo
            Do i = 1, onr
                sphere_avg = SUM(phi_avg(:,i)*tiweights(:))
                olddata(:,:,i) = olddata(:,:,i)-sphere_avg
            Enddo
            DeAllocate(tiweights)
        Endif
        DeAllocate(phi_avg)
    Endif

    If (to_cartesian) Then

        If (verbose) Then
            Write(6,*)'Interpolating on N^3 Cartesian grid with N=',ncube
            Write(6,*)''
        Endif
        Call Interp3d_Old(phi,oldtheta,oldr,newx,newy,newz)
        If (verbose) Then
            Write(6,*)'Interpolation complete.'
            Write(6,*)'Output data min: ', Minval(data)
            Write(6,*)'Output data max: ', Maxval(data)
            Write(6,*)''
            Write(6,*)'Writing output to: ', Trim(output_file)
            Write(6,*)''
        Endif
        If (double_precision_output) Then
            Open(interp_unit, file=Trim(output_file), form='unformatted', access='stream', status='unknown')
            Write(interp_unit) data
            Close(interp_unit)
        Else
            Allocate(buff(ncube,ncube,ncube))
            buff = Real(data)
            Open(interp_unit, file=Trim(output_file), form='unformatted', access='stream', status='unknown')
            !Write(interp_unit,rec=1) buff
            Write(interp_unit) buff
            Close(interp_unit)
            Deallocate(buff)
        EndIf
    Else
        chunk = 2
        !$OMP PARALLEL SHARED(olddata,data,oldtheta,oldr,theta,radius,nphi,nth,onr,nr) PRIVATE(ir,ip,otmp,ntmp)
        Allocate(otmp(nth,onr),ntmp(nth,nr))
        !$OMP DO SCHEDULE(DYNAMIC,chunk)
        Do ip=1,nphi
            otmp = olddata(ip,:,:)
            ntmp = 0d0
            Call Interp2d(ntmp,otmp,oldtheta,oldr,theta,radius)
            Do ir=1,nr
                data(ip,:,ir) = ntmp(:,nr-ir+1) !Reverse the radius
            EndDo
        EndDo
        !$OMP END DO NOWAIT
        Deallocate(otmp,ntmp)
        !$OMP END PARALLEL

        Print*, 'Writing ', Trim(output_file)
        reclen = nfloat*nr*nth*nphi

        If (nfloat .eq. 1) Then
            Allocate(buff(nphi,nth,nr))
            buff = Real(data)
            Open(interp_unit, file=Trim(output_file), form='unformatted', access='direct', status='unknown',recl=reclen)
            Write(interp_unit,rec=1) buff
            Close(interp_unit)
            Deallocate(buff)
        Else
            Open(interp_unit, file=Trim(output_file), form='unformatted', access='direct', status='unknown',recl=reclen)
            Write(interp_unit,rec=1) data
            Close(interp_unit)
        EndIf
    EndIf

  End Subroutine Interpolate_Old

  Subroutine Interp2d(newarr, arr, ox, oy, nx, ny)
    Implicit None
    Real*8, Dimension(:,:), Intent(InOut) :: newarr !Interpolated data
    Real*8, Dimension(:,:), Intent(In) :: arr !Old data
    Real*8, Dimension(:), Intent(In) :: ox, oy !Old grid
    Real*8, Dimension(:), Intent(In) :: nx, ny !New grid
    Integer :: onx, ony !Size of ox and oy
    Integer :: nnx, nny !Size of nx and ny
    Integer :: ix(4), iy(4)  !Integer storage array for locations of 16 nearest-neighbors on the old grid

    Real*8 :: xc, yc !Current location on the new grid
    Real*8 :: mu, mu2, mu3 !Normalized current grid coordinate
    Real*8 ::  akm1, ak, akp1, akp2 !Storage for old data
    Real*8 :: a0, a1, a2, a3 !Interpolant coefficients
    Real*8 :: ry(4) !Array temporary after x interpolation
    Integer :: i1, i2, ij, ii !loop variables
  
    nnx = Size(nx)
    nny = Size(ny)
    
    onx = Size(ox)
    ony = Size(oy)
    
    Do i2=1,nny
       yc = ny(i2)

       ii = 1
       Do While ((ii .lt. ony) .and. (yc .lt. oy(ii))) 
          ii = ii + 1
       EndDo
       iy(2) = Min(Max(ii-1,1),ony)
       iy(3) = Min(iy(2)+1,ony)
       iy(4) = Min(iy(2)+2,ony)
       iy(1) = Max(iy(2)-1,1)

       Do i1=1,nnx
          xc = nx(i1)        
          !For bicubic convolution interpolation, need nearest 2 points from old grid
          !  in each direction. 
          !Find second index, need linear search as grid is non-uniform starting at previous as guess
          ii = 1
          Do While ((ii .lt. onx) .and. (xc .gt. ox(ii))) 
             ii = ii + 1
          EndDo
          ix(2) = Min(Max(ii-1,1),onx)
          ix(3) = Min(ix(2)+1,onx)
          ix(4) = Min(ix(2)+2,onx)
          ix(1) = Max(ix(2)-1,1)

          mu = (xc-ox(ix(2)))/(ox(ix(3))-ox(ix(2)))
          mu2 = mu*mu
          mu3 = mu*mu2
          Do ij=1,4
             akm1 = arr(ix(1),iy(ij))
             ak   = arr(ix(2),iy(ij))
             akp1 = arr(ix(3),iy(ij))
             akp2 = arr(ix(4),iy(ij))
             
             a3 = ak 
             a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
             a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
             a2 = -0.5d0*akm1 + 0.5d0*akp1
             ry(ij) = a0*mu3+a1*mu2+a2*mu+a3           
          EndDo
          
          mu = (yc-oy(iy(2)))/(oy(iy(3))-oy(iy(2)))
          mu2 = mu*mu
          mu3 = mu*mu2
          
          akm1 = ry(1)
          ak   = ry(2)
          akp1 = ry(3)
          akp2 = ry(4)
          
          a3 = ak 
          a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
          a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
          a2 = -0.5d0*akm1 + 0.5d0*akp1
          newarr(i1,i2) = a0*mu3+a1*mu2+a2*mu+a3     
       EndDo
    EndDo
  End Subroutine Interp2d



  Subroutine Interp3d_Old(ox, oy, oz, nx, ny, nz)
    Real*8, Dimension(:), Intent(In) :: ox, oy, oz !Old grid
    Real*8, Dimension(:), Intent(In) :: nx, ny, nz !New grid
    Integer :: onx, ony, onz !Size of ox and oy
    Integer :: nnx, nny, nnz !Size of nx and ny
    Integer :: ix(4), iy(4), iz(4)  !Integer storage array for locations of 64 nearest-neighbors on the old grid
    Real*8 :: xc, yc, zc, rc, tc, pc !Current location on the new grid
    Real*8 :: mu, mu2, mu3 !Normalized current grid coordinate
    Real*8 :: akm1, ak, akp1, akp2 !Storage for old data
    Real*8 :: a0, a1, a2, a3 !Interpolant coefficients
    Real*8 :: ry(4,4), rz(4) !Array temporary after x interpolation
    Integer :: i1, i2, i3, ii, ij, ik, chunk, nelem !loop variables
    Integer :: ttest

    nnx = Size(nx)
    nny = Size(ny)
    nnz = Size(nz)
    
    onx = Size(ox)
    ony = Size(oy)
    onz = Size(oz)
    chunk = 2
    nelem = 0
    !$OMP PARALLEL DO SHARED(olddata,data,nnx,nny,nnz,onx,ony,onz,nx,ny,nz,ox,oy,oz,chunk,rmax,rmin) &
    !$OMP PRIVATE(i1,i2,i3,ii,ij,ik,zc,yc,xc,pc,tc,rc,ix,iy,iz,mu,mu2,mu3,akm1,ak,akp1,akp2,a0,a1,a2,a3,ry,rz,ttest) &
    !$OMP SCHEDULE(DYNAMIC,chunk)
    Do i3=1,nnz
       zc = nz(i3)
       Do i2=1,nny
          yc = ny(i2)
          Do i1=1,nnx
             xc = nx(i1)

             !Convert cart to sphere
             rc = sqrt(xc*xc+yc*yc+zc*zc)
             If ((rc .gt. rmax).or.(rc .lt. rmin)) Then
                data(i1,i2,i3) = 0d0
                nelem=nelem+1
             Else
                tc = acos(zc/rc)
                pc = atan2(yc,xc)

                !For tricubic convolution interpolation, need nearest 64 points from old grid
                !  in each direction. Find them in O(N) time.
                ii = 1
                Do While ((ii .lt. onz) .and. (rc .lt. oz(ii))) 
                   ii = ii + 1
                EndDo
                iz(2) = Min(Max(ii-1,1),onz)
                iz(3) = Min(iz(2)+1,onz)
                iz(4) = Min(iz(2)+2,onz)
                iz(1) = Max(iz(2)-1,1)

                ii = 1
                Do While ((ii .lt. ony) .and. (tc .lt. oy(ii)))  ! Note -- had to change this due to reversed theta grid (was tc .gt. oy(ii)) Nick F
                   ii = ii + 1
                EndDo
                iy(2) = Min(Max(ii-1,1),ony)
                iy(3) = Min(iy(2)+1,ony)
                iy(4) = Min(iy(2)+2,ony)
                iy(1) = Max(iy(2)-1,1)

                ii = 1
                Do While ((ii .lt. onx) .and. (pc .gt. ox(ii))) 
                   ii = ii + 1
                EndDo
                ix(2) = Min(Max(ii-1,1),onx)
                ix(3) = Min(ix(2)+1,onx)
                ix(4) = Min(ix(2)+2,onx)
                ix(1) = Max(ix(2)-1,1)
                
                mu = (pc-ox(ix(2)))/(ox(ix(3))-ox(ix(2)))
                mu2 = mu*mu
                mu3 = mu*mu2
                Do ik=1,4
                   Do ij=1,4
                      akm1 = olddata(ix(1),iy(ij),iz(ik))
                      ak   = olddata(ix(2),iy(ij),iz(ik))
                      akp1 = olddata(ix(3),iy(ij),iz(ik))
                      akp2 = olddata(ix(4),iy(ij),iz(ik))
                      
                      a3 = ak 
                      a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
                      a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
                      a2 = -0.5d0*akm1 + 0.5d0*akp1
                      ry(ij,ik) = a0*mu3+a1*mu2+a2*mu+a3           
                   EndDo
                EndDo
                
                mu = (tc-oy(iy(2)))/(oy(iy(3))-oy(iy(2)))
                mu2 = mu*mu
                mu3 = mu*mu2             
                Do ik=1,4             
                   akm1 = ry(1,ik)
                   ak   = ry(2,ik)
                   akp1 = ry(3,ik)
                   akp2 = ry(4,ik)
                   
                   a3 = ak 
                   a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
                   a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
                   a2 = -0.5d0*akm1 + 0.5d0*akp1
                   rz(ik) = a0*mu3+a1*mu2+a2*mu+a3
                EndDo
                
                mu = (rc-oz(iz(2)))/(oz(iz(3))-oz(iz(2)))
                mu2 = mu*mu
                mu3 = mu*mu2

                akm1 = rz(1)
                ak   = rz(2)
                akp1 = rz(3)
                akp2 = rz(4)
             
                a3 = ak 
                a0 = -0.5d0*akm1 + 1.5d0*a3 - 1.5d0*akp1 + 0.5d0*akp2
                a1 = akm1 - 2.5d0*a3 + 2d0*akp1 - 0.5d0*akp2
                a2 = -0.5d0*akm1 + 0.5d0*akp1
                data(i1,i2,i3) = a0*mu3+a1*mu2+a2*mu+a3
             EndIf
          EndDo
       EndDo
    EndDo
    !$OMP END PARALLEL DO 

    !Print*, 'Percent zeros:', dble(nelem)/dble(nnx*nny*nnz), '%'
  End Subroutine Interp3d_Old


End Module Interpolation
