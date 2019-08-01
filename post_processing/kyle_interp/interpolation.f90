Module Interpolation
  !Use IFPORT !If Intel

  Implicit None
  Logical :: cartesian = .False.
  Integer, Parameter :: TrueDouble = Selected_real_kind(14,99)
  Integer, Public :: iteration, iteration_step, initial_iteration, final_iteration 
  Integer, Public :: nr, nth, nphi, onr, nthrd, ncube, nfloat
  Real*8 :: dpi = 3.1415926535897932384626433832795028841972d0
  Real*8 :: rmin, rmax
  Real*8, Allocatable, Dimension(:) :: newx, newy, newz !Old grid
  Real*8, Allocatable, Dimension(:) :: oldtheta, oldr !Old grid
  Real*8, Allocatable, Dimension(:) :: phi, theta, radius !New grid
  Real*8, Allocatable, Dimension(:,:,:) :: olddata, data
  Character(256) :: perm_dir, iter_filename, interp_filename, grid_filename
  Character(24) :: quantity

Contains
  
  Subroutine Finalize_Interp()
    Deallocate(olddata, data, oldtheta, oldr, phi, theta, radius)
    If (cartesian) Then
       Deallocate(newx, newy, newz)
    EndIf
  End Subroutine Finalize_Interp

  Subroutine Make_Iteration_Filename()
    Implicit None
    Character(256) :: fiter
    Character(16) :: temp
    
    Write(temp,fmt='(i8.8)') iteration
    !Write(temp,fmt='(i7.7)') iteration
    fiter = Trim(perm_dir)//Trim(temp)
    grid_filename = Trim(fiter)//'_grid'
    iter_filename = Trim(fiter)//Trim(quantity)
    Write(6,*)'quantity: ', quantity
    If (cartesian) Then
       interp_filename = Trim(fiter)//Trim(quantity)//'_cube'
    Else
       interp_filename = Trim(fiter)//Trim(quantity)//'_interp'
    EndIf
   
  End Subroutine Make_Iteration_Filename

  Subroutine Read_Grid_ASH()
    Implicit None
    Integer :: grid_unit, ntot, reclen, ir, it
    Real(TrueDouble) :: tmptd(3), dth, dr, dphi
    Real(TrueDouble), Allocatable, Dimension(:) :: grid_info
    grid_unit = 11
             
    ! Get the grid
    Print*, 'Reading grid from ', Trim(grid_filename)
    Inquire(iolength=reclen) tmptd
    Open(grid_unit, file=Trim(grid_filename), form='unformatted', access='direct',status='old',recl=reclen)
    !Read grid
    Read(grid_unit,rec=1) tmptd
    onr = Int(tmptd(1))
    nth = Int(tmptd(2))
    nphi = Int(tmptd(3))
    Close(grid_unit)

    Print*, 'Have read old_nr=', onr, ', nth=', nth, ', nphi=', nphi

    ntot = onr + nth + 5
    Allocate(grid_info(ntot), oldr(onr), oldtheta(nth), phi(nphi), theta(nth), radius(nr))

    Inquire(iolength=reclen) grid_info
    Open(grid_unit, file=Trim(grid_filename), form='unformatted', access='direct',status='old',recl=reclen)

    Read(grid_unit,rec=1) grid_info

    ! Now get the actual radial and latitudinal grids
    oldr = grid_info(6:onr+5)
    oldtheta = acos(grid_info(onr+6:))
    Close(grid_unit)
    Deallocate(grid_info)
    
    rmin = Minval(oldr)
    rmax = Maxval(oldr)

    dphi = 2d0*dpi/Dble(nphi)
    dth = dpi/Dble(nth)
    dr = (rmax-rmin)/Dble(nr-1)
    Do it=0,nphi-1
       phi(it+1) = dphi*Dble(it) - dpi
    EndDo
    Do it=0,nth-1
       theta(it+1) = dth*Dble(it)+0.5d0*dth
    EndDo
    Do ir=0,nr-1
       radius(ir+1) = rmax - dr*Dble(ir)
    EndDo

    If (cartesian) Then
       Allocate(olddata(nphi,nth,onr),data(ncube,ncube,ncube),newx(ncube),newy(ncube),newz(ncube))
       dr = 2d0*rmax/Dble(ncube-1)
       Do it=0,ncube-1
          newx(it+1) = dr*Dble(it)-rmax
       EndDo
       newy = newx
       newz = newx
    Else
       Allocate(olddata(nphi,nth,onr),data(nphi,nth,nr))
    EndIf

  End Subroutine Read_Grid_ASH

  Subroutine Read_Grid()
    Implicit None
    Integer :: grid_unit, ntot,  ir, it, etag
    Real(TrueDouble) ::  dth, dr, dphi
    Real(TrueDouble), Allocatable, Dimension(:) ::  oldctheta
    grid_unit = 11
             
    ! Get the grid
    Print*, 'Reading grid from ', Trim(grid_filename)
    Open(grid_unit, file=Trim(grid_filename), form='unformatted', access='stream',status='old')
    Read(grid_unit) etag
    Read(grid_unit) onr
    Read(grid_unit) nth
    Read(grid_unit) nphi



    Print*, 'Have read old_nr=', onr, ', nth=', nth, ', nphi=', nphi, ', etag=',etag


    ntot = onr + nth + 5
    Allocate(oldr(onr), oldtheta(nth), oldctheta(nth), phi(nphi), theta(nth), radius(nr))
    Read(grid_unit) oldr
    Read(grid_unit) oldtheta
    !oldtheta=acos(oldctheta)
    DeAllocate(oldctheta)

    Close(grid_unit)
    
    rmin = Minval(oldr)
    rmax = Maxval(oldr)

    Write(6,*)'rmin/rmax: ', rmin, rmax
    Write(6,*)'oldtheta: ', oldr

    dphi = 2d0*dpi/Dble(nphi)
    dth = dpi/Dble(nth)
    dr = (rmax-rmin)/Dble(nr-1)
    Do it=0,nphi-1
       phi(it+1) = dphi*Dble(it) - dpi
    EndDo
    Do it=0,nth-1
       theta(it+1) = oldtheta(it+1) ! dth*Dble(it)+0.5d0*dth
    EndDo
    Do ir=0,nr-1
       radius(ir+1) = rmax - dr*Dble(ir)
    EndDo

    If (cartesian) Then
       Allocate(olddata(nphi,nth,onr),data(ncube,ncube,ncube),newx(ncube),newy(ncube),newz(ncube))
       dr = 2d0*rmax/Dble(ncube-1)
       Do it=0,ncube-1
          newx(it+1) = dr*Dble(it)-rmax
       EndDo
       newy = newx
       newz = newx
    Else
       Allocate(olddata(nphi,nth,onr),data(nphi,nth,nr))
    EndIf

  End Subroutine Read_Grid




  Subroutine Interpolate()
    Implicit None
    Integer :: interp_unit, ir, ip, chunk
    Integer*4 :: reclen
    Real*8, Allocatable, Dimension(:,:) :: otmp, ntmp
    Real*4, Allocatable, Dimension(:,:,:) :: buff
    interp_unit=12
    olddata = 0d0
    data = 0d0

    Print*, 'Reading ', Trim(iter_filename)
    reclen = nfloat*onr*nth*nphi
    If (nfloat .eq. 1) Then
       Allocate(buff(nphi,nth,onr))
       !Read in bulk data


       Open(interp_unit, file=Trim(iter_filename), form='unformatted', access='stream', status='unknown')
       Read(interp_unit) buff

       Close(interp_unit)
       olddata = Dble(buff)
       Deallocate(buff)
    Else
       !Read in bulk data


       Write(6,*)'Opening: ', iter_filename
       Open(interp_unit, file=Trim(iter_filename), form='unformatted', access='stream', status='unknown')
       Read(interp_unit) olddata

       Close(interp_unit)
    EndIf

    If (cartesian) Then
        Print*, 'Min(data)=', Minval(data), 'Max(data)=', Maxval(data)
        Print*, 'Min(olddata)=', Minval(olddata), 'Max(olddata)=', Maxval(olddata)

       Call Interp3d(phi,oldtheta,oldr,newx,newy,newz)

       Print*, 'Min(data)=', Minval(data), 'Max(data)=', Maxval(data)
       Print*, 'Writing ', Trim(interp_filename)
       reclen = nfloat*ncube**3

       If (nfloat .eq. 1) Then
          Allocate(buff(ncube,ncube,ncube))
          buff = Real(data)
          Inquire(iolength=reclen) buff
          Open(interp_unit, file=Trim(interp_filename), form='unformatted', access='direct', status='unknown',recl=reclen)
          Write(interp_unit,rec=1) buff
          Close(interp_unit)
          Deallocate(buff)
       Else
          Inquire(iolength=reclen) data
          Open(interp_unit, file=Trim(interp_filename), form='unformatted', access='direct', status='unknown',recl=reclen)
          Write(interp_unit,rec=1) data
          Close(interp_unit)
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

       Print*, 'Writing ', Trim(interp_filename)
       reclen = nfloat*nr*nth*nphi

       If (nfloat .eq. 1) Then
          Allocate(buff(nphi,nth,nr))
          buff = Real(data)
          Open(interp_unit, file=Trim(interp_filename), form='unformatted', access='direct', status='unknown',recl=reclen)
          Write(interp_unit,rec=1) buff
          Close(interp_unit)
          Deallocate(buff)
       Else
          Open(interp_unit, file=Trim(interp_filename), form='unformatted', access='direct', status='unknown',recl=reclen)
          Write(interp_unit,rec=1) data
          Close(interp_unit)
       EndIf
    EndIf

  End Subroutine Interpolate

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

  Subroutine Interp3d(ox, oy, oz, nx, ny, nz)
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
  
    nnx = Size(nx)
    nny = Size(ny)
    nnz = Size(nz)
    
    onx = Size(ox)
    ony = Size(oy)
    onz = Size(oz)
    chunk = 2
    nelem = 0
    !$OMP PARALLEL DO SHARED(olddata,data,nnx,nny,nnz,onx,ony,onz,nx,ny,nz,ox,oy,oz,chunk,rmax,rmin) &
    !$OMP PRIVATE(i1,i2,i3,ii,ij,ik,zc,yc,xc,pc,tc,rc,ix,iy,iz,mu,mu2,mu3,akm1,ak,akp1,akp2,a0,a1,a2,a3,ry,rz) &
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

    Print*, 'Percent zeros:', dble(nelem)/dble(nnx*nny*nnz), '%'
  End Subroutine Interp3d

End Module Interpolation
