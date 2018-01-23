#define DO_PSI Do t = my_theta%min, my_theta%max; Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max; Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
#define DDBUFF d2buffer%p3a
Module Diagnostics_Velocity_Diffusion
    Use Diagnostics_Base
    Implicit None


Contains
    !//////////////////////////////////////////////////////////////////////////////
    ! Note:  These diagnostic quantities were programmed shortly before 
    !        Rayleigh's release.  They were programmed with readability, and
    !        not efficiency or vectorization, in mind.

    Subroutine Compute_Velocity_Diffusion(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: mypi, amp,del2u, estress
        Real*8, Allocatable :: mu_visc(:), dmudr(:), ovstheta(:), ovs2theta(:)

        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2
        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta

        !Compute the dynamic viscosity mu=rho*nu (nu is kinematic viscosity)
        Allocate(mu_visc(1:N_R), dmudr(1:N_R))

        mu_visc = ref%density*nu
        dmudr = mu_visc*ref%density*(ref%dlnrho+dlnu)


        !////////////////////////////////////////////////////////
        ! Full viscous forces
        ! r-direction
        If (compute_quantity(viscous_force_r)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = DDBUFF(PSI,dvrdrdr)+Two_Over_R(r)*buffer(PSI,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvrdtdt)+cottheta(t)*buffer(PSI,dvrdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        buffer(PSI,vr) + & 
                        buffer(PSI,dvtdt)+buffer(PSI,vtheta)*cottheta(t) + &
                        ovstheta(t)*buffer(PSI,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(buffer(PSI,dvrdr)*ref%dlnrho(r)+ &
                        buffer(PSI,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = buffer(PSI,dvrdr)-One_Third*buffer(PSI,vr)*ref%dlnrho(r)

                qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


            END_DO

            Call Add_Quantity(qty)
        Endif

        !Theta-direction
        If (compute_quantity(viscous_force_theta)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_theta}
                del2u = DDBUFF(PSI,dvtdrdr)+Two_Over_R(r)*buffer(PSI,dvtdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvtdtdt)+cottheta(t)*buffer(PSI,dvtdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u +OneOverRSquared(r)*( -2.0d0*buffer(PSI,dvrdt) + &
                        ovs2theta(t)*(   buffer(PSI,vtheta) + &
                        2.0d0*costheta(t)*buffer(PSI,dvpdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*buffer(PSI,dvrdt)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(buffer(PSI,dvrdt)-buffer(PSI,vtheta) ) +buffer(PSI,dvtdr)

                qty(PSI) = dmudr(r)*estress +mu_visc(r)*del2u
            END_DO

            Call Add_Quantity(qty)
        Endif


        !Phi-direction
        If (compute_quantity(viscous_force_phi)) Then

            DO_PSI
                del2u = DDBUFF(PSI,dvpdrdr)+Two_Over_R(r)*buffer(PSI,dvpdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvpdtdt)+cottheta(t)*buffer(PSI,dvpdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvpdpdp)*ovs2theta(t)


                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*ovstheta(t)*buffer(PSI,dvrdp)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(ovstheta(t)*buffer(PSI,dvrdp)-buffer(PSI,vphi) )+ &
                          buffer(PSI,dvpdr)

                qty(PSI) =dmudr(r)*estress + mu_visc(r)*del2u
            END_DO

            Call Add_Quantity(qty)
        Endif


        DeAllocate(mu_visc, dmudr)
        DeAllocate(ovstheta,ovs2theta)
    End Subroutine Compute_Velocity_Diffusion

End Module Diagnostics_Velocity_Diffusion
