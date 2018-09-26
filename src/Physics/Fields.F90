!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

Module Fields
    Use Parallel_Framework
    Use ProblemSize
    Use Controls
    Use Spherical_Buffer
    !///////////////////////////////////////////////////////////
    !  The scalable framework of this program is centered around
    !  a large buffer that holds multiple fields.
    !  The buffer is reshaped during each transpose, and fields are
    !  added and removed as the buffer moves amongst the different
    !  configuration spaces.

    !  The field type is meta data that allows us to know where
    !  a field is in the buffer at a particular time, or if it's
    !  even there at all
    Implicit None
    Type Field_Indexer
        Integer :: c1a_counter = 0
        Integer :: c1b_counter = 0
        Integer :: c2a_counter = 0
        Integer :: c2b_counter = 0
        Integer :: c3a_counter = 0
        Integer :: c3b_counter = 0


        Integer :: global_counter = 0

        Contains
            Procedure :: Add_Field

    End Type

    Type(Field_Indexer) :: wsp_indices, co_indices

    !////////////////////////////////////////////////////////////////////////////
    ! Variable locations in the global field meta data buffer
    Integer :: wvar, pvar, tvar, zvar
    Integer :: dpdr1  ! We reserve dpdr for the final location at output time
    Integer :: dwdr,  d3wdr3, dtdr, dzdr, d2zdr2, d2tdr2, d2wdr2
    Integer :: vr, vtheta, vphi,dtdt,dvrdt,dvtdr,dvpdr,dvrdr
    Integer :: dvrdp, dvtdp, dvpdp, dtdp

    Integer :: dvtdt, dvpdt

    Integer :: avar, dadr, d2adr2 ! Toroidal magnetic streamfunction and its radial derivatives
    Integer :: cvar, dcdr, d2cdr2 !Poloidal magnetic streamfunction and its radial derivatives

    Integer :: br,btheta,bphi, curlbr,curlbtheta,curlbphi
    Integer :: emfr,emftheta,emfphi
    !///////////////////////////////////////////////////////////////////////////



    !==============================================================================
    !        Integer parameters used to reference specific equations and variables
    !                FOR LINEAR Solve
    !==============================================================================
    Integer, parameter :: weq = 1,  peq = 2,  teq = 3
    Integer, parameter :: zeq = 4,  ceq = 5,  aeq = 6




    Type(SphericalBuffer) :: wsp ! Primary workspace for the entire run
    Integer :: wsfcount(3,2)     ! wsp field count for each configuration

    Type(SphericalBuffer) :: cobuffer ! Workspace for holding output-related variables
                                      ! not already contained in wsp
    Integer :: cbfcount(3,2)    ! cobuffer field count for each configuration
    !======================================================================
    !  Finally, we have some diagnostic field variables
    !  These indices are used to reference the values held in diag_fields s2a and p3a
    !  Integer :: dbrdr_ia, dbtdr_ia, dbpdr_ia, dbrdt_ia,avar_ia, dpdr_ia
    !  Or maybe use "db" for diagnostics_buffer
    Integer :: dpdr_cb, dpdt_cb
    Integer :: dbrdr_cb, dbtdr_cb, dbpdr_cb, dbrdt_cb,avar_cb
    !  These indices are used to reference the expanded workspace buffer during
    !  output iterations (for configurations s3a and p3a only)
    !  They are set in Diagnostics_Base.F90
    Integer :: dpdr, dpdt, dpdp  ! pressure derivatives
    Integer :: dbrdr, dbrdt, dbrdp, dbtdr, dbtdt, dbtdp, dbpdr, dbpdt, dbpdp
    Integer :: output_nextra  ! wsp%p3a hold output_nextra more fields at output time
    Integer :: vindex(1:12), bindex(1:12)  ! the indices within the p3a buffer
    !  (at output time) corresponding to v,B, and their derivatives

Contains

    Subroutine Add_Field(self,ind,cfg)
        Implicit None
        ! This routine keeps a tally of field counts
        ! for each configuration and assigns a value
        ! to ind (a field reference integer)
        Class(Field_Indexer)    :: self
        Integer, intent(InOut)  :: ind
        Character*3, Intent(In) :: cfg

        self%global_counter = self%global_counter+1
        ind = self%global_counter
        Select Case(cfg)
            !Forward Leg Configurations
            Case('p1a')
                self%c1a_counter = self%c1a_counter+1
                self%c2a_counter = self%c2a_counter+1
                self%c3a_counter = self%c3a_counter+1
            Case('p2a')
                self%c2a_counter = self%c2a_counter+1
                self%c3a_counter = self%c3a_counter+1
            Case('p3a')
                self%c3a_counter = self%c3a_counter+1

            !Return Leg Configurations
            Case('p3b')
                self%c3b_counter = self%c3b_counter+1
                self%c2b_counter = self%c2b_counter+1
                self%c1b_counter = self%c1b_counter+1
            Case('p2b')
                self%c2b_counter = self%c2b_counter+1
                self%c1b_counter = self%c1b_counter+1
            Case('p1b')
                self%c1b_counter = self%c1b_counter+1
        End Select

    End Subroutine Add_Field



    Subroutine Initialize_Field_Structure()
        Implicit None
        Character*3 :: config


        !This routine serves two purposes:
        ! 1 - It assigns values to each of our field reference integers
        ! 2 - It tallies the number of fields added to the buffer in
        !     in each configuration, so that the necessary buffer size
        !     may be computed.

        !  Note:  wsp_indices is a simple object used for keeping track of
        !  of the field indexing.  This indexing be done manually
        !  (albeit carefully) without wsp_indices.  The advantage here
        !  is that additional fields can be added at any stage without
        !  manually shifting around all the indices.

        !  Once a field is added using Add_Field, its memory space
        !  is assumed to persist throughout that leg of the loop.


        !//////////////////////////////////
        !  First, we need an accounting of all fields that will be stored
        !  (even temporarily) in the p1a buffer.  We will assume that
        !  These fields persist out to p3a for now


        config = 'p1a'
        !//////////////////////////////
        !  It is rather important that the primary fields
        !  we solve for are added first
        !  Their numbers should be 1 through nfields and it is OK
        !  if they are not used in all configurations.  We want the
        !  equation numbering and the field numbering to agree - always
        !  Add the primary fields - this should ALWAYS come first
        Call wsp_indices%Add_Field(Wvar , config)
        Call wsp_indices%Add_Field(Pvar , config)
        Call wsp_indices%Add_Field(Tvar , config)
        Call wsp_indices%Add_Field(Zvar , config)
        If (magnetism) Then
          Call wsp_indices%Add_Field(cvar , config)
          Call wsp_indices%Add_field(avar , config)
        Endif

        Call wsp_indices%Add_Field(d3Wdr3 , config)
        Call wsp_indices%Add_Field(dPdr1  , config)
        Call wsp_indices%Add_field(d2Zdr2 , config)
        Call wsp_indices%Add_Field(d2Wdr2 , config)


        If (magnetism) Then
          Call wsp_indices%Add_field(dcdr   , config)
          Call wsp_indices%Add_field(dadr   , config)
          Call wsp_indices%Add_field(d2adr2 , config)
        Endif


        !//////////////////////////////////////////////////////////
        !  Next, we want to account for fields that we build in s2a/p2a (many are d by dtheta fields)
        config = 'p2a'
        Call wsp_indices%Add_Field(vtheta , config)
        Call wsp_indices%Add_Field(vphi   , config)
        Call wsp_indices%Add_Field(dvtdr  , config)
        Call wsp_indices%Add_Field(dvpdr  , config)
        If (magnetism) Then

            Call wsp_indices%Add_Field(curlbr     ,config)

        Endif
        !///////////////////////////////////////////////////////
        !  Next we have fields of the d_by_dphi variety
        !  that we add to the p3a buffer
        config = 'p3a'
        Call wsp_indices%Add_Field(dvrdp,config)
        Call wsp_indices%Add_field(dvtdp,config)
        Call wsp_indices%Add_field(dvpdp,config)
        Call wsp_indices%Add_Field(dvtdt,config)
        Call wsp_indices%Add_field(dvpdt,config)
        Call wsp_indices%Add_field(dtdp,config)

        !//////////////////////////////////////////////////////////
        !   Throughout the forward loop, many variables are replaced
        !   with new variables (e.g., vr overwrites W to save memory).
        !   Those overwrites are handled here, as they do not modify
        !   the buffer size.

        d2Tdr2 = dPdr1
        dWdr   = d3Wdr3
        dTdr   = dPdr1 !replaces d2tdr2, which replaced dpdr1
        dZdr   = d2Zdr2

        vr    = wvar
        dvrdr = dwdr
        dtdt  = d2wdr2
        dvrdt = dzdr

        If (magnetism) Then
            d2cdr2 = d2adr2
            Br     = cvar
            Btheta = avar ! might need to rethink this one
            Bphi   = dcdr
            curlbtheta = d2cdr2
            curlbphi   = dadr
        Endif

        !///////////////////////////////////////////////////////////////
        ! Now that all fields have been added for the forward leg of the loop,
        ! we have a count of how large the buffer needs to be in each
        ! configuration

        wsfcount(1,1) = wsp_indices%c1a_counter
        wsfcount(2,1) = wsp_indices%c2a_counter
        wsfcount(3,1) = wsp_indices%c3a_counter

        !//////////////////////////////////////////////
        ! The Return Trip is much simpler.
        ! We have computed the explicit portion of RHS
        ! for each equation in physical space, and
        ! we need to send it back to spectral space.

        ! These following code should pretty much never be modified by the user.
        if (.not. magnetism) then
            wsfcount(1,2) = 4        ! four RHS's go back for the solve
            wsfcount(2,2) = 4
            wsfcount(3,2) = 4
        else
            emfr     = avar
            emftheta = cvar
            emfphi   = avar+1
            wsfcount(1,2) = 7        ! seven RHS's go back for the solve (1 field is differentiated and combined at the end)
            wsfcount(2,2) = 7
            wsfcount(3,2) = 7
        endif
        !Write(6,*)'Fields initialized'
        !Write(6,*)'c1a_counter is: ', wsp_indices%c1a_counter
        !Write(6,*)'c2a_counter is: ', wsp_indices%c2a_counter
        !Write(6,*)'c3a_counter is: ', wsp_indices%c3a_counter

    End Subroutine Initialize_Field_Structure
    Subroutine Initialize_Diagnostic_Indices()
        Implicit None
        Character*3 :: config
        vindex(1)  = vr
        vindex(2)  = vtheta
        vindex(3)  = vphi
        vindex(4)  = dvrdr
        vindex(5)  = dvrdt
        vindex(6)  = dvrdp
        vindex(7)  = dvtdr
        vindex(8)  = dvtdt
        vindex(9)  = dvtdp
        vindex(10) = dvpdr
        vindex(11) = dvpdt
        vindex(12) = dvpdp


        !here we have some indices within the cobuffer
        ! Config p1a
        config = 'p1a'
        Call co_indices%Add_Field(dpdr_cb,config)


        ! Config p2a
        config = 'p2a'
        Call co_indices%Add_Field(dpdt_cb,config)
        If (magnetism) Then
            Call co_indices%Add_Field(dbrdr_cb , config)
            Call co_indices%Add_Field(dbtdr_cb , config)
            Call co_indices%Add_Field(dbpdr_cb , config)
            Call co_indices%Add_Field(dbrdt_cb , config)
            Call co_indices%Add_Field( avar_cb , config)
        Endif
        cbfcount(:,:) = 0
        cbfcount(1,1) = co_indices%c1a_counter
        cbfcount(2,1) = co_indices%c2a_counter
        cbfcount(3,1) = co_indices%c3a_counter

        !  Finally, we adjust the field count for the p3a buffer
        !  at output time.  Note that the basic field count for
        !  non-output iterations was already set in
        !  Initialize_Field_Structure.
        !
        !  We now continue to increment the counters.
        config = 'p3a'
        Call wsp_indices%Add_Field(dpdr , config)
        Call wsp_indices%Add_Field(dpdt , config)
        Call wsp_indices%Add_Field(dpdp , config)


        If (magnetism) Then

            Call wsp_indices%Add_Field(dbrdr , config)
            Call wsp_indices%Add_Field(dbtdr , config)
            Call wsp_indices%Add_field(dbpdr , config)

            Call wsp_indices%Add_Field(dbrdt , config)
            Call wsp_indices%Add_Field(dbtdt , config)
            Call wsp_indices%Add_field(dbpdt , config)

            Call wsp_indices%Add_Field(dbrdp , config)
            Call wsp_indices%Add_Field(dbtdp , config)
            Call wsp_indices%Add_field(dbpdp , config)

            bindex(1)  = br
            bindex(2)  = btheta
            bindex(3)  = bphi
            bindex(4)  = dbrdr
            bindex(5)  = dbrdt
            bindex(6)  = dbrdp
            bindex(7)  = dbtdr
            bindex(8)  = dbtdt
            bindex(9)  = dbtdp
            bindex(10) = dbpdr
            bindex(11) = dbpdt
            bindex(12) = dbpdp
        Endif
        !Write(6,*)'NEW CHECK: ', wsp%nf3a, wsp_indices%c3a_counter
        output_nextra = wsp_indices%c3a_counter- wsfcount(3,1)
        !Write(6,*)'CHECK:  ', output_nextra, wsp%nf3a, wsfcount(3,1)
    End Subroutine Initialize_Diagnostic_Indices

End Module Fields
