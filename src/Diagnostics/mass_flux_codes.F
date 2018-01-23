    Integer, Parameter :: rhovoffset = voffset+200 ! :OFFSET CODE: 
    !------------ Mass Flux ---------------------!
    ! Full
    Integer, Parameter :: rhov_r      = rhovoffset+1 ! :tex: $\mathrm{f}_1v_r$
    Integer, Parameter :: rhov_theta  = rhovoffset+2 ! :tex: $\mathrm{f}_1v_\theta$
    Integer, Parameter :: rhov_phi    = rhovoffset+3 ! :tex: $\mathrm{f}_1v_\phi$
    ! Fluctuating
    Integer, Parameter :: rhovp_r     = rhovoffset+4 ! :tex: $\mathrm{f}_1v_r'$ 
    Integer, Parameter :: rhovp_theta = rhovoffset+5 ! :tex: $\mathrm{f}_1v_\theta'$
    Integer, Parameter :: rhovp_phi   = rhovoffset+6 ! :tex: $\mathrm{f}_1v_\phi'$
    ! Mean
    Integer, Parameter :: rhovm_r     = rhovoffset+7 ! :tex: $\mathrm{f}_1\overline{v_r}$
    Integer, Parameter :: rhovm_theta = rhovoffset+8 ! :tex: $\mathrm{f}_1\overline{v_\theta}$
    Integer, Parameter :: rhovm_phi   = rhovoffset+9 ! :tex: $\mathrm{f}_1\overline{v_\phi}$
