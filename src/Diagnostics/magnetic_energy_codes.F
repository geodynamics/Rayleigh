    !///////////////////////////////////////////////////////////
    !           Magnetic Energies
    Integer, Parameter :: meoffset = joffset+100 ! :OFFSET CODE:

    Integer, Parameter :: magnetic_energy = meoffset+1 ! :tex: $\frac{1}{2}c_4\boldsymbol{B}^2$
    Integer, Parameter :: radial_me       = meoffset+2 ! :tex: $\frac{1}{2}c_4B_r^2$
    Integer, Parameter :: theta_me        = meoffset+3 ! :tex: $\frac{1}{2}c_4B_\theta^2$
    Integer, Parameter :: phi_me        = meoffset+4 ! :tex: $\frac{1}{2}c_4B_\phi^2$

    Integer, Parameter :: mmagnetic_energy = meoffset+5 ! :tex: $\frac{1}{2}c_4\overline{\boldsymbol{B}}^2$
    Integer, Parameter :: radial_mme       = meoffset+6 ! :tex: $\frac{1}{2}c_4\overline{B_r}^2$
    Integer, Parameter :: theta_mme        = meoffset+7 ! :tex: $\frac{1}{2}c_4\overline{B_\theta}^2$
    Integer, Parameter :: phi_mme        = meoffset+8 ! :tex: $\frac{1}{2}c_4\overline{B_\phi}^2$

    Integer, Parameter :: pmagnetic_energy = meoffset+9  ! :tex: $\frac{1}{2}c_4\boldsymbol{B'}^2$
    Integer, Parameter :: radial_pme       = meoffset+10 ! :tex: $\frac{1}{2}c_4{B_r'}^2$
    Integer, Parameter :: theta_pme        = meoffset+11 ! :tex: $\frac{1}{2}c_4{B_\theta'}^2$
    Integer, Parameter :: phi_pme        = meoffset+12 ! :tex: $\frac{1}{2}c_4{B_\phi'}^2$

