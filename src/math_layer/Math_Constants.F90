Module Math_Constants
    Implicit None
    Real*8 :: One_Third = 1.0d0/3.0d0
    Real*8 :: Pi  = 3.1415926535897932384626433832795028841972d0
    Real*8 :: four_pi, over_eight_pi, two_pi
    Real*8 :: Half = 0.5d0
    Real*8 :: Zero = 0.0d0
    Real*8 :: one=1.0d0
Contains
    Subroutine Set_Math_Constants()
        Implicit None
        four_pi = 4.0d0*pi
        two_pi = 2.0d0*pi
        over_eight_pi = 1.0d0/(8.0d0*pi)
    End Subroutine Set_Math_Constants
End Module Math_Constants
