Module Structures
	!//////////////////////////////////////////////////////////
	! This module contains various data structures that do not
	!   have associated methods.
    Type, Public :: rmcontainer1d
        Real*8, Allocatable :: data(:)
    End Type rmcontainer1d

	Type, Public :: rmcontainer
		Real*8, Allocatable :: data(:,:)
	End Type rmcontainer

	Type, Public :: rmcontainer3d
		Real*8, Allocatable :: data(:,:,:)
	End Type rmcontainer3d

	Type, Public :: rmcontainer4d
		Real*8, Allocatable :: data(:,:,:,:)
	End Type rmcontainer4d

End Module Structures
