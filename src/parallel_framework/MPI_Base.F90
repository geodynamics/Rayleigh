Module MPI_Base
	Type communicator
		Integer :: comm	! The mpi handle for this group
		Integer :: np	! The number of processors in this group
		Integer :: rank ! A processes's local rank within this group
	End Type communicator
  
	Include 'mpif.h'
	Private ::  mpi_null_delete_fn, mpi_dup_fn, mpi_null_copy_fn
	Public :: mpi_wtime,mpi_wtick
End Module MPI_Base
