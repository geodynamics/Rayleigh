PROGRAM HELLO
 
 INTEGER NTHREADS, TID, OMP_GET_NUM_THREADS,
 + OMP_GET_THREAD_NUM
 
C Fork a team of threads giving them their own copies of variables
!$OMP PARALLEL PRIVATE(NTHREADS, TID)


C Obtain thread number
 TID = OMP_GET_THREAD_NUM()
 PRINT *, 'Hello World from thread = ', TID

C Only master thread does this
 IF (TID .EQ. 0) THEN
  NTHREADS = OMP_GET_NUM_THREADS()
 PRINT *, 'Number of threads = ', NTHREADS
 END IF

C All threads join master thread and disband
!$OMP END PARALLEL

 END

