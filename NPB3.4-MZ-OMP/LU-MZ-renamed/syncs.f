c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_init( ny, iam, mthreadnum, isync )

c---------------------------------------------------------------------
c   Initialize synchronization variables
c---------------------------------------------------------------------
      use lu_data

      implicit none

      integer ny
      integer isync(0:ny), mthreadnum, iam

c---------------------------------------------------------------------
c---------------------------------------------------------------------

!$    integer, external :: omp_get_thread_num
!$    integer, external :: omp_get_num_threads

      mthreadnum = 0
!$    mthreadnum = omp_get_num_threads() - 1
      if (mthreadnum .gt. ny) mthreadnum = ny
      iam = 0
!$    iam = omp_get_thread_num()
      if (iam.le.mthreadnum) isync(iam) = 0

c      if (global_npb_verbose .ge. 1) then
c            write (*, 100) iam, mthreadnum
c      endif 

100   format("sync_init : Thread ", i5, " total ", i5)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_left( nxmax, ny, nz, v, 
     &                      iam, mthreadnum, isync )

c---------------------------------------------------------------------
c   Thread synchronization for pipeline operation
c---------------------------------------------------------------------

      implicit none

      integer nxmax, ny, nz
      double precision  v(5, nxmax, ny, nz)

      integer isync(0:ny), mthreadnum, iam

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer neigh, iv

c      write (*, 100) "Entering", iam
      if (iam .gt. 0 .and. iam .le. mthreadnum) then
         neigh = iam - 1
!$omp atomic read
         iv = isync(neigh)
         do while (iv .eq. 0)
!$omp atomic read
            iv = isync(neigh)
         end do
!$omp atomic write
         isync(neigh) = 0
      endif
!$omp flush(isync,v)
c      write (*, 100) "Leaving", iam
100   format(A, " sync_left with thread", i5)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine sync_right( nxmax, ny, nz, v,
     &                       iam, mthreadnum, isync )

c---------------------------------------------------------------------
c   Thread synchronization for pipeline operation
c---------------------------------------------------------------------

      implicit none

      integer nxmax, ny, nz
      double precision  v(5, nxmax, ny, nz)

      integer isync(0:ny), mthreadnum, iam

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      integer iv

c      write (*, 100) "Entering", iam
!$omp flush(isync,v)
      if (iam .lt. mthreadnum) then
!$omp atomic read
         iv = isync(iam)
         do while (iv .eq. 1)
!$omp atomic read
            iv = isync(iam)
         end do
!$omp atomic write
         isync(iam) = 1
      endif
c      write (*, 100) "Leaving", iam
100   format(A, " sync_right with thread ", i5)

      return
      end
