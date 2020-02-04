c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine ssor(u, rsd, frct, qs, rho_i, tv, a, b, c, d,
     $                au, bu, cu, du, nx, nxmax, ny, nz, isync)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   to perform pseudo-time stepping SSOR iterations
c   for five nonlinear pde's.
c---------------------------------------------------------------------

      use lu_data
      implicit none

      integer          nx, nxmax, ny, nz, isync(0:global_problem_size)
      double precision u(5,nxmax,ny,nz), rsd(5,nxmax,ny,nz), 
     $                 frct(5,nxmax,ny,nz), qs(nxmax,ny,nz), 
     $                 rho_i(nxmax,ny,nz), tv(5,2:nxmax-1,ny),
     $                 a (5,5,2:nxmax-1,ny), b (5,5,2:nxmax-1,ny), 
     $                 c (5,5,2:nxmax-1,ny), d (5,5,2:nxmax-1,ny),
     $                 au(5,5,2:nxmax-1,ny), bu(5,5,2:nxmax-1,ny), 
     $                 cu(5,5,2:nxmax-1,ny), du(5,5,2:nxmax-1,ny)

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer mthreadnum, iam
      integer i, j, k, m
      double precision  tmp
      external timer_read
      double precision timer_read
      integer thread_num, othread

!$    integer, external :: omp_get_thread_num

 
!$    othread = omp_get_thread_num()
c      write (*, 101) thread_num
101   format("ssor (begin): I am thread ", i5)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,i,j,k,tmp,mthreadnum,iam,
!$OMP& thread_num)
!$OMP MASTER
      if (global_timeron) call timer_start(global_t_rhs)
!$OMP END MASTER
      do k = 2, nz-1
!$OMP DO SCHEDULE(STATIC)
         do j = 2, ny-1
            do i = 2, nx-1
               do m = 1, 5
                  rsd(m,i,j,k) = global_dt * rsd(m,i,j,k)
               end do
            end do
         end do
!$OMP END DO nowait
      end do
!$OMP MASTER
      if (global_timeron) call timer_stop(global_t_rhs)
!$OMP END MASTER


      call sync_init( global_problem_size, iam, mthreadnum, isync )
!$OMP BARRIER

!$    thread_num = omp_get_thread_num()
c      write (*, 100) thread_num
100   format("ssor (after sync_init): I am thread ", i5)
!$omp barrier
      do k = 2, nz-1
c---------------------------------------------------------------------
c   form the lower triangular part of the jacobian matrix
c---------------------------------------------------------------------
!$OMP MASTER
         if (global_timeron) call timer_start(global_t_jacld)
!$OMP END MASTER
         write (*, 103) thread_num, othread, k, nx,ny,nz
103      format ("[beforejacld] Thread ", i1, " (outer thread ", i1,
     >  ")] Looping", 
     >           " with k = ", i5, " and n is ", i5,i5,i5)

         call jacld(k, u, rho_i, qs, a, b, c, d, nx, nxmax, ny, nz)

!$OMP MASTER
         if (global_timeron) call timer_stop(global_t_jacld)
 
c---------------------------------------------------------------------
c   perform the lower triangular solution
c---------------------------------------------------------------------
         if (global_timeron) call timer_start(global_t_blts)
!$OMP END MASTER

         call sync_left( nxmax, ny, nz, rsd, iam, mthreadnum, isync )

         call blts( nx, nxmax, ny, nz, k, global_omega, rsd, a, b, c,
     >              d)
         write (*, 104) thread_num, othread, k, nx,ny,nz
104      format ("[afterblts] Thread ", i1, " (outer thread ", i1,
     >  ")] Looping", 
     >           " with k = ", i5, " and n is ", i5,i5,i5)
     
         call sync_right( nxmax, ny, nz, rsd, iam, mthreadnum, isync )

!$OMP MASTER
         if (global_timeron) call timer_stop(global_t_blts)
!$OMP END MASTER
      end do
 
!$OMP BARRIER


      do k = nz-1, 2, -1
c---------------------------------------------------------------------
c   form the strictly upper triangular part of the jacobian matrix
c---------------------------------------------------------------------
!$OMP MASTER
         if (global_timeron) call timer_start(global_t_jacu)
!$OMP END MASTER

         call jacu(k, u, rho_i, qs, au, bu, cu, du, 
     $             nx, nxmax, ny, nz)

!$OMP MASTER
         if (global_timeron) call timer_stop(global_t_jacu)

c---------------------------------------------------------------------
c   perform the upper triangular solution
c---------------------------------------------------------------------
         if (global_timeron) call timer_start(global_t_buts)
!$OMP END MASTER

         call sync_left( nxmax, ny, nz, rsd, iam, mthreadnum, isync )

         call buts( nx, nxmax, ny, nz, k, global_omega, rsd, tv, 
     $              du, au, bu, cu)

         call sync_right( nxmax, ny, nz, rsd, iam, mthreadnum, isync )

!$OMP MASTER
         if (global_timeron) call timer_stop(global_t_buts)
!$OMP END MASTER
      end do

!$OMP BARRIER


c---------------------------------------------------------------------
c   update the variables
c---------------------------------------------------------------------

      tmp = 1.0d0 / ( global_omega * ( 2.0d0 - global_omega ) ) 

!$OMP MASTER
      if (global_timeron) call timer_start(global_t_add)
!$OMP END MASTER
      do k = 2, nz-1
!$OMP DO SCHEDULE(STATIC)
         do j = 2, ny-1
            do i = 2, nx-1
               do m = 1, 5
                  u(m,i,j,k) = u(m,i,j,k) + tmp * rsd(m,i,j,k)
               end do
            end do
         end do
!$OMP END DO nowait
      end do
!$OMP MASTER
      if (global_timeron) call timer_stop(global_t_add)
!$OMP END MASTER
!$OMP END PARALLEL

!$    thread_num = omp_get_thread_num()
c      write (*, 102) thread_num
102   format("ssor (end): I am thread ", i5)

c---------------------------------------------------------------------
c   compute the steady-state residuals
c---------------------------------------------------------------------
      if (global_timeron) call timer_start(global_t_rhs)
      call rhs(u, rsd, frct, qs, rho_i, nx, nxmax, ny, nz)
      if (global_timeron) call timer_stop(global_t_rhs)


      return
      end
