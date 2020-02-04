!-------------------------------------------------------------------------!
!                                                                         !
!        N  A  S     P A R A L L E L     B E N C H M A R K S  3.4         !
!                                                                         !
!          O p e n M P    M U L T I - Z O N E    V E R S I O N            !
!                                                                         !
!                           L U - M Z - O M P                             !
!                                                                         !
!-------------------------------------------------------------------------!
!                                                                         !
!    This benchmark is an OpenMP version of the NPB LU code.              !
!    Refer to NAS Technical Reports 95-020 for details.                   !
!                                                                         !
!    Permission to use, copy, distribute and modify this software         !
!    for any purpose with or without fee is hereby granted.  We           !
!    request, however, that all derived work reference the NAS            !
!    Parallel Benchmarks 3.4. This software is provided "as is"           !
!    without express or implied warranty.                                 !
!                                                                         !
!    Information on NPB 3.4, including the technical report, the          !
!    original specifications, source code, results and information        !
!    on how to submit new results, is available at:                       !
!                                                                         !
!           http://www.nas.nasa.gov/Software/NPB/                         !
!                                                                         !
!    Send comments or suggestions to  npb@nas.nasa.gov                    !
!                                                                         !
!          NAS Parallel Benchmarks Group                                  !
!          NASA Ames Research Center                                      !
!          Mail Stop: T27A-1                                              !
!          Moffett Field, CA   94035-1000                                 !
!                                                                         !
!          E-mail:  npb@nas.nasa.gov                                      !
!          Fax:     (650) 604-3957                                        !
!                                                                         !
!-------------------------------------------------------------------------!

c---------------------------------------------------------------------
c
c Authors: S. Weeratunga
c          V. Venkatakrishnan
c          E. Barszcz
c          M. Yarrow
C          R.F. Van der Wijngaart
C          H. Jin
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
      program LU_MZ
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   driver for the performance evaluation of the solver for
c   five coupled parabolic/elliptic partial differential equations.
c
c---------------------------------------------------------------------

      use lu_data
      use lu_fields
      use ompnpb

      implicit none

      integer num_zones

      logical verified
      double precision mflops, nsur, navg, n3

      integer i, zone, step, iz, tot_threads, itimer, nthreads
      double precision t, tmax, timer_read, trecs(global_t_last),
     $                 tsum(global_t_last), tming(global_t_last),
     >                 tmaxg(global_t_last),
     $                 rsdnm(5), rsdnm_aux(5), errnm(5), errnm_aux(5),
     $                 frc, frc_aux
      external timer_read
      character t_names(global_t_last)*8
      integer thread_num

!$    integer, external :: omp_get_thread_num


c---------------------------------------------------------------------
c   read input data
c---------------------------------------------------------------------
      call read_input(tot_threads, itimer)

      if (global_timeron) then
         t_names(global_t_total)  = 'total'
         t_names(global_t_rhsx)   = 'rhsx'
         t_names(global_t_rhsy)   = 'rhsy'
         t_names(global_t_rhsz)   = 'rhsz'
         t_names(global_t_rhs)    = 'rhs'
         t_names(global_t_jacld)  = 'jacld'
         t_names(global_t_blts)   = 'blts'
         t_names(global_t_jacu)   = 'jacu'
         t_names(global_t_buts)   = 'buts'
         t_names(global_t_add)    = 'add'
         t_names(global_t_l2norm) = 'l2norm'
         t_names(global_t_rdis1)  = 'qbc_copy'
         t_names(global_t_rdis2)  = 'qbc_comm'
      endif

c---------------------------------------------------------------------
c   set up domain sizes
c---------------------------------------------------------------------
      call zone_setup(global_nx, global_nxmax, global_ny, global_nz)

      num_zones = global_max_zones
      call setup_omp(num_zones, global_nx, global_ny, global_nz,
     >               tot_threads)
      write (*, *) global_myid
      call zone_starts(num_zones, global_nx, global_nxmax, global_ny,
     >                 global_nz)

c---------------------------------------------------------------------
c      allocate space for field arrays
c---------------------------------------------------------------------
       call alloc_field_space

c---------------------------------------------------------------------
c   set up coefficients
c---------------------------------------------------------------------
      call setcoeff()


      if (global_timeron) then
         do i = 1, global_t_last
            tsum(i) = 0.d0
            tming(i) = huge(0.d0)
            tmaxg(i) = 0.d0
         end do
      endif

c---------------------------------------------------------------------
c   start of the outer parallel region
c---------------------------------------------------------------------
!$omp parallel
!$omp& private(iz,i,zone,step,t,tmax,trecs,nthreads,global_isync,
!$omp&  global_a,global_b,global_c,global_d,global_au,global_bu,
!$omp&  global_cu,global_du,global_tv,global_phi1,global_phi2,
!$omp&  errnm_aux,rsdnm_aux,frc_aux,
!$omp&  global_proc_num_zones,global_proc_zone_id,thread_num)
!$omp&  if(global_nested.ne.2)

      call init_omp(num_zones, global_proc_zone_id,
     >              global_proc_num_zones)

      nthreads = global_proc_num_threads(global_myid+1)
c$    call omp_set_num_threads(nthreads)
!$    thread_num = omp_get_thread_num()

      write(*, 1000) thread_num
1000  format("I am thread ", i5)

      write (*, *)
      do i = 1, global_t_last
         call timer_clear(i)
      end do

      do iz = 1, global_proc_num_zones
        zone = global_proc_zone_id(iz)

c---------------------------------------------------------------------
c   set the boundary values for dependent variables
c---------------------------------------------------------------------
        call setbv(global_u(global_start5(zone)),
     $             global_nx(zone), global_nxmax(zone), global_ny(zone),
     >             global_nz(zone))

c---------------------------------------------------------------------
c   set the initial values for dependent variables
c---------------------------------------------------------------------
        call setiv(global_u(global_start5(zone)),
     $             global_nx(zone), global_nxmax(zone), global_ny(zone),
     >             global_nz(zone))

c---------------------------------------------------------------------
c   compute the forcing term based on prescribed exact solution
c---------------------------------------------------------------------
        call erhs(global_frct(global_start5(zone)),
     >            global_rsd(global_start5(zone)),
     $            global_nx(zone), global_nxmax(zone), global_ny(zone),
     >            global_nz(zone))

c---------------------------------------------------------------------
c   compute the steady-state residuals
c---------------------------------------------------------------------
        call rhs(global_u(global_start5(zone)),
     >           global_rsd(global_start5(zone)), 
     $           global_frct(global_start5(zone)),
     >           global_qs(global_start1(zone)), 
     $           global_rho_i(global_start1(zone)), 
     $           global_nx(zone), global_nxmax(zone), global_ny(zone),
     >           global_nz(zone))

      end do

c---------------------------------------------------------------------
c   initialize a,b,c,d to zero (guarantees that page tables have been
c   formed, if applicable on given architecture, before timestepping).
c   extra working arrays au, bu, cu, du are used in the OpenMP version
c   to align/touch data pages properly in the upper triangular solver.
c---------------------------------------------------------------------
      zone = global_proc_zone_id(1)
      call init_workarray(global_nx(zone), global_nxmax(zone),
     >                    global_ny(zone),
     $                    global_a, global_b, global_c, global_d,
     >                    global_au, global_bu, global_cu, global_du, 
     >                    global_tv)

c---------------------------------------------------------------------
c   perform one SSOR iteration to touch all data pages
c---------------------------------------------------------------------
      call exch_qbc(global_u, global_qbc, global_nx, global_nxmax,
     >              global_ny, global_nz, 
     &              global_proc_zone_id, global_proc_num_zones)

      do iz = 1, global_proc_num_zones
        zone = global_proc_zone_id(iz)
        call ssor(global_u(global_start5(zone)),
     >            global_rsd(global_start5(zone)), 
     $            global_frct(global_start5(zone)),
     >            global_qs(global_start1(zone)), 
     $            global_rho_i(global_start1(zone)), global_tv, 
     $            global_a, global_b, global_c, global_d, global_au,
     >            global_bu, global_cu, global_du, 
     $            global_nx(zone), global_nxmax(zone), global_ny(zone),
     >            global_nz(zone),
     &            global_isync)
      end do

c---------------------------------------------------------------------
c   reset the boundary and initial values
c---------------------------------------------------------------------
      do iz = 1, global_proc_num_zones
        zone = global_proc_zone_id(iz)

        call setbv(global_u(global_start5(zone)),
     $             global_nx(zone), global_nxmax(zone), global_ny(zone),
     >             global_nz(zone))

        call setiv(global_u(global_start5(zone)),
     $             global_nx(zone), global_nxmax(zone), global_ny(zone),
     >             global_nz(zone))

c---------------------------------------------------------------------
c   compute the steady-state residuals
c---------------------------------------------------------------------
        call rhs(global_u(global_start5(zone)),
     >           global_rsd(global_start5(zone)), 
     $           global_frct(global_start5(zone)),
     >           global_qs(global_start1(zone)), 
     $           global_rho_i(global_start1(zone)), 
     $           global_nx(zone), global_nxmax(zone), global_ny(zone),
     >           global_nz(zone))

      end do

c---------------------------------------------------------------------
c   begin pseudo-time stepping iterations
c---------------------------------------------------------------------

      do i = 1, global_t_last
         call timer_clear(i)
      end do

!$omp barrier
      call timer_start(1)

c---------------------------------------------------------------------
c   the timestep loop
c---------------------------------------------------------------------
      do step = 1, global_itmax

!$omp master
        if (mod(step,20) .eq. 0 .or. step .eq. 1 .or.
     >        step .eq. global_itmax) then
           write( *, 200) step
 200       format(' Time step ', i4)
        endif
!$omp end master

        call exch_qbc(global_u, global_qbc, global_nx, global_nxmax,
     >                global_ny, global_nz, 
     &                global_proc_zone_id, global_proc_num_zones)

c---------------------------------------------------------------------
c   perform the SSOR iterations
c---------------------------------------------------------------------

        do iz = 1, global_proc_num_zones
          zone = global_proc_zone_id(iz)
          call ssor(global_u(global_start5(zone)),
     >              global_rsd(global_start5(zone)), 
     $              global_frct(global_start5(zone)),
     >              global_qs(global_start1(zone)), 
     $              global_rho_i(global_start1(zone)), global_tv, 
     $              global_a, global_b, global_c, global_d, global_au,
     >              global_bu, global_cu, global_du, 
     $              global_nx(zone), global_nxmax(zone),
     >              global_ny(zone), global_nz(zone),
     &              global_isync)
        end do

      end do

!$omp master
      do i = 1, 5
         rsdnm(i) = 0.d0
         errnm(i) = 0.d0
      end do
      frc = 0.d0
!$omp end master
!$omp barrier

c---------------------------------------------------------------------
c   compute the max-norms of newton iteration residuals
c---------------------------------------------------------------------
      if (global_timeron) call timer_start(global_t_l2norm)
      do iz = 1, global_proc_num_zones
        zone = global_proc_zone_id(iz)
        call l2norm(global_rsd(global_start5(zone)), rsdnm_aux, 
     $              global_nx(zone), global_nxmax(zone),
     >              global_ny(zone), global_nz(zone))
        do i = 1, 5
!$omp atomic
          rsdnm(i) = rsdnm(i) + rsdnm_aux(i)
        end do
      end do

      if (global_timeron) call timer_stop(global_t_l2norm)

!$omp barrier
      call timer_stop(1)
      tmax = timer_read(1)

c---------------------------------------------------------------------
c   compute the solution error and surface integral
c---------------------------------------------------------------------
      do iz = 1, global_proc_num_zones
        zone = global_proc_zone_id(iz)
        call error(global_u(global_start5(zone)), errnm_aux,
     $             global_nx(zone), global_nxmax(zone), global_ny(zone),
     >             global_nz(zone))
        call pintgr(global_u(global_start5(zone)), global_phi1,
     >              global_phi2, frc_aux,
     $              global_nx(zone), global_nxmax(zone),
     >              global_ny(zone), global_nz(zone))
        do i = 1, 5
!$omp atomic
          errnm(i) = errnm(i) + errnm_aux(i)
        end do
!$omp atomic
        frc = frc + frc_aux
      end do


c---------------------------------------------------------------------
c   verification test
c---------------------------------------------------------------------
!$omp barrier
!$omp master
      call verify ( rsdnm, errnm, frc, verified )


      global_maxtime = tmax
      mflops = 0.d0

      if (global_maxtime .ne. 0.d0) then
        do zone = 1, num_zones
          n3 = dble(global_nx(zone))*global_ny(zone)*global_nz(zone)
          navg = (global_nx(zone) + global_ny(zone) +
     >            global_nz(zone))/3.d0

          nsur = (global_nx(zone)*global_ny(zone) +
     >            global_nx(zone)*global_nz(zone) +
     >            global_ny(zone)*global_nz(zone))/3.d0

          mflops = mflops + float(global_itmax)*1.0d-6 *
     >       (1984.77d0 * n3 - 10923.3d0 * nsur
     >         + 27770.9d0 * navg - 144010.d0)
     >       / global_maxtime
        end do
      endif

      call print_results('LU-MZ', global_class, global_gx_size,
     >                   global_gy_size, global_gz_size, 
     >                   global_itmax, global_maxtime, mflops, 
     >                   global_num_othreads, tot_threads,
     >                   '          floating point', verified, 
     >                   global_npbversion, global_compiletime, 
     >                   global_cs1, global_cs2, global_cs3, 
     >                   global_cs4, global_cs5, global_cs6, 
     >                   '(none)')

!$omp end master
!$omp barrier

c---------------------------------------------------------------------
c      More timers
c---------------------------------------------------------------------
c      if (.not.global_timeron) goto 999

      do i=1, global_t_last
         trecs(i) = timer_read(i)
      end do
      tmax = global_maxtime
      if (tmax .eq. 0.0) tmax = 1.0

      do i=1, global_t_last
!$omp atomic
         tsum(i) = tsum(i) + trecs(i)
!$omp atomic
         tming(i) = min(tming(i), trecs(i))
!$omp atomic
         tmaxg(i) = max(tmaxg(i), trecs(i))
      end do
!$omp barrier

!$omp master
!$    write(*, 700) global_num_othreads
!$    do i = 1, global_t_last
!$       tsum(i) = tsum(i) / global_num_othreads
c!$       write(*, 710) i, t_names(i), tming(i), tmaxg(i), tsum(i)
!$    end do
!$omp end master
 700  format(' #othrs =', i6, 11x, 'minimum', 5x, 'maximum', 
     >       5x, 'average')
 710  format(' timer ', i2, '(', A8, ') :', 3(2x,f10.4))

!$    if (itimer .lt. 2) goto 999

!$omp barrier
!$omp critical (ptime)
!$    write(*,800) global_myid, nthreads
 800  format(/' global_myid =',i5,'   num_ithreads =',i4)
      write(*,805)
 805  format('  SECTION   Time (secs)')
      do i=1, global_t_last
         write(*,810) t_names(i), trecs(i), trecs(i)*100./tmax
         if (i.eq.global_t_rhs) then
            t = trecs(global_t_rhsx) + trecs(global_t_rhsy) +
     >          trecs(global_t_rhsz)
            write(*,820) 'sub-rhs', t, t*100./tmax
            t = trecs(i) - t
            write(*,820) 'rest-rhs', t, t*100./tmax
         elseif (i.eq.global_t_rdis2) then
            t = trecs(global_t_rdis1) + trecs(global_t_rdis2)
            write(*,820) 'exch_qbc', t, t*100./tmax
         endif
 810     format(2x,a8,':',f9.3,'  (',f6.2,'%)')
 820     format(5x,'--> total ',a8,':',f9.3,'  (',f6.2,'%)')
      end do
!$omp end critical (ptime)

 999  continue

!$omp end parallel

      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine init_workarray(nx, nxmax, ny, a, b, c, d, 
     $                          au, bu, cu, du, tv)
      implicit none

c---------------------------------------------------------------------
c   initialize a,b,c,d to zero (guarantees that page tables have been
c   formed, if applicable on given architecture, before timestepping).
c   extra working arrays au, bu, cu, du are used in the OpenMP version
c   to align/touch data pages properly in the upper triangular solver.
c---------------------------------------------------------------------

      integer nx, nxmax, ny
      double precision a (25,2:nxmax-1,ny), b (25,2:nxmax-1,ny),
     $                 c (25,2:nxmax-1,ny), d (25,2:nxmax-1,ny),
     $                 au(25,2:nxmax-1,ny), bu(25,2:nxmax-1,ny),
     $                 cu(25,2:nxmax-1,ny), du(25,2:nxmax-1,ny),
     $                 tv( 5,2:nxmax-1,ny)

      integer i, j, m

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,i,j)
!$OMP DO SCHEDULE(STATIC)
      do j = 2, ny-1
        do i = 2, nx-1
          do m = 1, 25
            a(m,i,j) = 0.d0
            b(m,i,j) = 0.d0
            c(m,i,j) = 0.d0
            d(m,i,j) = 0.d0
          end do
        end do
      end do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
      do j = ny-1, 2, -1
        do i = nx-1, 2, -1
          do m = 1, 25
            au(m,i,j) = 0.d0
            bu(m,i,j) = 0.d0
            cu(m,i,j) = 0.d0
            du(m,i,j) = 0.d0
          end do
          do m = 1, 5
            tv(m,i,j) = 0.d0
          end do
        end do
      end do
!$OMP END DO nowait
!$OMP END PARALLEL

      return
      end
