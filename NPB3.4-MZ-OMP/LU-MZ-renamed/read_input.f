
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine read_input(tot_threads, itimer)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use lu_data
      use ompnpb

      implicit none

      integer tot_threads, itimer

      integer fstatus


      write(*, 1000) 

      call check_timer_flag( itimer )

      open (unit=2,file='inputlu-mz.data',status='old',
     >      access='sequential',form='formatted', iostat=fstatus)

      if (fstatus .eq. 0) then

         write(*,*) 'Reading from input file inputlu-mz.data'

         read (2,*)
         read (2,*)
         read (2,*) global_ipr, global_inorm
         read (2,*)
         read (2,*)
         read (2,*) global_itmax
         read (2,*)
         read (2,*)
         read (2,*) global_dt
         read (2,*)
         read (2,*)
         read (2,*) global_omega
         read (2,*)
         read (2,*)
         read (2,*)
     >global_tolrsd(1),global_tolrsd(2),global_tolrsd(3),
     >global_tolrsd(4),global_tolrsd(5)
         read (2,*,err=20,end=20)
         read (2,*,err=20,end=20)
         read (2,*,err=20,end=20) itimer
   20    close(2)

         if (global_itmax .eq. 0)  global_itmax = global_itmax_default
         if (global_dt .eq. 0.d0)  global_dt    = global_dt_default

      else
         global_ipr   = global_ipr_default
         global_inorm = global_inorm_default
         global_itmax = global_itmax_default
         global_dt    = global_dt_default
         global_omega = global_omega_default
         global_tolrsd(1) = global_tolrsd1_def
         global_tolrsd(2) = global_tolrsd2_def
         global_tolrsd(3) = global_tolrsd3_def
         global_tolrsd(4) = global_tolrsd4_def
         global_tolrsd(5) = global_tolrsd5_def
      endif

      write(*, 1001) global_x_zones, global_y_zones
      write(*, 1002) global_gx_size, global_gy_size, global_gz_size
      write(*, 1003) global_itmax, global_dt

 1000 format(//,' NAS Parallel Benchmarks (NPB3.4-MZ OpenMP)',
     >          ' - LU-MZ Benchmark', /)
 1001 format(' Number of zones: ', i3, ' x ', i3)
 1002 format(' Total mesh size: ', i5, ' x ', i5, ' x ', i3)
 1003 format(' Iterations: ', i3, '    global_dt: ', F10.6/)


      global_timeron = (itimer .gt. 0)
      call env_setup(tot_threads)

      return
      end


