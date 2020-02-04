c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  lu_data module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module lu_data

c---------------------------------------------------------------------
c   npbparams.h defines parameters that depend on the class and 
c   number of nodes
c---------------------------------------------------------------------

      include 'npbparams.h'

      integer(global_kind2), parameter :: 
     &        global_max_xysize=(global_gx_size+global_x_zones) * 
     >                           global_gy_size,
     >        global_proc_max_size=global_max_xysize*global_gz_size,
     &        global_proc_max_size5=global_proc_max_size*5,
     &        global_proc_max_bcsize=global_max_xysize*20

c---------------------------------------------------------------------
c   parameters which can be overridden in runtime config file
c   isiz1,isiz2,isiz3 give the maximum size
c   ipr = 1 to print out verbose information
c   omega = 2.0 is correct for all classes
c   tolrsd is tolerance levels for steady state residuals
c---------------------------------------------------------------------
      integer global_ipr_default
      parameter (global_ipr_default = 1)
      double precision global_omega_default
      parameter (global_omega_default = 1.2d0)
      double precision global_tolrsd1_def, global_tolrsd2_def,
     >                 global_tolrsd3_def, global_tolrsd4_def, 
     >                 global_tolrsd5_def
      parameter (global_tolrsd1_def=1.0e-08, 
     >          global_tolrsd2_def=1.0e-08, global_tolrsd3_def=1.0e-08,
     >          global_tolrsd4_def=1.0e-08, global_tolrsd5_def=1.0e-08)

      double precision global_c1, global_c2, global_c3, global_c4,
     >                 global_c5
      parameter( global_c1 = 1.40d+00, global_c2 = 0.40d+00,
     >           global_c3 = 1.00d-01, global_c4 = 1.00d+00,
     >           global_c5 = 1.40d+00 )

c---------------------------------------------------------------------
c   grid
c---------------------------------------------------------------------
      double precision  global_dxi, global_deta, global_dzeta
      double precision  global_tx1, global_tx2, global_tx3
      double precision  global_ty1, global_ty2, global_ty3
      double precision  global_tz1, global_tz2, global_tz3

c---------------------------------------------------------------------
c   dissipation
c---------------------------------------------------------------------
      double precision global_dx1, global_dx2, global_dx3, global_dx4,
     >                 global_dx5
      double precision global_dy1, global_dy2, global_dy3, global_dy4,
     >                 global_dy5
      double precision global_dz1, global_dz2, global_dz3, global_dz4,
     >                 global_dz5
      double precision global_dssp

      integer   global_max_zones
      parameter (global_max_zones=global_x_zones*global_y_zones)
      integer   global_x_start(global_x_zones),
     >          global_x_end(global_x_zones), 
     >          global_x_size(global_x_zones),
     >          global_y_start(global_y_zones),
     >          global_y_end(global_y_zones), 
     >          global_y_size(global_y_zones),
     >          global_iz_west (global_max_zones), 
     >          global_iz_east(global_max_zones),
     >          global_iz_south(global_max_zones),
     >          global_iz_north(global_max_zones)

      integer(global_kind2) :: global_start1(global_max_zones),
     >                         global_start5(global_max_zones),
     $                         global_qstart_west (global_max_zones),
     >                         global_qstart_east (global_max_zones),
     $                         global_qstart_south(global_max_zones),
     >                         global_qstart_north(global_max_zones)

c---------------------------------------------------------------------
c   output control parameters
c---------------------------------------------------------------------
      integer global_ipr, global_inorm, global_npb_verbose

c---------------------------------------------------------------------
c   newton-raphson iteration control parameters
c---------------------------------------------------------------------
      integer global_itmax, global_invert
      double precision  global_dt, global_omega, global_tolrsd(5),
     >                  global_ttotal

c---------------------------------------------------------------------
c   coefficients of the exact solution
c---------------------------------------------------------------------
      double precision global_ce(5,13)

c---------------------------------------------------------------------
c   1-d working arrays
c---------------------------------------------------------------------
      double precision  global_flux(5,global_problem_size), 
     >                  global_utmp(6,global_problem_size),
     >                  global_rtmp(5,global_problem_size)
!$omp threadprivate( global_flux, global_utmp, global_rtmp )

c---------------------------------------------------------------------
c   timers
c---------------------------------------------------------------------
      integer global_t_rhsx,global_t_rhsy,global_t_rhsz,global_t_rhs,
     >        global_t_jacld,global_t_blts,global_t_jacu,
     >        global_t_buts,global_t_add,global_t_l2norm,global_t_rdis1,
     >        global_t_rdis2,global_t_last,global_t_total
      parameter (global_t_total = 1)
      parameter (global_t_rhsx = 2)
      parameter (global_t_rhsy = 3)
      parameter (global_t_rhsz = 4)
      parameter (global_t_rhs = 5)
      parameter (global_t_jacld = 6)
      parameter (global_t_blts = 7)
      parameter (global_t_jacu = 8)
      parameter (global_t_buts = 9)
      parameter (global_t_add = 10)
      parameter (global_t_l2norm = 11)
      parameter (global_t_rdis1 = 12)
      parameter (global_t_rdis2 = 13)
      parameter (global_t_last = 13)
      logical global_timeron
      double precision global_maxtime

      end module lu_data


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  lu_fields module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module lu_fields

      use lu_data, only : global_max_zones, global_problem_size

      integer   global_nx(global_max_zones), 
     >          global_nxmax(global_max_zones),
     >          global_ny(global_max_zones), 
     $          global_nz(global_max_zones),
     &          global_proc_zone_id(global_max_zones), 
     >          global_proc_num_zones

      integer   global_isync(0:global_problem_size)

c---------------------------------------------------------------------
c   Define all field arrays as one-dimenional arrays, to be reshaped
c---------------------------------------------------------------------

      double precision , allocatable ::
     >                 global_u     (:),
     >                 global_rsd   (:),
     >                 global_frct  (:),
     >                 global_qs    (:),
     >                 global_rho_i (:),
     >                 global_qbc   (:)

c---------------------------------------------------------------------
c   2D auxiliary arrays are dimensioned to accommodate the largest
c   zone cross section
c---------------------------------------------------------------------

      double precision, allocatable ::
     $                 global_a   (:), 
     $                 global_b   (:), 
     $                 global_c   (:), 
     $                 global_d   (:),
     $                 global_au  (:), 
     $                 global_bu  (:), 
     $                 global_cu  (:), 
     $                 global_du  (:),
     $                 global_tv  (:),
     $                 global_phi1(:),
     $                 global_phi2(:)

      end module lu_fields


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  Allocate space for field arrays
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine alloc_field_space

      use lu_data, only :
     &            global_proc_max_size,
     &            global_proc_max_size5,
     &            global_proc_max_bcsize
      use lu_fields

      implicit none

      integer ios

      allocate (
     >          global_u     (global_proc_max_size5),
     >          global_rsd   (global_proc_max_size5),
     >          global_frct  (global_proc_max_size5),
     >          global_qs    (global_proc_max_size ),
     >          global_rho_i (global_proc_max_size ),
     >          global_qbc   (global_proc_max_bcsize),
     >          stat=ios )

      if (ios .eq. 0) allocate (
     $          global_a (25*global_problem_size*global_problem_size), 
     $          global_b (25*global_problem_size*global_problem_size), 
     $          global_c (25*global_problem_size*global_problem_size), 
     $          global_d (25*global_problem_size*global_problem_size),
     $          global_au(25*global_problem_size*global_problem_size), 
     $          global_bu(25*global_problem_size*global_problem_size), 
     $          global_cu(25*global_problem_size*global_problem_size), 
     $          global_du(25*global_problem_size*global_problem_size),
     $          global_tv( 5*global_problem_size*global_problem_size),
     $          global_phi1 (global_problem_size*global_problem_size),
     $          global_phi2 (global_problem_size*global_problem_size),
     >          stat=ios )

      if (ios .ne. 0) then
         write(*,*) 'Error encountered in allocating space'
         stop
      endif

      return
      end

