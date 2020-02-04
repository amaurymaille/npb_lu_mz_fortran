
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine setcoeff

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      use lu_data
      implicit none

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c   set up coefficients
c---------------------------------------------------------------------

c     The following three settings are based on average cell size,
c     not actual cell size

      global_dxi   = 1.0d0 / ( dble(global_gx_size)/dble(global_x_zones)
     >- 1 )
      global_deta  = 1.0d0 / ( dble(global_gy_size)/dble(global_y_zones)
     >- 1 )
      global_dzeta = 1.0d0 / ( dble(global_gz_size - 1) )

      global_tx1 = 1.0d0 / ( global_dxi * global_dxi )
      global_tx2 = 1.0d0 / ( 2.0d0 * global_dxi )
      global_tx3 = 1.0d0 / global_dxi

      global_ty1 = 1.0d0 / ( global_deta * global_deta )
      global_ty2 = 1.0d0 / ( 2.0d0 * global_deta )
      global_ty3 = 1.0d0 / global_deta

      global_tz1 = 1.0d0 / ( global_dzeta * global_dzeta )
      global_tz2 = 1.0d0 / ( 2.0d0 * global_dzeta )
      global_tz3 = 1.0d0 / global_dzeta

c---------------------------------------------------------------------
c   diffusion coefficients
c---------------------------------------------------------------------
      global_dx1 = 0.75d0
      global_dx2 = global_dx1
      global_dx3 = global_dx1
      global_dx4 = global_dx1
      global_dx5 = global_dx1

      global_dy1 = 0.75d0
      global_dy2 = global_dy1
      global_dy3 = global_dy1
      global_dy4 = global_dy1
      global_dy5 = global_dy1

      global_dz1 = 1.00d0
      global_dz2 = global_dz1
      global_dz3 = global_dz1
      global_dz4 = global_dz1
      global_dz5 = global_dz1

c---------------------------------------------------------------------
c   fourth difference dissipation
c---------------------------------------------------------------------
      global_dssp = ( max (global_dx1, global_dy1, global_dz1 ) ) /
     >4.0d0

c---------------------------------------------------------------------
c   coefficients of the exact solution to the first pde
c---------------------------------------------------------------------
      global_ce(1,1) = 2.0d0
      global_ce(1,2) = 0.0d0
      global_ce(1,3) = 0.0d0
      global_ce(1,4) = 4.0d0
      global_ce(1,5) = 5.0d0
      global_ce(1,6) = 3.0d0
      global_ce(1,7) = 5.0d-1
      global_ce(1,8) = 2.0d-2
      global_ce(1,9) = 1.0d-2
      global_ce(1,10) = 3.0d-2
      global_ce(1,11) = 5.0d-1
      global_ce(1,12) = 4.0d-1
      global_ce(1,13) = 3.0d-1

c---------------------------------------------------------------------
c   coefficients of the exact solution to the second pde
c---------------------------------------------------------------------
      global_ce(2,1) = 1.0d0
      global_ce(2,2) = 0.0d0
      global_ce(2,3) = 0.0d0
      global_ce(2,4) = 0.0d0
      global_ce(2,5) = 1.0d0
      global_ce(2,6) = 2.0d0
      global_ce(2,7) = 3.0d0
      global_ce(2,8) = 1.0d-2
      global_ce(2,9) = 3.0d-2
      global_ce(2,10) = 2.0d-2
      global_ce(2,11) = 4.0d-1
      global_ce(2,12) = 3.0d-1
      global_ce(2,13) = 5.0d-1

c---------------------------------------------------------------------
c   coefficients of the exact solution to the third pde
c---------------------------------------------------------------------
      global_ce(3,1) = 2.0d0
      global_ce(3,2) = 2.0d0
      global_ce(3,3) = 0.0d0
      global_ce(3,4) = 0.0d0
      global_ce(3,5) = 0.0d0
      global_ce(3,6) = 2.0d0
      global_ce(3,7) = 3.0d0
      global_ce(3,8) = 4.0d-2
      global_ce(3,9) = 3.0d-2
      global_ce(3,10) = 5.0d-2
      global_ce(3,11) = 3.0d-1
      global_ce(3,12) = 5.0d-1
      global_ce(3,13) = 4.0d-1

c---------------------------------------------------------------------
c   coefficients of the exact solution to the fourth pde
c---------------------------------------------------------------------
      global_ce(4,1) = 2.0d0
      global_ce(4,2) = 2.0d0
      global_ce(4,3) = 0.0d0
      global_ce(4,4) = 0.0d0
      global_ce(4,5) = 0.0d0
      global_ce(4,6) = 2.0d0
      global_ce(4,7) = 3.0d0
      global_ce(4,8) = 3.0d-2
      global_ce(4,9) = 5.0d-2
      global_ce(4,10) = 4.0d-2
      global_ce(4,11) = 2.0d-1
      global_ce(4,12) = 1.0d-1
      global_ce(4,13) = 3.0d-1

c---------------------------------------------------------------------
c   coefficients of the exact solution to the fifth pde
c---------------------------------------------------------------------
      global_ce(5,1) = 5.0d0
      global_ce(5,2) = 4.0d0
      global_ce(5,3) = 3.0d0
      global_ce(5,4) = 2.0d0
      global_ce(5,5) = 1.0d-1
      global_ce(5,6) = 4.0d-1
      global_ce(5,7) = 3.0d-1
      global_ce(5,8) = 5.0d-2
      global_ce(5,9) = 4.0d-2
      global_ce(5,10) = 3.0d-2
      global_ce(5,11) = 1.0d-1
      global_ce(5,12) = 3.0d-1
      global_ce(5,13) = 2.0d-1

      return
      end


