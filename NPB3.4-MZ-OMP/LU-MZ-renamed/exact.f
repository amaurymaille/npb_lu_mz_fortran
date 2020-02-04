c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine exact( i, j, k, u000ijk, nx, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c   compute the exact solution at (i,j,k)
c
c---------------------------------------------------------------------

      use lu_data
      implicit none

c---------------------------------------------------------------------
c  input parameters
c---------------------------------------------------------------------
      integer i, j, k, nx, ny, nz
      double precision u000ijk(*)

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer m
      double precision xi, eta, zeta

      xi  = ( dble ( i - 1 ) ) / ( nx - 1 )
      eta  = ( dble ( j - 1 ) ) / ( ny - 1 )
      zeta = ( dble ( k - 1 ) ) / ( nz - 1 )


      do m = 1, 5
         u000ijk(m) =  global_ce(m,1)
     >        + (global_ce(m,2)
     >        + (global_ce(m,5)
     >        + (global_ce(m,8)
     >        +  global_ce(m,11) * xi) * xi) * xi) * xi
     >        + (global_ce(m,3)
     >        + (global_ce(m,6)
     >        + (global_ce(m,9)
     >        +  global_ce(m,12) * eta) * eta) * eta) * eta
     >        + (global_ce(m,4)
     >        + (global_ce(m,7)
     >        + (global_ce(m,10)
     >        +  global_ce(m,13) * zeta) * zeta) * zeta) * zeta
      end do

      return
      end
