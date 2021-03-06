
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine rhs(u, rsd, frct, qs, rho_i, nx, nxmax, ny, nz)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c   compute the right hand sides
c---------------------------------------------------------------------

      use lu_data
      implicit none

      integer          nx, nxmax, ny, nz
      double precision u(5,nxmax,ny,nz), rsd(5,nxmax,ny,nz), 
     $                 frct(5,nxmax,ny,nz), qs(nxmax,ny,nz),
     $                 rho_i(nxmax,ny,nz)

c---------------------------------------------------------------------
c  local variables
c---------------------------------------------------------------------
      integer i, j, k, m
      double precision  q
      double precision  tmp
      double precision  u21, u31, u41
      double precision  u21i, u31i, u41i, u51i
      double precision  u21j, u31j, u41j, u51j
      double precision  u21k, u31k, u41k, u51k
      double precision  u21im1, u31im1, u41im1, u51im1
      double precision  u21jm1, u31jm1, u41jm1, u51jm1
      double precision  u21km1, u31km1, u41km1, u51km1


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(u51im1,u41im1,u31im1,u21im1,
!$OMP& u51i,u41i,u31i,u21i,tmp,m,q,u21,i,j,k,
!$OMP& u51jm1,u41jm1,u31jm1,u21jm1,u51j,u41j,u31j,u21j,u31,
!$OMP& u51km1,u41km1,u31km1,u21km1,u51k,u41k,u31k,u21k,u41)

      do k = 1, nz
!$OMP DO SCHEDULE(STATIC)
         do j = 1, ny
            do i = 1, nx
               do m = 1, 5
                  rsd(m,i,j,k) = - frct(m,i,j,k)
               end do
               tmp = 1.0d+00 / u(1,i,j,k)
               rho_i(i,j,k) = tmp
               qs(i,j,k) = 0.50d+00 * (  u(2,i,j,k) * u(2,i,j,k)
     >                         + u(3,i,j,k) * u(3,i,j,k)
     >                         + u(4,i,j,k) * u(4,i,j,k) )
     >                      * tmp
            end do
         end do
!$OMP END DO nowait
      end do
!$OMP BARRIER

!$OMP MASTER
      if (global_timeron) call timer_start(global_t_rhsx)
!$OMP END MASTER
c---------------------------------------------------------------------
c   xi-direction flux differences
c---------------------------------------------------------------------

      do k = 2, nz - 1
!$OMP DO SCHEDULE(STATIC)
         do j = 2, ny-1
            do i = 1, nx
               global_flux(1,i) = u(2,i,j,k)
               u21 = u(2,i,j,k) * rho_i(i,j,k)

               q = qs(i,j,k)

               global_flux(2,i) = u(2,i,j,k) * u21 + global_c2 * 
     >                        ( u(5,i,j,k) - q )
               global_flux(3,i) = u(3,i,j,k) * u21
               global_flux(4,i) = u(4,i,j,k) * u21
               global_flux(5,i) = ( global_c1 * u(5,i,j,k) - global_c2 *
     >q ) * u21
            end do

            do i = 2, nx-1
               do m = 1, 5
                  rsd(m,i,j,k) =  rsd(m,i,j,k)
     >                 - global_tx2 * ( global_flux(m,i+1) -
     >global_flux(m,i-1) )
               end do
            end do

            do i = 2, nx
               tmp = rho_i(i,j,k)

               u21i = tmp * u(2,i,j,k)
               u31i = tmp * u(3,i,j,k)
               u41i = tmp * u(4,i,j,k)
               u51i = tmp * u(5,i,j,k)

               tmp = rho_i(i-1,j,k)

               u21im1 = tmp * u(2,i-1,j,k)
               u31im1 = tmp * u(3,i-1,j,k)
               u41im1 = tmp * u(4,i-1,j,k)
               u51im1 = tmp * u(5,i-1,j,k)

               global_flux(2,i) = (4.0d+00/3.0d+00) * global_tx3 *
     >(u21i-u21im1)
               global_flux(3,i) = global_tx3 * ( u31i - u31im1 )
               global_flux(4,i) = global_tx3 * ( u41i - u41im1 )
               global_flux(5,i) = 0.50d+00 * ( 1.0d+00 -
     >global_c1*global_c5 )
     >              * global_tx3 * ( ( u21i  **2 + u31i  **2 + u41i  **2
     >)
     >                      - ( u21im1**2 + u31im1**2 + u41im1**2 ) )
     >              + (1.0d+00/6.0d+00)
     >              * global_tx3 * ( u21i**2 - u21im1**2 )
     >              + global_c1 * global_c5 * global_tx3 * ( u51i -
     >u51im1 )
            end do

            do i = 2, nx-1
               rsd(1,i,j,k) = rsd(1,i,j,k)
     >              + global_dx1 * global_tx1 * (           
     >u(1,i-1,j,k)
     >                             - 2.0d+00 * u(1,i,j,k)
     >                             +           u(1,i+1,j,k) )
               rsd(2,i,j,k) = rsd(2,i,j,k)
     >          + global_tx3 * global_c3 * global_c4 * (
     >global_flux(2,i+1) - global_flux(2,i) )
     >              + global_dx2 * global_tx1 * (           
     >u(2,i-1,j,k)
     >                             - 2.0d+00 * u(2,i,j,k)
     >                             +           u(2,i+1,j,k) )
               rsd(3,i,j,k) = rsd(3,i,j,k)
     >          + global_tx3 * global_c3 * global_c4 * (
     >global_flux(3,i+1) - global_flux(3,i) )
     >              + global_dx3 * global_tx1 * (           
     >u(3,i-1,j,k)
     >                             - 2.0d+00 * u(3,i,j,k)
     >                             +           u(3,i+1,j,k) )
               rsd(4,i,j,k) = rsd(4,i,j,k)
     >          + global_tx3 * global_c3 * global_c4 * (
     >global_flux(4,i+1) - global_flux(4,i) )
     >              + global_dx4 * global_tx1 * (           
     >u(4,i-1,j,k)
     >                             - 2.0d+00 * u(4,i,j,k)
     >                             +           u(4,i+1,j,k) )
               rsd(5,i,j,k) = rsd(5,i,j,k)
     >          + global_tx3 * global_c3 * global_c4 * (
     >global_flux(5,i+1) - global_flux(5,i) )
     >              + global_dx5 * global_tx1 * (           
     >u(5,i-1,j,k)
     >                             - 2.0d+00 * u(5,i,j,k)
     >                             +           u(5,i+1,j,k) )
            end do

c---------------------------------------------------------------------
c   Fourth-order dissipation
c---------------------------------------------------------------------
            do m = 1, 5
               rsd(m,2,j,k) = rsd(m,2,j,k)
     >           - global_dssp * ( + 5.0d+00 * u(m,2,j,k)
     >                      - 4.0d+00 * u(m,3,j,k)
     >                      +           u(m,4,j,k) )
               rsd(m,3,j,k) = rsd(m,3,j,k)
     >           - global_dssp * ( - 4.0d+00 * u(m,2,j,k)
     >                      + 6.0d+00 * u(m,3,j,k)
     >                      - 4.0d+00 * u(m,4,j,k)
     >                      +           u(m,5,j,k) )
            end do

            do i = 4, nx - 3
               do m = 1, 5
                  rsd(m,i,j,k) = rsd(m,i,j,k)
     >              - global_dssp * (            u(m,i-2,j,k)
     >                        - 4.0d+00 * u(m,i-1,j,k)
     >                        + 6.0d+00 * u(m,i,j,k)
     >                        - 4.0d+00 * u(m,i+1,j,k)
     >                        +           u(m,i+2,j,k) )
               end do
            end do


            do m = 1, 5
               rsd(m,nx-2,j,k) = rsd(m,nx-2,j,k)
     >           - global_dssp * (             u(m,nx-4,j,k)
     >                      - 4.0d+00 * u(m,nx-3,j,k)
     >                      + 6.0d+00 * u(m,nx-2,j,k)
     >                      - 4.0d+00 * u(m,nx-1,j,k)  )
               rsd(m,nx-1,j,k) = rsd(m,nx-1,j,k)
     >           - global_dssp * (             u(m,nx-3,j,k)
     >                      - 4.0d+00 * u(m,nx-2,j,k)
     >                      + 5.0d+00 * u(m,nx-1,j,k) )
            end do

         end do
!$OMP END DO nowait
      end do
!$OMP BARRIER

!$OMP MASTER
      if (global_timeron) call timer_stop(global_t_rhsx)

      if (global_timeron) call timer_start(global_t_rhsy)
!$OMP END MASTER
c---------------------------------------------------------------------
c   eta-direction flux differences
c---------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
      do k = 2, nz - 1
         do i = 2, nx-1
            do j = 1, ny
               global_utmp(1,j) = u(1,i,j,k)
               global_utmp(2,j) = u(2,i,j,k)
               global_utmp(3,j) = u(3,i,j,k)
               global_utmp(4,j) = u(4,i,j,k)
               global_utmp(5,j) = u(5,i,j,k)
               global_utmp(6,j) = rho_i(i,j,k)
            end do

            do j = 1, ny
               global_flux(1,j) = global_utmp(3,j)
               u31 = global_utmp(3,j) * global_utmp(6,j)

               q = qs(i,j,k)

               global_flux(2,j) = global_utmp(2,j) * u31 
               global_flux(3,j) = global_utmp(3,j) * u31 + global_c2 *
     >(global_utmp(5,j)-q)
               global_flux(4,j) = global_utmp(4,j) * u31
               global_flux(5,j) = ( global_c1 * global_utmp(5,j) -
     >global_c2 * q ) * u31
            end do

            do j = 2, ny-1
               do m = 1, 5
                  rsd(m,i,j,k) =  rsd(m,i,j,k)
     >                   - global_ty2 * ( global_flux(m,j+1) -
     >global_flux(m,j-1) )
               end do
            end do

            do j = 2, ny
               tmp = global_utmp(6,j)

               u21j = tmp * global_utmp(2,j)
               u31j = tmp * global_utmp(3,j)
               u41j = tmp * global_utmp(4,j)
               u51j = tmp * global_utmp(5,j)

               tmp = global_utmp(6,j-1)
               u21jm1 = tmp * global_utmp(2,j-1)
               u31jm1 = tmp * global_utmp(3,j-1)
               u41jm1 = tmp * global_utmp(4,j-1)
               u51jm1 = tmp * global_utmp(5,j-1)

               global_flux(2,j) = global_ty3 * ( u21j - u21jm1 )
               global_flux(3,j) = (4.0d+00/3.0d+00) * global_ty3 *
     >(u31j-u31jm1)
               global_flux(4,j) = global_ty3 * ( u41j - u41jm1 )
               global_flux(5,j) = 0.50d+00 * ( 1.0d+00 -
     >global_c1*global_c5 )
     >              * global_ty3 * ( ( u21j  **2 + u31j  **2 + u41j  **2
     >)
     >                      - ( u21jm1**2 + u31jm1**2 + u41jm1**2 ) )
     >              + (1.0d+00/6.0d+00)
     >              * global_ty3 * ( u31j**2 - u31jm1**2 )
     >              + global_c1 * global_c5 * global_ty3 * ( u51j -
     >u51jm1 )
            end do

            do j = 2, ny-1

               rsd(1,i,j,k) = rsd(1,i,j,k)
     >              + global_dy1 * global_ty1 * (           
     >global_utmp(1,j-1)
     >                             - 2.0d+00 * global_utmp(1,j)
     >                             +           global_utmp(1,j+1) )

               rsd(2,i,j,k) = rsd(2,i,j,k)
     >          + global_ty3 * global_c3 * global_c4 * (
     >global_flux(2,j+1) - global_flux(2,j) )
     >              + global_dy2 * global_ty1 * (           
     >global_utmp(2,j-1)
     >                             - 2.0d+00 * global_utmp(2,j)
     >                             +           global_utmp(2,j+1) )

               rsd(3,i,j,k) = rsd(3,i,j,k)
     >          + global_ty3 * global_c3 * global_c4 * (
     >global_flux(3,j+1) - global_flux(3,j) )
     >              + global_dy3 * global_ty1 * (           
     >global_utmp(3,j-1)
     >                             - 2.0d+00 * global_utmp(3,j)
     >                             +           global_utmp(3,j+1) )

               rsd(4,i,j,k) = rsd(4,i,j,k)
     >          + global_ty3 * global_c3 * global_c4 * (
     >global_flux(4,j+1) - global_flux(4,j) )
     >              + global_dy4 * global_ty1 * (           
     >global_utmp(4,j-1)
     >                             - 2.0d+00 * global_utmp(4,j)
     >                             +           global_utmp(4,j+1) )

               rsd(5,i,j,k) = rsd(5,i,j,k)
     >          + global_ty3 * global_c3 * global_c4 * (
     >global_flux(5,j+1) - global_flux(5,j) )
     >              + global_dy5 * global_ty1 * (           
     >global_utmp(5,j-1)
     >                             - 2.0d+00 * global_utmp(5,j)
     >                             +           global_utmp(5,j+1) )

            end do

c---------------------------------------------------------------------
c   fourth-order dissipation
c---------------------------------------------------------------------
            do m = 1, 5
               rsd(m,i,2,k) = rsd(m,i,2,k)
     >           - global_dssp * ( + 5.0d+00 * global_utmp(m,2)
     >                      - 4.0d+00 * global_utmp(m,3)
     >                      +           global_utmp(m,4) )
               rsd(m,i,3,k) = rsd(m,i,3,k)
     >           - global_dssp * ( - 4.0d+00 * global_utmp(m,2)
     >                      + 6.0d+00 * global_utmp(m,3)
     >                      - 4.0d+00 * global_utmp(m,4)
     >                      +           global_utmp(m,5) )
            end do

            do j = 4, ny - 3
               do m = 1, 5
                  rsd(m,i,j,k) = rsd(m,i,j,k)
     >              - global_dssp * (            global_utmp(m,j-2)
     >                        - 4.0d+00 * global_utmp(m,j-1)
     >                        + 6.0d+00 * global_utmp(m,j)
     >                        - 4.0d+00 * global_utmp(m,j+1)
     >                        +           global_utmp(m,j+2) )
               end do
            end do

            do m = 1, 5
               rsd(m,i,ny-2,k) = rsd(m,i,ny-2,k)
     >           - global_dssp * (             global_utmp(m,ny-4)
     >                      - 4.0d+00 * global_utmp(m,ny-3)
     >                      + 6.0d+00 * global_utmp(m,ny-2)
     >                      - 4.0d+00 * global_utmp(m,ny-1)  )
               rsd(m,i,ny-1,k) = rsd(m,i,ny-1,k)
     >           - global_dssp * (             global_utmp(m,ny-3)
     >                      - 4.0d+00 * global_utmp(m,ny-2)
     >                      + 5.0d+00 * global_utmp(m,ny-1) )
            end do

         end do
      end do
!$OMP END DO

!$OMP MASTER
      if (global_timeron) call timer_stop(global_t_rhsy)

      if (global_timeron) call timer_start(global_t_rhsz)
!$OMP END MASTER
c---------------------------------------------------------------------
c   zeta-direction flux differences
c---------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
      do j = 2, ny-1
         do i = 2, nx-1
            do k = 1, nz
               global_utmp(1,k) = u(1,i,j,k)
               global_utmp(2,k) = u(2,i,j,k)
               global_utmp(3,k) = u(3,i,j,k)
               global_utmp(4,k) = u(4,i,j,k)
               global_utmp(5,k) = u(5,i,j,k)
               global_utmp(6,k) = rho_i(i,j,k)
            end do

            do k = 1, nz
               global_flux(1,k) = global_utmp(4,k)
               u41 = global_utmp(4,k) * global_utmp(6,k)

               q = qs(i,j,k)

               global_flux(2,k) = global_utmp(2,k) * u41 
               global_flux(3,k) = global_utmp(3,k) * u41 
               global_flux(4,k) = global_utmp(4,k) * u41 + global_c2 *
     >(global_utmp(5,k)-q)
               global_flux(5,k) = ( global_c1 * global_utmp(5,k) -
     >global_c2 * q ) * u41
            end do

            do k = 2, nz - 1
               do m = 1, 5
                  global_rtmp(m,k) =  rsd(m,i,j,k)
     >                - global_tz2 * ( global_flux(m,k+1) -
     >global_flux(m,k-1) )
               end do
            end do

            do k = 2, nz
               tmp = global_utmp(6,k)

               u21k = tmp * global_utmp(2,k)
               u31k = tmp * global_utmp(3,k)
               u41k = tmp * global_utmp(4,k)
               u51k = tmp * global_utmp(5,k)

               tmp = global_utmp(6,k-1)

               u21km1 = tmp * global_utmp(2,k-1)
               u31km1 = tmp * global_utmp(3,k-1)
               u41km1 = tmp * global_utmp(4,k-1)
               u51km1 = tmp * global_utmp(5,k-1)

               global_flux(2,k) = global_tz3 * ( u21k - u21km1 )
               global_flux(3,k) = global_tz3 * ( u31k - u31km1 )
               global_flux(4,k) = (4.0d+00/3.0d+00) * global_tz3 *
     >(u41k-u41km1)
               global_flux(5,k) = 0.50d+00 * ( 1.0d+00 -
     >global_c1*global_c5 )
     >              * global_tz3 * ( ( u21k  **2 + u31k  **2 + u41k  **2
     >)
     >                      - ( u21km1**2 + u31km1**2 + u41km1**2 ) )
     >              + (1.0d+00/6.0d+00)
     >              * global_tz3 * ( u41k**2 - u41km1**2 )
     >              + global_c1 * global_c5 * global_tz3 * ( u51k -
     >u51km1 )
            end do

            do k = 2, nz - 1
               global_rtmp(1,k) = global_rtmp(1,k)
     >              + global_dz1 * global_tz1 * (           
     >global_utmp(1,k-1)
     >                             - 2.0d+00 * global_utmp(1,k)
     >                             +           global_utmp(1,k+1) )
               global_rtmp(2,k) = global_rtmp(2,k)
     >          + global_tz3 * global_c3 * global_c4 * (
     >global_flux(2,k+1) - global_flux(2,k) )
     >              + global_dz2 * global_tz1 * (           
     >global_utmp(2,k-1)
     >                             - 2.0d+00 * global_utmp(2,k)
     >                             +           global_utmp(2,k+1) )
               global_rtmp(3,k) = global_rtmp(3,k)
     >          + global_tz3 * global_c3 * global_c4 * (
     >global_flux(3,k+1) - global_flux(3,k) )
     >              + global_dz3 * global_tz1 * (           
     >global_utmp(3,k-1)
     >                             - 2.0d+00 * global_utmp(3,k)
     >                             +           global_utmp(3,k+1) )
               global_rtmp(4,k) = global_rtmp(4,k)
     >          + global_tz3 * global_c3 * global_c4 * (
     >global_flux(4,k+1) - global_flux(4,k) )
     >              + global_dz4 * global_tz1 * (           
     >global_utmp(4,k-1)
     >                             - 2.0d+00 * global_utmp(4,k)
     >                             +           global_utmp(4,k+1) )
               global_rtmp(5,k) = global_rtmp(5,k)
     >          + global_tz3 * global_c3 * global_c4 * (
     >global_flux(5,k+1) - global_flux(5,k) )
     >              + global_dz5 * global_tz1 * (           
     >global_utmp(5,k-1)
     >                             - 2.0d+00 * global_utmp(5,k)
     >                             +           global_utmp(5,k+1) )
            end do

c---------------------------------------------------------------------
c   fourth-order dissipation
c---------------------------------------------------------------------
            do m = 1, 5
               rsd(m,i,j,2) = global_rtmp(m,2)
     >           - global_dssp * ( + 5.0d+00 * global_utmp(m,2)
     >                      - 4.0d+00 * global_utmp(m,3)
     >                      +           global_utmp(m,4) )
               rsd(m,i,j,3) = global_rtmp(m,3)
     >           - global_dssp * ( - 4.0d+00 * global_utmp(m,2)
     >                      + 6.0d+00 * global_utmp(m,3)
     >                      - 4.0d+00 * global_utmp(m,4)
     >                      +           global_utmp(m,5) )
            end do

            do k = 4, nz - 3
               do m = 1, 5
                  rsd(m,i,j,k) = global_rtmp(m,k)
     >              - global_dssp * (            global_utmp(m,k-2)
     >                        - 4.0d+00 * global_utmp(m,k-1)
     >                        + 6.0d+00 * global_utmp(m,k)
     >                        - 4.0d+00 * global_utmp(m,k+1)
     >                        +           global_utmp(m,k+2) )
               end do
            end do

            do m = 1, 5
               rsd(m,i,j,nz-2) = global_rtmp(m,nz-2)
     >           - global_dssp * (             global_utmp(m,nz-4)
     >                      - 4.0d+00 * global_utmp(m,nz-3)
     >                      + 6.0d+00 * global_utmp(m,nz-2)
     >                      - 4.0d+00 * global_utmp(m,nz-1)  )
               rsd(m,i,j,nz-1) = global_rtmp(m,nz-1)
     >           - global_dssp * (             global_utmp(m,nz-3)
     >                      - 4.0d+00 * global_utmp(m,nz-2)
     >                      + 5.0d+00 * global_utmp(m,nz-1) )
            end do
         end do
      end do
!$OMP END DO nowait
!$OMP MASTER
      if (global_timeron) call timer_stop(global_t_rhsz)
!$OMP END MASTER
!$OMP END PARALLEL

      return
      end
