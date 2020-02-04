       subroutine zone_setup(nx, nxmax, ny, nz)

       use lu_data
       use ompnpb

       implicit none

       integer nx(*), nxmax(*), ny(*), nz(*)

       integer           i,  j, zone_no
       integer           id_west, id_east, jd_south, jd_north
       double precision  x_r, y_r, x_smallest, y_smallest

       if (dabs(global_ratio-1.d0) .gt. 1.d-10) then

c        compute zone stretching only if the prescribed zone size ratio 
c        is substantially larger than unity       

         x_r   = dexp(dlog(global_ratio)/(global_x_zones-1))
         y_r   = dexp(dlog(global_ratio)/(global_y_zones-1))
         x_smallest = dble(global_gx_size)*(x_r-1.d0) / 
     >               (x_r**global_x_zones-1.d0)
         y_smallest = dble(global_gy_size)*(y_r-1.d0) / 
     >               (y_r**global_y_zones-1.d0)

c        compute tops of intervals, using a slightly tricked rounding
c        to make sure that the intervals are increasing monotonically
c        in size

         do i = 1, global_x_zones
            global_x_end(i) = x_smallest*(x_r**i-1.d0)/
     >                        (x_r-1.d0)+0.45d0
         end do

         do j = 1, global_y_zones
            global_y_end(j) = y_smallest*(y_r**j-1.d0)/
     >                        (y_r-1.d0)+0.45d0
         end do
 
       else

c        compute essentially equal sized zone dimensions

         do i = 1, global_x_zones
           global_x_end(i)   = (i*global_gx_size)/global_x_zones
         end do

         do j = 1, global_y_zones
           global_y_end(j)   = (j*global_gy_size)/global_y_zones
         end do

       endif

       global_x_start(1) = 1
       do i = 1, global_x_zones
          if (i .ne. global_x_zones) global_x_start(i+1) =
     >                               global_x_end(i) + 1
          global_x_size(i)  = global_x_end(i) - global_x_start(i) + 1
       end do

       global_y_start(1) = 1
       do j = 1, global_y_zones
          if (j .ne. global_y_zones) global_y_start(j+1) =
     >                               global_y_end(j) + 1
          global_y_size(j) = global_y_end(j) - global_y_start(j) + 1
       end do

       if (global_npb_verbose .gt. 1) write (*,98)
 98    format(/' Zone sizes:')

       do j = 1, global_y_zones
         do i = 1, global_x_zones
           zone_no = (i-1)+(j-1)*global_x_zones+1
           nx(zone_no) = global_x_size(i)
           nxmax(zone_no) = nx(zone_no) + 1 - mod(nx(zone_no),2)
           ny(zone_no) = global_y_size(j)
           nz(zone_no) = global_gz_size

           id_west  = mod(i-2+global_x_zones,global_x_zones)
           id_east  = mod(i,          global_x_zones)
           jd_south = mod(j-2+global_y_zones,global_y_zones)
           jd_north = mod(j,          global_y_zones)
           global_iz_west (zone_no) = id_west + (j-1)*global_x_zones + 1
           global_iz_east (zone_no) = id_east + (j-1)*global_x_zones + 1
           global_iz_south(zone_no) = (i-1) + jd_south*global_x_zones +
     >                                1
           global_iz_north(zone_no) = (i-1) + jd_north*global_x_zones +
     >                                1

           if (global_npb_verbose .gt. 1) then
             write (*,99) zone_no, nx(zone_no), ny(zone_no), 
     $                    nz(zone_no), nxmax(zone_no)
             write (*, 100) zone_no, global_iz_west(zone_no), 
     >                      global_iz_east(zone_no),
     >                      global_iz_south(zone_no),
     >                      global_iz_north(zone_no)
           endif
         end do
       end do

 99    format(i5,':  ',i5,' x',i5,' x',i5, ", nxmax: ", i5)
100    format("Zone ", i5, ": west = ", i5, ", east = ", i5, 
     >        ", south = ", i5, ", north = ", i5)

       return
       end


       subroutine zone_starts(num_zones, nx, nxmax, ny, nz)

       use lu_data
       use ompnpb

       implicit none

       integer   num_zones
       integer   nx(*), nxmax(*), ny(*), nz(*)

       integer   zone, zone_size
       integer   x_face_size, y_face_size

c ... index start for u & qbc
       do zone = 1, num_zones
          zone_size = nxmax(zone)*ny(zone)*nz(zone)
          x_face_size = (ny(zone)-2)*(nz(zone)-2)*5
          y_face_size = (nx(zone)-2)*(nz(zone)-2)*5

          if (zone .eq. 1) then
             global_qstart_west(zone) = 1
             global_start1(zone) = 1
             global_start5(zone) = 1
          endif
          global_qstart_east(zone)  = global_qstart_west(zone) +
     >                                x_face_size
          global_qstart_south(zone) = global_qstart_east(zone) +
     >                                x_face_size
          global_qstart_north(zone) = global_qstart_south(zone)+
     >                                y_face_size
          if (zone .ne. num_zones) then
             global_qstart_west(zone+1) = global_qstart_north(zone) +
     $                                    y_face_size
             global_start1(zone+1) = global_start1(zone) + zone_size
             global_start5(zone+1) = global_start5(zone) + zone_size*5
          else
             if (global_start1(zone) + zone_size-1 .gt.
     >           global_proc_max_size) then
                write(*,50) zone,global_proc_max_size,
     >                      global_start1(zone)+zone_size-1
                stop
             endif
          endif
   50     format(' Error in size: zone',i5,' global_proc_max_size',i10,
     &          ' access_size',i10)
       enddo

       if (global_npb_verbose .gt. 1) then
          do zone = 1, num_zones
             write(*,10) zone,global_start1(zone),global_start5(zone)
          enddo
       endif
   10  format(' zone=',i5,' global_start1=',i10,' global_start5=',i10)

       return
       end
