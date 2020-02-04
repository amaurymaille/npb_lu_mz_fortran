c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c  ompnpb module
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      module ompnpb
c
      use lu_data, only : global_max_zones
c
      integer   global_zone_proc_id(global_max_zones),        ! othread_id for >each zone
     &          global_proc_zone_count(global_max_zones),     ! #zones assigned >to othread
     &          global_proc_num_threads(global_max_zones),    ! #ithreads for >each othread
     &          global_proc_group(global_max_zones)           ! group_id for >each othread
      double precision global_proc_zone_size(global_max_zones)
c
      integer          global_myid, global_root, global_num_othreads,
     >                 global_num_threads, 
     &                 global_mz_bload, global_max_threads,
     >                 global_nested
!$omp threadprivate(global_myid, global_root)

      end module ompnpb

