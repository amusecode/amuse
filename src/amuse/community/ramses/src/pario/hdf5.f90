! Copyright (c) 2007-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef HDF5

! TODO: que faire ds le cas ou la taille du tableau lu a change depuis
! l'ecriture? (ngridmax par exemple)

module ramses_hdf5
  use amr_commons,only:MPI_COMM_COMP,myid,pario_version,ramses_version
  use hdf5
  use mpi

  implicit none

  integer::status
  integer(hid_t)::mode_id, mode_indep_id

  interface hdf5_write_var
     module procedure hdf5_write_int,    hdf5_write_real4, &
                      hdf5_write_real8,  hdf5_write_real8_1D, &
                      hdf5_write_int_2D, hdf5_write_int8_2D, &
                      hdf5_write_char
  end interface

  interface hdf5_write_distrib_var
     module procedure hdf5_write_distrib_int,     hdf5_write_distrib_real4, &
                      hdf5_write_distrib_real8,   hdf5_write_distrib_int_1D, &
                      hdf5_write_distrib_int_2D,  hdf5_write_distrib_real8_1D, &
                      hdf5_write_distrib_real8_2D
  end interface

  interface hdf5_get_var
     module procedure hdf5_get_int,      hdf5_get_real8, &
                      hdf5_get_real8_1D, hdf5_get_int_2D, &
                      hdf5_get_int8_2D,  hdf5_get_char
  end interface

  interface hdf5_get_distrib_var
     module procedure hdf5_get_distrib_int,     hdf5_get_distrib_real8, &
                      hdf5_get_distrib_int_1D,  hdf5_get_distrib_int_2D, &
                      hdf5_get_distrib_real8_1D,hdf5_get_distrib_real8_2D
  end interface


  contains
    subroutine hdf5_init
      call h5open_f(status)
    end subroutine hdf5_init

    subroutine hdf5_finalize
      call h5close_f(status)
    end subroutine hdf5_finalize

    subroutine hdf5_create_file(fileloc,f_id)
      use amr_parameters,only:MAXLINE
      use hints,only:hints_set

      integer(hid_t),intent(out)::f_id
      character(LEN=*),intent(in)::fileloc

      integer(hid_t)::attr_id,grp_id,prp_id
      integer(hid_t)::comm_dt_id,comm_ds_id
      integer(hsize_t),dimension(1)::dim
      integer::info,ierr

      integer::majnum, minnum, relnum
      character(len=MAXLINE)::hdf5_version
      character(len=MAXLINE)::hostname
      character(len=8)::datec
      character(len=10)::timec

      call H5Pcreate_f(H5P_FILE_ACCESS_F,prp_id,status)
      call MPI_Info_create(info,status)
      !IBM_largeblock_io=true => bad performance for this code (x3 increaze in write time)
      !call MPI_Info_set(info,"IBM_largeblock_io","true",status)
      call hints_set(info)
      call H5Pset_fapl_mpio_f(prp_id,MPI_COMM_COMP,info,status)
      ! Using H5Pset_fapl_mpiposix_f instead of H5Pset_fapl_mpio_f seems to work
      ! but is extremely slow (equivalent to individual mode for H5Pset_fapl_mpio_f)
      !call H5Pset_fapl_mpiposix_f(prp_id,MPI_COMM_COMP,.true.,status)

      !H5Pset_alignment and H5Pset_cache: no effect on Blue Gene/P (20090825)
      !call H5Pset_alignment_f(prp_id,int(1024,kind=hsize_t),int(262144,kind=hsize_t),status)
      !call H5Pset_cache_f(prp_id,1000,int(10000,kind=size_t),int(16777216,kind=size_t),0.5,status)
      call H5Fcreate_f(fileloc,H5F_ACC_TRUNC_F,f_id,status,access_prp=prp_id)
      call H5Pclose_f(prp_id,status)
      call MPI_Info_free(info,ierr)

      ! Add version numbers as attributes
      call H5Gopen_f(f_id, "/",grp_id,status)

      call H5Tcopy_f(H5T_NATIVE_CHARACTER,comm_dt_id,status)
      call H5Tset_size_f(comm_dt_id,int(len(ramses_version),kind=size_t),status)
      call H5Screate_f(H5S_SCALAR_F,comm_ds_id,status)
      call H5Acreate_f(grp_id,'RAMSES version',comm_dt_id,comm_ds_id,attr_id,status)
      dim(1)=1
      call H5Awrite_f(attr_id,comm_dt_id,ramses_version,dim,status)
      call H5Aclose_f(attr_id,status)
      call H5Sclose_f(comm_ds_id,status)
      call H5Tclose_f(comm_dt_id,status)

      call H5Tcopy_f(H5T_NATIVE_CHARACTER,comm_dt_id,status)
      call H5Tset_size_f(comm_dt_id,int(len(pario_version),kind=size_t),status)
      call H5Screate_f(H5S_SCALAR_F,comm_ds_id,status)
      call H5Acreate_f(grp_id,'Parallel IO version',comm_dt_id,comm_ds_id,attr_id,status)
      dim(1)=1
      call H5Awrite_f(attr_id,comm_dt_id,pario_version,dim,status)
      call H5Aclose_f(attr_id,status)
      call H5Sclose_f(comm_ds_id,status)
      call H5Tclose_f(comm_dt_id,status)

      call H5get_libversion_f(majnum,minnum,relnum,status)
      dim(1)=3
      call H5Screate_simple_f(1,dim,comm_ds_id,status)
      call H5Acreate_f(grp_id,'HDF5 version',H5T_NATIVE_INTEGER,comm_ds_id,attr_id,status)
      call H5Awrite_f(attr_id,H5T_NATIVE_INTEGER,(/majnum,minnum,relnum/),dim,status)
      call H5Aclose_f(attr_id,status)
      call H5Sclose_f(comm_ds_id,status)

      ! Add hostname as attribute
      hostname=' ' ! Necessary to have a valid Fortran string
      call gethname(hostname)
      call H5Tcopy_f(H5T_NATIVE_CHARACTER,comm_dt_id,status)
      call H5Tset_size_f(comm_dt_id,int(len_trim(hostname),kind=size_t),status)
      call H5Screate_f(H5S_SCALAR_F,comm_ds_id,status)
      call H5Acreate_f(grp_id,'hostname',comm_dt_id,comm_ds_id,attr_id,status)
      dim(1)=1
      call H5Awrite_f(attr_id,comm_dt_id,trim(hostname),dim,status)
      call H5Aclose_f(attr_id,status)
      call H5Sclose_f(comm_ds_id,status)
      call H5Tclose_f(comm_dt_id,status)

      ! Add date and time as attributes
      call date_and_time(date=datec,time=timec)
      call H5Tcopy_f(H5T_NATIVE_CHARACTER,comm_dt_id,status)
      call H5Tset_size_f(comm_dt_id,int(len(datec),kind=size_t),status)
      call H5Screate_f(H5S_SCALAR_F,comm_ds_id,status)
      call H5Acreate_f(grp_id,'date',comm_dt_id,comm_ds_id,attr_id,status)
      dim(1)=1
      call H5Awrite_f(attr_id,comm_dt_id,datec,dim,status)
      call H5Aclose_f(attr_id,status)
      call H5Sclose_f(comm_ds_id,status)
      call H5Tclose_f(comm_dt_id,status)

      call H5Tcopy_f(H5T_NATIVE_CHARACTER,comm_dt_id,status)
      call H5Tset_size_f(comm_dt_id,int(len(timec),kind=size_t),status)
      call H5Screate_f(H5S_SCALAR_F,comm_ds_id,status)
      call H5Acreate_f(grp_id,'time',comm_dt_id,comm_ds_id,attr_id,status)
      dim(1)=1
      call H5Awrite_f(attr_id,comm_dt_id,timec,dim,status)
      call H5Aclose_f(attr_id,status)
      call H5Sclose_f(comm_ds_id,status)
      call H5Tclose_f(comm_dt_id,status)


      call H5Gclose_f(grp_id,status)


      !Create property list for access mode (used when switching between
      !independent and collective MPI-I/O modes
      call H5Pcreate_f(H5P_DATASET_XFER_F,mode_id,status) 
      call H5Pcreate_f(H5P_DATASET_XFER_F,mode_indep_id,status)

      !Put mode_id in collective mode
      call H5Pset_dxpl_mpio_f(mode_id,H5FD_MPIO_COLLECTIVE_F,status)
      !Put mode_indep_id in independent mode
      call H5Pset_dxpl_mpio_f(mode_indep_id,H5FD_MPIO_INDEPENDENT_F,status)

    end subroutine hdf5_create_file

    subroutine hdf5_open_file(fileloc,f_id)
      use hints,only:hints_set

      integer(hid_t),intent(out)::f_id
      character(LEN=*),intent(in)::fileloc

      integer :: ierr,info
      integer(hid_t)::prp_id

      call H5Pcreate_f(H5P_FILE_ACCESS_F,prp_id,status)
      call MPI_Info_create(info,status)
      call hints_set(info)
      call H5Pset_fapl_mpio_f(prp_id,MPI_COMM_COMP,info,status)
      ! Using H5Pset_fapl_mpiposix_f instead of H5Pset_fapl_mpio_f seems to work
      ! but is extremely slow (equivalent to individual mode for H5Pset_fapl_mpio_f)
      !call H5Pset_fapl_mpiposix_f(prp_id,MPI_COMM_COMP,.true.,status)
      call H5Fopen_f(fileloc,H5F_ACC_RDONLY_F,f_id,status,access_prp=prp_id)
      call H5Pclose_f(prp_id,status)
      call MPI_Info_free(info,ierr)

      !Create property list for access mode (used when switching between
      !independent and collective MPI-I/O modes
      call H5Pcreate_f(H5P_DATASET_XFER_F,mode_id,status) 

      !Put mode_id in collective mode
      call H5Pset_dxpl_mpio_f(mode_id,H5FD_MPIO_COLLECTIVE_F,status)

    end subroutine hdf5_open_file

    subroutine hdf5_close_file(f_id)
      integer(hid_t)::f_id

      call H5Pclose_f(mode_id,status)
      call H5Fclose_f(f_id,status)
    end subroutine hdf5_close_file

    subroutine hdf5_write_int(f_id,name,var)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,intent(in)::var

      integer(hid_t)::dset_id,dsp_id,prp_id
      integer(hsize_t),dimension(1)::dim

      dim(1)=1
      call H5Pcreate_f(H5P_DATASET_CREATE_F,prp_id,status)
      !Compact storage does not work if only 1 process writes (tested on lin2 with hdf5 1.8.7)
      !call H5Pset_layout_f(prp_id,H5D_COMPACT_F,status) !Use compact storage (faster)
      call H5Screate_simple_f(1,dim,dsp_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,dsp_id,dset_id,status,prp_id)
      ! Only process 1 do the write but the previous calls are collective (especially H5Dcreate_f)
      if(myid==1) call H5Dwrite_f(dset_id,H5T_NATIVE_INTEGER,var,dim,status,xfer_prp=mode_indep_id)
      call H5Dclose_f(dset_id,status)
      call H5Sclose_f(dsp_id,status)
      call H5Pclose_f(prp_id,status)

    end subroutine hdf5_write_int

    subroutine hdf5_write_real4(f_id,name,var)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=4),intent(in)::var

      integer(hid_t)::dset_id,dsp_id,prp_id
      integer(hsize_t),dimension(1)::dim

      dim(1)=1
      !call H5Pcreate_f(H5P_DATASET_CREATE_F,prp_id,status)
      !Compact storage does not work if only 1 process writes (tested on lin2 with hdf5 1.8.7)
      !call H5Pset_layout_f(prp_id,H5D_COMPACT_F,status) !Use compact storage (faster)
      call H5Screate_simple_f(1,dim,dsp_id,status)
      !call H5Dcreate_f(f_id,name,H5T_NATIVE_REAL,dsp_id,dset_id,status,prp_id)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_REAL,dsp_id,dset_id,status)
      ! Only process 1 do the write but the previous calls are collective (especially H5Dcreate_f)
      if(myid==1) call H5Dwrite_f(dset_id,H5T_NATIVE_REAL,var,dim,status,xfer_prp=mode_indep_id)
      call H5Dclose_f(dset_id,status)
      call H5Sclose_f(dsp_id,status)
      !call H5Pclose_f(prp_id,status)

    end subroutine hdf5_write_real4

    subroutine hdf5_write_real8(f_id,name,var)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),intent(in)::var

      integer(hid_t)::dset_id,dsp_id
      integer(hsize_t),dimension(1)::dim

      dim(1)=1
      !call H5Pcreate_f(H5P_DATASET_CREATE_F,prp_id,status)
      !Compact storage does not work if only 1 process writes (tested on lin2 with hdf5 1.8.7)
      !call H5Pset_layout_f(prp_id,H5D_COMPACT_F,status) !Use compact storage (faster)
      call H5Screate_simple_f(1,dim,dsp_id,status)
      !call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,dsp_id,dset_id,status,prp_id)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,dsp_id,dset_id,status)
      ! Only process 1 do the write but the previous calls are collective (especially H5Dcreate_f)
      if(myid==1) call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,dim,status,xfer_prp=mode_indep_id)
      call H5Dclose_f(dset_id,status)
      call H5Sclose_f(dsp_id,status)
      !call H5Pclose_f(prp_id,status)

    end subroutine hdf5_write_real8

    subroutine hdf5_write_real8_1D(f_id,name,var,count)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(in)::var
      integer,optional,intent(in)::count

      integer(hid_t)::dset_id,dsp_id
      integer(hsize_t),dimension(1)::dim

      if(present(count))then
         dim(1)=count
      else
         dim(1)=size(var,1)
      end if

      call H5Screate_simple_f(1,dim,dsp_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,dsp_id,dset_id,status)
      ! Only process 1 do the write but the previous calls are collective (especially H5Dcreate_f)
      if(myid==1) call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,dim,status,xfer_prp=mode_indep_id)
      call H5Dclose_f(dset_id,status)
      call H5Sclose_f(dsp_id,status)

    end subroutine hdf5_write_real8_1D

    subroutine hdf5_write_int_2D(f_id,name,var,count)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:,:),intent(in)::var
      integer,dimension(2),optional,intent(in)::count

      integer(hid_t)::dset_id,dsp_id
      integer(hsize_t),dimension(2)::dim

      if(present(count))then
         dim(1)=count(1)
         dim(2)=count(2)
      else
         dim(1)=size(var,1)
         dim(2)=size(var,2)
      end if

      call H5Screate_simple_f(2,dim,dsp_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,dsp_id,dset_id,status)
      ! Only process 1 do the write but the previous calls are collective (especially H5Dcreate_f)
      if(myid==1) call H5Dwrite_f(dset_id,H5T_NATIVE_INTEGER,var,dim,status,xfer_prp=mode_indep_id)
      call H5Dclose_f(dset_id,status)
      call H5Sclose_f(dsp_id,status)

    end subroutine hdf5_write_int_2D

    subroutine hdf5_write_int8_2D(f_id,name,var,count)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer(kind=8),dimension(:,:),intent(in)::var
      integer,dimension(2),optional,intent(in)::count

      integer(hid_t)::dset_id,dsp_id,int8_id
      integer(hsize_t),dimension(2)::dim

      if(present(count))then
         dim(1)=count(1)
         dim(2)=count(2)
      else
         dim(1)=size(var,1)
         dim(2)=size(var,2)
      end if

      !TODO:write 64 bit integers (do not seem to work for the moment?)
      !Create int*8 datatype
      call H5Tcopy_f(H5T_NATIVE_INTEGER,int8_id,status)
      call H5Tset_precision_f(int8_id,int(64,kind=size_t),status)

      call H5Screate_simple_f(2,dim,dsp_id,status)
      !call H5Dcreate_f(f_id,name,int8_id,dsp_id,dset_id,status)
      !call H5Dwrite_f(dset_id,int8_id,var,dim,status,xfer_prp=mode_indep_id)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,dsp_id,dset_id,status)
      ! Only process 1 do the write but the previous calls are collective (especially H5Dcreate_f)
      if(myid==1) call H5Dwrite_f(dset_id,H5T_NATIVE_DOUBLE,reshape(transfer(var(:,:),1.0_8,dim(1)*dim(2)),shape(var)),&
                      dim,status,xfer_prp=mode_indep_id)
      call H5Dclose_f(dset_id,status)
      call H5Sclose_f(dsp_id,status)
      !call H5Tclose_f(int8_id,status)

    end subroutine hdf5_write_int8_2D

    subroutine hdf5_write_char(f_id,name,var)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      character(len=*),intent(in)::var

      integer(hid_t)::dset_id,dsp_id
      integer(hsize_t),dimension(1)::dim

      dim(1)=len(var)
      !call H5Pcreate_f(H5P_DATASET_CREATE_F,prp_id,status)
      !Compact storage does not work if only 1 process writes (tested on lin2 with hdf5 1.8.7)
      !call H5Pset_layout_f(prp_id,H5D_COMPACT_F,status) !Use compact storage (faster)
      call H5Screate_simple_f(1,dim,dsp_id,status)
      !call H5Dcreate_f(f_id,name,H5T_NATIVE_CHARACTER,dsp_id,dset_id,status,prp_id)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_CHARACTER,dsp_id,dset_id,status)
      ! Only process 1 do the write but the previous calls are collective (especially H5Dcreate_f)
      if(myid==1) call H5Dwrite_f(dset_id,H5T_NATIVE_CHARACTER,var,dim,status,xfer_prp=mode_indep_id)
      call H5Dclose_f(dset_id,status)
      call H5Sclose_f(dsp_id,status)
      !call H5Pclose_f(prp_id,status)

    end subroutine hdf5_write_char



    subroutine hdf5_write_distrib_int(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,intent(in)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::dims,starts,counts

      dims(1)=1
      dims(2)=ncpu
      call H5Screate_simple_f(1,dims,ds_id,status) !1D in memory
      call H5Screate_simple_f(2,dims,fs_id,status) !2D in file
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(1)=dims(1)
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,var,dims,status,file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_distrib_int

    subroutine hdf5_write_distrib_real4(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=4),intent(in)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::dims,starts,counts

      dims(1)=1
      dims(2)=ncpu
      call H5Screate_simple_f(1,dims,ds_id,status) !1D in memory
      call H5Screate_simple_f(2,dims,fs_id,status) !2D in file
      call H5Dcreate_f(f_id,name,H5T_NATIVE_REAL,fs_id,da_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(1)=dims(1)
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_REAL,var,dims,status,file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_distrib_real4

    subroutine hdf5_write_distrib_real8(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),intent(in)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::dims,starts,counts

      dims(1)=1
      dims(2)=ncpu
      call H5Screate_simple_f(1,dims,ds_id,status) !1D in memory
      call H5Screate_simple_f(2,dims,fs_id,status) !2D in file
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,fs_id,da_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(1)=dims(1)
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,var,dims,status,file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_distrib_real8

    subroutine hdf5_write_distrib_int_1D(f_id,name,var,count)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(in)::var
      integer,optional,intent(in)::count

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::dims,starts,counts

      if(present(count))then
         dims(1)=count
      else
          dims(1)=size(var,1)
     end if
      dims(2)=ncpu
      call H5Screate_simple_f(1,dims,ds_id,status) !1D in memory
      call H5Screate_simple_f(2,dims,fs_id,status) !2D in file
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(1)=dims(1)
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,var,dims,status,file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_distrib_int_1D

    subroutine hdf5_write_distrib_int_2D(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:,:),intent(in)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(3)::dims,starts,counts

      dims(1)=size(var,1)
      dims(2)=size(var,2)
      dims(3)=ncpu
      call H5Screate_simple_f(2,dims,ds_id,status) !2D in memory
      call H5Screate_simple_f(3,dims,fs_id,status) !3D in file
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)
      starts(1)=0
      starts(2)=0
      starts(3)=myid-1
      counts(1)=dims(1)
      counts(2)=dims(2)
      counts(3)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,var,dims,status,file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_distrib_int_2D

    subroutine hdf5_write_distrib_real8_1D(f_id,name,var,count)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(in)::var
      integer,optional,intent(in)::count

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::dims,starts,counts

      if(present(count))then
         dims(1)=count
      else
          dims(1)=size(var,1)
      end if
      dims(2)=ncpu
      call H5Screate_simple_f(1,dims,ds_id,status) !1D in memory
      call H5Screate_simple_f(2,dims,fs_id,status) !2D in file
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,fs_id,da_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(1)=dims(1)
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,var,dims,status,file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_distrib_real8_1D

    subroutine hdf5_write_distrib_real8_2D(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(in)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(3)::dims,starts,counts

      dims(1)=size(var,1)
      dims(2)=size(var,2)
      dims(3)=ncpu
      call H5Screate_simple_f(2,dims,ds_id,status) !2D in memory
      call H5Screate_simple_f(3,dims,fs_id,status) !3D in file
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,fs_id,da_id,status)
      starts(1)=0
      starts(2)=0
      starts(3)=myid-1
      counts(1)=dims(1)
      counts(2)=dims(2)
      counts(3)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,var,dims,status,file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_distrib_real8_2D

    subroutine hdf5_get_int(f_id,name,var)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,intent(out)::var

      integer(hid_t)::dset_id
      integer(hsize_t),dimension(1)::dim

      dim(1)=1
      call H5Dopen_f(f_id,name,dset_id,status)
      call H5Dread_f(dset_id,H5T_NATIVE_INTEGER,var,dim,status,xfer_prp=mode_id)
      call H5Dclose_f(dset_id,status)
    end subroutine hdf5_get_int

    subroutine hdf5_get_real8(f_id,name,var)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),intent(out)::var

      integer(hid_t)::dset_id
      integer(hsize_t),dimension(1)::dim

      dim(1)=1
      call H5Dopen_f(f_id,name,dset_id,status)
      call H5Dread_f(dset_id,H5T_NATIVE_DOUBLE,var,dim,status,xfer_prp=mode_id)
      call H5Dclose_f(dset_id,status)
    end subroutine hdf5_get_real8

    subroutine hdf5_get_real8_1D(f_id,name,var,count)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(out)::var
      integer,optional,intent(in)::count

      integer(hid_t)::dset_id
      integer(hsize_t),dimension(1)::dim

       if(present(count))then
         dim(1)=count
      else
         dim(1)=size(var,1)
      end if
      call H5Dopen_f(f_id,name,dset_id,status)
      call H5Dread_f(dset_id,H5T_NATIVE_DOUBLE,var,dim,status,xfer_prp=mode_id)
      call H5Dclose_f(dset_id,status)
    end subroutine hdf5_get_real8_1D

    subroutine hdf5_get_int_2D(f_id,name,var,counts)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:,:),intent(out)::var
      integer,optional,dimension(2),intent(in)::counts

      integer(hid_t)::dset_id
      integer(hsize_t),dimension(2)::dims

       if(present(counts))then
         dims(1)=counts(1)
         dims(2)=counts(2)
      else
         dims(1)=size(var,1)
         dims(2)=size(var,2)
      end if
      call H5Dopen_f(f_id,name,dset_id,status)
      call H5Dread_f(dset_id,H5T_NATIVE_INTEGER,var,dims,status,xfer_prp=mode_id)
      call H5Dclose_f(dset_id,status)
    end subroutine hdf5_get_int_2D

    subroutine hdf5_get_int8_2D(f_id,name,var,counts)
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer(kind=8),dimension(:,:),intent(out)::var
      integer,optional,dimension(2),intent(in)::counts

      integer(hid_t)::dset_id
      integer(hsize_t),dimension(2)::dims
      real(kind=8),dimension(:,:),allocatable::var_dbl

       if(present(counts))then
         dims(1)=counts(1)
         dims(2)=counts(2)
      else
         dims(1)=size(var,1)
         dims(2)=size(var,2)
      end if

      allocate(var_dbl(dims(1),dims(2)))

      call H5Dopen_f(f_id,name,dset_id,status)
      call H5Dread_f(dset_id,H5T_NATIVE_DOUBLE,var_dbl,dims,status,xfer_prp=mode_id)
      call H5Dclose_f(dset_id,status)

      var = reshape(transfer(var_dbl,1_8,size(var)),shape(var))

      deallocate(var_dbl)
    end subroutine hdf5_get_int8_2D

    subroutine hdf5_get_char(f_id,name,var)
      !TODO: secure length
      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      character(len=*),intent(out)::var

      integer(hid_t)::dset_id
      integer(hsize_t),dimension(1)::dim

      dim(1)=len(var)
      call H5Dopen_f(f_id,name,dset_id,status)
      call H5Dread_f(dset_id,H5T_NATIVE_CHARACTER,var,dim,status,xfer_prp=mode_id)
      call H5Dclose_f(dset_id,status)
    end subroutine hdf5_get_char



    subroutine hdf5_get_distrib_int(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,intent(out)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::starts,counts
      integer(hsize_t),dimension(1)::dims

      dims(1)=1
      counts(1)=1
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_INTEGER,var,dims,status,&
                     file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_get_distrib_int

    subroutine hdf5_get_distrib_real8(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),intent(out)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::starts,counts
      integer(hsize_t),dimension(1)::dims

      dims(1)=1
      counts(1)=1
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_DOUBLE,var,dims,status,&
                     file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_get_distrib_real8

    subroutine hdf5_get_distrib_int_1D(f_id,name,var,count)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(out)::var
      integer,optional,intent(in)::count

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::starts,counts
      integer(hsize_t),dimension(1)::dims

      if(present(count))then
         dims(1)=count
         counts(1)=count
      else
         dims(1)=size(var,1)
         counts(1)=size(var,1)
      end if
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_INTEGER,var,dims,status,&
                     file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_get_distrib_int_1D

    subroutine hdf5_get_distrib_int_2D(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:,:),intent(out)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(3)::starts,counts
      integer(hsize_t),dimension(2)::dims

      dims(1)=size(var,1)
      dims(2)=size(var,2)
      call H5Screate_simple_f(2,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=0
      starts(2)=0
      starts(3)=myid-1
      counts(1)=size(var,1)
      counts(2)=size(var,2)
      counts(3)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_INTEGER,var,dims,status,&
                     mem_space_id=ds_id,file_space_id=fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_get_distrib_int_2D

    subroutine hdf5_get_distrib_real8_1D(f_id,name,var,count)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(out)::var
      integer,optional,intent(in)::count

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(2)::starts,counts
      integer(hsize_t),dimension(1)::dims

      if(present(count))then
         dims(1)=count
         counts(1)=count
      else
         dims(1)=size(var,1)
         counts(1)=size(var,1)
      end if
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=0
      starts(2)=myid-1
      counts(2)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_DOUBLE,var,dims,status,&
                     file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_get_distrib_real8_1D

    subroutine hdf5_get_distrib_real8_2D(f_id,name,var)
      use amr_commons,only:myid,ncpu

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(out)::var

      integer(hid_t)::da_id,ds_id,fs_id
      integer(hsize_t),dimension(3)::starts,counts
      integer(hsize_t),dimension(2)::dims

      dims(1)=size(var,1)
      dims(2)=size(var,2)
      call H5Screate_simple_f(2,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=0
      starts(2)=0
      starts(3)=myid-1
      counts(1)=size(var,1)
      counts(2)=size(var,2)
      counts(3)=1
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_DOUBLE,var,dims,status,&
                     mem_space_id=ds_id,file_space_id=fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_get_distrib_real8_2D

    subroutine hdf5_write_multilevel_int_1D(f_id,name,var)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(in)::var

      integer(kind=mpi_offset_kind)::sz
      integer(hsize_t),dimension(1)::dims,starts,counts,start_ds
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,idx_old,igrid,ilevel,ibound,ncache,istart,i
      integer,dimension(buffer_size)::iig

      ! Create memory dataspace
      dims(1)=buffer_size
      call H5Screate_simple_f(1,dims,ds_id,status)

      ! Create file dataspace
      dims(1)=size_max
      call H5Screate_simple_f(1,dims,fs_id,status)

      ! Create dataset
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)

      idx=1
      sz=0
      start_ds(1)=0
      call H5Sselect_none_f(fs_id,status)
      ! Loop on the refinement levels
      do ilevel=1,nlevelmax
         idx_old=idx
         do ibound=1,nboundary+ncpu
            if(ibound<=ncpu)then
               ncache=numbl(ibound,ilevel)
               istart=headl(ibound,ilevel)
            else
               ncache=numbb(ibound-ncpu,ilevel)
               istart=headb(ibound-ncpu,ilevel)
            end if
            if(ncache>0)then
               igrid=istart
               do i=idx,idx-1+ncache
                  iig(i)=var(igrid)
                  igrid=next(igrid)
               end do
               sz=sz+ncache
               idx=idx+ncache
            end if
        end do
        starts(1)=level_offsets(ilevel)
        counts(1)=idx-idx_old
        ! Add current level to the file dataspace selection
        call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

        ! Write level(s) if outputs for this level is true
        if(outputs(ilevel))then
          counts(1)=sz
          if (sz>0) then
            call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,start_ds,counts,status)
          else
            call H5Sselect_none_f(ds_id,status)
            call H5Sselect_none_f(fs_id,status)
          endif
          call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,iig,counts,status,ds_id,fs_id,xfer_prp=mode_id)
          if(status<0)print *,'Problem for writing in hdf5_write_multilevel_int_1D (status=',status,')'
          if(ilevel<nlevelmax)then
            idx=1
            sz=0
            call H5Sselect_none_f(fs_id,status)
          endif
        endif
      end do
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_multilevel_int_1D

    subroutine hdf5_write_multilevel_real8_2D(f_id,name,var,ndim)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(in)::var
      integer,intent(in)::ndim

      integer(kind=mpi_offset_kind)::sz
      integer(hsize_t),dimension(2)::dims,starts,counts
      integer(hsize_t),dimension(1)::start_ds
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,idx_old,igrid,ilevel,ibound,ncache,istart,i,idim
      real(kind=8),dimension(buffer_size)::xdp

      ! Create memory dataspace
      dims(1)=buffer_size
      call H5Screate_simple_f(1,dims,ds_id,status)

      ! Create file dataspace
      dims(1)=size_max
      dims(2)=ndim
      call H5Screate_simple_f(2,dims,fs_id,status)

      ! Create dataset
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,fs_id,da_id,status)

      counts(2)=1
      start_ds(1)=0
      ! Loop on the dimensions
      do idim=1,ndim
        idx=1
        sz=0
        start_ds(1)=0
        starts(2)=idim-1
        call H5Sselect_none_f(fs_id,status)
        ! Loop on the refinement levels
        do ilevel=1,nlevelmax
          idx_old=idx
          do ibound=1,nboundary+ncpu
            if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
            else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
            end if
            if(ncache>0)then
              igrid=istart
              do i=idx,idx-1+ncache
                xdp(i)=var(igrid,idim)
                igrid=next(igrid)
              end do
              sz=sz+ncache
              idx=idx+ncache
            end if
          end do
          starts(1)=level_offsets(ilevel)
          counts(1)=idx-idx_old
          ! Add current level to the file dataspace selection
          call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

          ! Write level(s) if outputs for this level is true
          if(outputs(ilevel))then
            counts(1)=sz
            if (sz>0) then
              call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,start_ds,counts,status)
            else
              call H5Sselect_none_f(ds_id,status)
              call H5Sselect_none_f(fs_id,status)
            endif
            call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,xdp,counts,status,ds_id,fs_id,xfer_prp=mode_id)
            if(status<0)print *,'Problem for writing in hdf5_write_multilevel_real8_2D (status=',status,')'
            if(ilevel<nlevelmax)then
              idx=1
              sz=0
              call H5Sselect_none_f(fs_id,status)
            endif
          endif
        end do
      end do
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_multilevel_real8_2D

    subroutine hdf5_write_multilevel_int_2D(f_id,name,var,ndim)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:,:),intent(in)::var
      integer,intent(in)::ndim

      integer(kind=mpi_offset_kind)::sz
      integer(hsize_t),dimension(2)::dims,starts,counts
      integer(hsize_t),dimension(1)::start_ds
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,idx_old,igrid,ilevel,ibound,ncache,istart,i,idim
      integer,dimension(buffer_size)::iig

      ! Create memory dataspace
      dims(1)=buffer_size
      call H5Screate_simple_f(1,dims,ds_id,status)

      ! Create file dataspace
      dims(1)=size_max
      dims(2)=ndim
      call H5Screate_simple_f(2,dims,fs_id,status)

      ! Create dataset
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)

      counts(2)=1
      start_ds(1)=0
      ! Loop on the dimensions
      do idim=1,ndim
        idx=1
        sz=0
        start_ds(1)=0
        starts(2)=idim-1
        call H5Sselect_none_f(fs_id,status)
        ! Loop on the refinement levels
        do ilevel=1,nlevelmax
            idx_old=idx
            do ibound=1,nboundary+ncpu
               if(ibound<=ncpu)then
                  ncache=numbl(ibound,ilevel)
                  istart=headl(ibound,ilevel)
               else
                  ncache=numbb(ibound-ncpu,ilevel)
                  istart=headb(ibound-ncpu,ilevel)
               end if
               if(ncache>0)then
                  igrid=istart
                  do i=idx,idx-1+ncache
                     iig(i)=var(igrid,idim)
                     igrid=next(igrid)
                  end do
                  sz=sz+ncache
                  idx=idx+ncache
              end if
            end do
          starts(1)=level_offsets(ilevel)
          counts(1)=idx-idx_old
          ! Add current level to the file dataspace selection
          call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

          ! Write level(s) if outputs for this level is true
          if(outputs(ilevel))then
            counts(1)=sz
            if (sz>0) then
              call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,start_ds,counts,status)
            else
              call H5Sselect_none_f(ds_id,status)
              call H5Sselect_none_f(fs_id,status)
            endif
            call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,iig,counts,status,ds_id,fs_id,xfer_prp=mode_id)
            if(status<0)print *,'Problem for writing in hdf5_write_multilevel_int_2D (status=',status,')'
            if(ilevel<nlevelmax)then
              idx=1
              sz=0
              call H5Sselect_none_f(fs_id,status)
            endif
          endif
        end do
      end do
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_multilevel_int_2D

    subroutine hdf5_write_multilevel_ncoarse_int_1D(f_id,name,var,ndims)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl,myid
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(in)::var
      integer,intent(in)::ndims

      integer(kind=mpi_offset_kind)::sz
      integer(hsize_t),dimension(2)::dims,starts,counts
      integer(hsize_t),dimension(1)::start_ds
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,idx_old,igrid,ilevel,ibound,ncache,istart,i,ind,iskip
      integer,dimension(buffer_size)::iig


      ! Create memory dataspace
      dims(1)=buffer_size
      call H5Screate_simple_f(1,dims,ds_id,status)

      ! Create file dataspace
      dims(1)=size_max
      dims(2)=ndims
      call H5Screate_simple_f(2,dims,fs_id,status)

      ! Create dataset
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)

      counts(2)=1
      start_ds(1)=0
      ! Loop on the var dimensions
      do ind=1,ndims
        iskip=ncoarse+(ind-1)*ngridmax
        idx=1
        sz=0
        start_ds(1)=0
        starts(2)=ind-1
        call H5Sselect_none_f(fs_id,status)
        ! Loop on the refinement levels
        do ilevel=1,nlevelmax
          idx_old=idx
          do ibound=1,nboundary+ncpu
            if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
              else
                ncache=numbb(ibound-ncpu,ilevel)
                istart=headb(ibound-ncpu,ilevel)
              end if
              if(ncache>0)then
                igrid=istart
                do i=idx,idx-1+ncache
                  iig(i)=var(igrid+iskip)
                  igrid=next(igrid)
                end do
                sz=sz+ncache
                idx=idx+ncache
              end if
            end do
            starts(1)=level_offsets(ilevel)
            counts(1)=idx-idx_old
            ! Add current level to the file dataspace selection
            call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

            ! Write level(s) if outputs for this level is true
            if(outputs(ilevel))then
              counts(1)=sz
              if (sz>0) then
                call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,start_ds,counts,status)
              else
                call H5Sselect_none_f(ds_id,status)
                call H5Sselect_none_f(fs_id,status)
              endif
              call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,iig,counts,status,ds_id,fs_id,xfer_prp=mode_id)
              if(status<0)print *,'Problem for writing in hdf5_write_multilevel_ncoarse_int_1D (status=',status,')'
              if(ilevel<nlevelmax)then
                idx=1
                sz=0
                call H5Sselect_none_f(fs_id,status)
              endif
            endif
          end do
      end do
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)

    end subroutine hdf5_write_multilevel_ncoarse_int_1D


    subroutine hdf5_read_multilevel_int_1D(f_id,name,var)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(out)::var

      integer(hsize_t),dimension(1)::dims,starts,startm,counts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,igrid,ilevel,ilevel2,ilevelstart,ibound,ncache,istart,i,sz
      integer,dimension(buffer_size)::iig

      dims(1)=buffer_size
      startm(1)=0
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)

      call H5Sselect_none_f(fs_id,status)
      ilevelstart=1

      sz=0
      do ilevel=1,nlevelmax
        sz=sz+size_loc_by_level(ilevel)
        starts(1)=level_offsets(ilevel)
        counts(1)=size_loc_by_level(ilevel)
        call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

        if(outputs(ilevel))then
          counts(1)=sz
          if (sz>0) then
            call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,startm,counts,status)
          else
            call H5Sselect_none_f(ds_id,status)
            call H5Sselect_none_f(fs_id,status)
          endif
          call H5Dread_f(da_id,H5T_NATIVE_INTEGER,iig,counts,status,&
                         file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
          if(status<0)print *,'Problem for reading in hdf5_read_multilevel_int_1D (status=',status,')'

          idx=1
          do ilevel2=ilevelstart,ilevel
            do ibound=1,nboundary+ncpu
              if(ibound<=ncpu)then
                ncache=numbl(ibound,ilevel2)
                istart=headl(ibound,ilevel2)
              else
                ncache=numbb(ibound-ncpu,ilevel2)
                istart=headb(ibound-ncpu,ilevel2)
              end if
              if(ncache>0)then
                igrid=istart
                do i=idx,idx-1+ncache
                  var(igrid)=iig(i)
                  igrid=next(igrid)
                end do
                idx=idx+ncache
              end if
            end do
          end do

          if(ilevel<nlevelmax)then
            sz=0
            ilevelstart=ilevel+1
            call H5Sselect_none_f(fs_id,status)
          endif
        endif
      end do

      call H5Sclose_f(fs_id,status)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_read_multilevel_int_1D

    subroutine hdf5_read_multilevel_int_2D(f_id,name,var,ndims,ngridmax2)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:,:),intent(out)::var
      integer,intent(in)::ndims
      integer,intent(in)::ngridmax2

      integer(hsize_t),dimension(1)::dims,startm
      integer(hsize_t),dimension(2)::starts,counts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,igrid,ilevel,ilevel2,ilevelstart,ind,ibound,ncache,istart,i,sz
      integer,dimension(buffer_size)::iig
      integer::tmp1,tmp2


      counts(2)=1
      startm(1)=0
      dims(1)=buffer_size
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)


      do ind=1,ndims
        sz=0
        starts(2)=ind-1
        ilevelstart=1
        call H5Sselect_none_f(fs_id,status)
        do ilevel=1,nlevelmax
          sz=sz+size_loc_by_level(ilevel)
          starts(1)=level_offsets(ilevel)
          counts(1)=size_loc_by_level(ilevel)
          call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

          if(outputs(ilevel))then
            counts(1)=sz
            if (sz>0) then
              call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,startm,counts,status)
            else
              call H5Sselect_none_f(ds_id,status)
              call H5Sselect_none_f(fs_id,status)
            endif
            call H5Dread_f(da_id,H5T_NATIVE_INTEGER,iig,counts,status,&
                           file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
            if(status<0)print *,'Problem for reading in hdf5_read_multilevel_int_2D (status=',status,')'

            idx=1
            do ilevel2=ilevelstart,ilevel
              do ibound=1,nboundary+ncpu
                if(ibound<=ncpu)then
                  ncache=numbl(ibound,ilevel2)
                  istart=headl(ibound,ilevel2)
                else
                  ncache=numbb(ibound-ncpu,ilevel2)
                  istart=headb(ibound-ncpu,ilevel2)
                end if
                if(ncache>0)then
                  if(ngridmax.ne.ngridmax2.and.ilevel>1)then
                    do i=idx,idx-1+ncache
                      tmp1=(iig(i)-ncoarse-1)/ngridmax2
                      tmp2=iig(i)-ncoarse-tmp1*ngridmax2
                      iig(i)=ncoarse+tmp1*ngridmax+tmp2
                    end do
                  end if
                  igrid=istart
                  do i=idx,idx-1+ncache
                    var(igrid,ind)=iig(i)
                    igrid=next(igrid)
                  end do
                  idx=idx+ncache
                end if
              end do
            end do

            if(ilevel<nlevelmax)then
              sz=0
              ilevelstart=ilevel+1
              call H5Sselect_none_f(fs_id,status)
            endif
          endif
        end do
      end do

      call H5Sclose_f(fs_id,status)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_read_multilevel_int_2D

    subroutine hdf5_read_multilevel_real8_2D(f_id,name,var,ndims)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(out)::var
      integer,intent(in)::ndims

      integer(hsize_t),dimension(1)::dims,startm
      integer(hsize_t),dimension(2)::counts,starts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,ind,igrid,ilevel,ilevel2,ilevelstart,ibound,ncache,istart,i,sz
      real(kind=8),dimension(buffer_size)::xxg

      counts(2)=1
      startm(1)=0
      dims(1)=buffer_size
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)

      do ind=1,ndims
        sz=0
        starts(2)=ind-1
        ilevelstart=1
        call H5Sselect_none_f(fs_id,status)
        do ilevel=1,nlevelmax
          sz=sz+size_loc_by_level(ilevel)
          starts(1)=level_offsets(ilevel)
          counts(1)=size_loc_by_level(ilevel)
          call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

          if(outputs(ilevel))then
            counts(1)=sz
            if (sz>0) then
              call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,startm,counts,status)
            else
              call H5Sselect_none_f(ds_id,status)
              call H5Sselect_none_f(fs_id,status)
            endif
            call H5Dread_f(da_id,H5T_NATIVE_DOUBLE,xxg,counts,status,&
                           file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
            if(status<0)print *,'Problem for reading in hdf5_read_multilevel_real8_2D (status=',status,')'

            idx=1
            do ilevel2=ilevelstart,ilevel
              do ibound=1,nboundary+ncpu
                if(ibound<=ncpu)then
                  ncache=numbl(ibound,ilevel2)
                  istart=headl(ibound,ilevel2)
                else
                  ncache=numbb(ibound-ncpu,ilevel2)
                  istart=headb(ibound-ncpu,ilevel2)
                end if
                if(ncache>0)then
                  igrid=istart
                  do i=idx,idx-1+ncache
                    var(igrid,ind)=xxg(i)
                    igrid=next(igrid)
                  end do
                  idx=idx+ncache
                end if
              end do
            end do

            if(ilevel<nlevelmax)then
              sz=0
              ilevelstart=ilevel+1
              call H5Sselect_none_f(fs_id,status)
            endif
          endif
        end do
      end do

      call H5Sclose_f(fs_id,status)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_read_multilevel_real8_2D

    subroutine hdf5_read_multilevel_ncoarse_int_1D(f_id,name,var,ndims)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(out)::var
      integer,intent(in)::ndims

      integer(hsize_t),dimension(1)::dims,startm
      integer(hsize_t),dimension(2)::starts,counts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idx,igrid,ilevel,ilevel2,ilevelstart,ibound,ind,ncache,istart,i,iskip,sz
      integer,dimension(buffer_size)::iig

      counts(2)=1
      dims(1)=buffer_size
      startm(1)=0
      call H5Screate_simple_f(1,dims,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      do ind=1,ndims
        starts(2)=ind-1
        iskip=ncoarse+(ind-1)*ngridmax
        sz=0
        ilevelstart=1
        call H5Sselect_none_f(fs_id,status)
        do ilevel=1,nlevelmax
          sz=sz+size_loc_by_level(ilevel)
          starts(1)=level_offsets(ilevel)
          counts(1)=size_loc_by_level(ilevel)
          call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_OR_F,starts,counts,status)

          if(outputs(ilevel))then
            counts(1)=sz
            if (sz>0) then
              call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,startm,counts,status)
            else
              call H5Sselect_none_f(ds_id,status)
              call H5Sselect_none_f(fs_id,status)
            endif
            call H5Dread_f(da_id,H5T_NATIVE_INTEGER,iig,counts,status,&
                           file_space_id=fs_id,mem_space_id=ds_id,xfer_prp=mode_id)
            if(status<0)print *,'Problem for reading in hdf5_read_multilevel_ncoarse_int_1D (status=',status,')'

            idx=1
            do ilevel2=ilevelstart,ilevel
              do ibound=1,nboundary+ncpu
                if(ibound<=ncpu)then
                  ncache=numbl(ibound,ilevel2)
                  istart=headl(ibound,ilevel2)
                else
                  ncache=numbb(ibound-ncpu,ilevel2)
                  istart=headb(ibound-ncpu,ilevel2)
                end if
                if(ncache>0)then
                  igrid=istart
                  do i=idx,idx-1+ncache
                    var(igrid+iskip)=iig(i)
                    igrid=next(igrid)
                  end do
                  idx=idx+ncache
                end if
              end do
            end do

            if(ilevel<nlevelmax)then
              sz=0
              ilevelstart=ilevel+1
              call H5Sselect_none_f(fs_id,status)
            endif
          endif
        end do
      end do

      call H5Sclose_f(fs_id,status)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_read_multilevel_ncoarse_int_1D

    subroutine hdf5_write_part_int_1D(f_id,name,var,npart,npart_tot,offset)
      use amr_parameters,only:ndim
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(in)::var
      integer,intent(in)::npart
      integer(kind=hsize_t),intent(in)::npart_tot
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(1)::dims,starts,counts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::ipart,i
      integer,dimension(npart)::ii

      dims(1)=npart
      call H5Screate_simple_f(1,dims,ds_id,status)
      dims(1)=npart_tot
      call H5Screate_simple_f(1,dims,fs_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)
      counts(1)=npart
      starts(1)=0
      call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,starts,counts,status)
      starts(1)=offset
      ipart=0
      do i=1,npartmax
        if(levelp(i)>0)then
          ipart=ipart+1
          ii(ipart)=var(i)
        end if
      end do
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,ii,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_part_int_1D

    subroutine hdf5_write_part_real8_1D(f_id,name,var,npart,npart_tot,offset)
      use amr_parameters,only:ndim
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(in)::var
      integer,intent(in)::npart
      integer(kind=hsize_t),intent(in)::npart_tot
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(1)::dims,starts,counts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::ipart,i
      real(kind=8),dimension(npart)::xdp

      dims(1)=npart
      call H5Screate_simple_f(1,dims,ds_id,status)
      dims(1)=npart_tot
      call H5Screate_simple_f(1,dims,fs_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,fs_id,da_id,status)
      counts(1)=npart
      starts(1)=0
      call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,starts,counts,status)
      starts(1)=offset
      ipart=0
      do i=1,npartmax
        if(levelp(i)>0)then
          ipart=ipart+1
          xdp(ipart)=var(i)
        end if
      end do
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,xdp,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_part_real8_1D

    subroutine hdf5_write_part_real8_2D(f_id,name,var,npart,npart_tot,offset)
      use amr_parameters,only:ndim
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(in)::var
      integer,intent(in)::npart
      integer(kind=hsize_t),intent(in)::npart_tot
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(2)::dims,starts,counts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::ipart,i,idim
      real(kind=8),dimension(npart)::xdp

      dims(1)=npart
      call H5Screate_simple_f(1,dims,ds_id,status)
      dims(1)=npart_tot
      dims(2)=ndim
      call H5Screate_simple_f(2,dims,fs_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,fs_id,da_id,status)
      counts(1)=npart
      counts(2)=1
      starts(1)=0
      call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,starts,counts,status)
      starts(1)=offset
      do idim=1,ndim
         ipart=0
         do i=1,npartmax
            if(levelp(i)>0)then
               ipart=ipart+1
               xdp(ipart)=var(i,idim)
            end if
         end do
         starts(2)=idim-1
         call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
         call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,xdp,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      end do
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_part_real8_2D

    subroutine hdf5_read_part_int_1D(f_id,name,var,npart,offset)
      use amr_parameters,only:ndim
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(out)::var
      integer,intent(in)::npart
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(1)::starts,counts
      integer(hid_t)::da_id,ds_id,fs_id

      counts(1)=npart
      call H5Screate_simple_f(1,counts,ds_id,status)

      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=offset
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_INTEGER,var,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_read_part_int_1D

    subroutine hdf5_read_part_real8_1D(f_id,name,var,npart,offset)
      use amr_parameters,only:ndim
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(out)::var
      integer,intent(in)::npart
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(1)::starts,counts
      integer(hid_t)::da_id,ds_id,fs_id

      counts(1)=npart
      call H5Screate_simple_f(1,counts,ds_id,status)

      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      starts(1)=offset
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dread_f(da_id,H5T_NATIVE_DOUBLE,var,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_read_part_real8_1D

    subroutine hdf5_read_part_real8_2D(f_id,name,var,npart,offset)
      use amr_parameters,only:ndim
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(out)::var
      integer,intent(in)::npart
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(2)::starts,counts
      integer(hid_t)::da_id,ds_id,fs_id
      integer::idim

      counts(1)=npartmax
      counts(2)=ndim
      call H5Screate_simple_f(2,counts,ds_id,status)
      call H5Dopen_f(f_id,name,da_id,status)
      call H5Dget_space_f(da_id,fs_id,status)
      counts(1)=npart
      counts(2)=1
      do idim=1,ndim
        starts(1)=0
        starts(2)=idim-1
         call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,starts,counts,status)
         starts(1)=offset
         call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
         call H5Dread_f(da_id,H5T_NATIVE_DOUBLE,var,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      end do
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_read_part_real8_2D

    subroutine hdf5_write_irregular_int_1D(f_id,name,var,nelt,nelt_tot,offset)
      use amr_parameters,only:ndim
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(in)::var
      integer,intent(in)::nelt
      integer(kind=hsize_t),intent(in)::nelt_tot
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(1)::dims,starts,counts
      integer(hid_t)::da_id,ds_id,fs_id

      dims(1)=nelt
      call H5Screate_simple_f(1,dims,ds_id,status)
      dims(1)=nelt_tot
      call H5Screate_simple_f(1,dims,fs_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_INTEGER,fs_id,da_id,status)
      counts(1)=nelt
      starts(1)=0
      call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,starts,counts,status)
      starts(1)=offset
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_INTEGER,var,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_irregular_int_1D

    subroutine hdf5_write_irregular_real8_1D(f_id,name,var,nelt,nelt_tot,offset)
      use amr_parameters,only:ndim
      use pm_parameters,only:npartmax

      integer(hid_t),intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(in)::var
      integer,intent(in)::nelt
      integer(kind=hsize_t),intent(in)::nelt_tot
      integer(kind=hsize_t),intent(in)::offset

      integer(hsize_t),dimension(1)::dims,starts,counts
      integer(hid_t)::da_id,ds_id,fs_id

      dims(1)=nelt
      call H5Screate_simple_f(1,dims,ds_id,status)
      dims(1)=nelt_tot
      call H5Screate_simple_f(1,dims,fs_id,status)
      call H5Dcreate_f(f_id,name,H5T_NATIVE_DOUBLE,fs_id,da_id,status)
      counts(1)=nelt
      starts(1)=0
      call H5Sselect_hyperslab_f(ds_id,H5S_SELECT_SET_F,starts,counts,status)
      starts(1)=offset
      call H5Sselect_hyperslab_f(fs_id,H5S_SELECT_SET_F,starts,counts,status)
      call H5Dwrite_f(da_id,H5T_NATIVE_DOUBLE,var,counts,status,ds_id,fs_id,xfer_prp=mode_id)
      call H5Dclose_f(da_id,status)
      call H5Sclose_f(fs_id,status)
      call H5Sclose_f(ds_id,status)
    end subroutine hdf5_write_irregular_real8_1D
end module ramses_hdf5
#endif
