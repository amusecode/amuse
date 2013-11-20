! Copyright (c) 2007-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef PNETCDF
module ramses_pnetcdf
  use amr_commons,only:MPI_COMM_COMP,myid,pario_version,ramses_version
  use mpi
  use pnetcdf
  implicit none

  integer::status
  integer(kind=mpi_offset_kind),parameter::one=1
  integer(kind=mpi_offset_kind),dimension(1),parameter::one1=1

  interface pnetcdf_get_var
     module procedure pnetcdf_get_var_int4,     pnetcdf_get_var_real8,&
                      pnetcdf_get_var_real8_1D, pnetcdf_get_var_int4_2D,&
                      pnetcdf_get_var_int8_2D,  pnetcdf_get_var_char
  end interface pnetcdf_get_var

  interface pnetcdf_get_vara
     module procedure pnetcdf_get_vara_int4,     pnetcdf_get_vara_real8, &
                      pnetcdf_get_vara_int4_1D,  pnetcdf_get_vara_int4_2D, &
                      pnetcdf_get_vara_real8_1D, pnetcdf_get_vara_real8_2D
  end interface pnetcdf_get_vara

  interface pnetcdf_put_var
     module procedure pnetcdf_put_var_int4,     pnetcdf_put_var_real4, &
                      pnetcdf_put_var_real8,    pnetcdf_put_var_real8_1D, &
                      pnetcdf_put_var_int4_2D,  pnetcdf_put_var_int8_2D, &
                      pnetcdf_put_var_char
  end interface pnetcdf_put_var

  interface pnetcdf_put_distrib_var
     module procedure pnetcdf_put_distrib_var_int4,    pnetcdf_put_distrib_var_real4, &
                      pnetcdf_put_distrib_var_real8,   pnetcdf_put_distrib_var_int4_1D, &
                      pnetcdf_put_distrib_var_int4_2D, pnetcdf_put_distrib_var_real8_1D, &
                      pnetcdf_put_distrib_var_real8_2D
  end interface pnetcdf_put_distrib_var

  contains
    subroutine create_pnetcdf_file(fileloc,f_id)
      use amr_parameters,only:MAXLINE
      use hints,only:hints_set

      integer,intent(out)::f_id
      character(LEN=*),intent(in)::fileloc

      integer::infos,ierr,mode
      character(len=80)::pnetcdf_version
      character(len=MAXLINE)::hostname
      character(len=8)::datec
      character(len=10)::timec


      ! mode = NF_CLOBBER+NF_64BIT_OFFSET 
      mode = NF_CLOBBER+NF_64BIT_DATA !NF_64BIT_DATA is only supported from PnetCDF 1.1.0
      call MPI_Info_create(infos,ierr)
      call hints_set(infos)
      status = nfmpi_create(MPI_COMM_COMP,fileloc,mode,infos,f_id)
      call CHKERR(status)
      call MPI_Info_free(infos,ierr)

      ! Add version numbers as attributes
      status = nfmpi_put_att_text(f_id,NF_GLOBAL,"RAMSES_version",&
               int(len(ramses_version),kind=mpi_offset_kind),ramses_version);call CHKERR(status);
      status = nfmpi_put_att_text(f_id,NF_GLOBAL,"Parallel_IO_version",&
               int(len(pario_version),kind=mpi_offset_kind),pario_version);call CHKERR(status);
      pnetcdf_version=nfmpi_inq_libvers()
      status = nfmpi_put_att_text(f_id,NF_GLOBAL,"Parallel-NetCDF_version",&
               int(len_trim(pnetcdf_version),kind=mpi_offset_kind),trim(pnetcdf_version));
               call CHKERR(status);

      ! Add hostname as attribute
      hostname=' ' ! Necessary to have a valid Fortran string
      if(myid==1) call gethname(hostname)
      call MPI_Bcast(hostname,len(hostname),MPI_CHARACTER,0,MPI_COMM_COMP,ierr)
      status = nfmpi_put_att_text(f_id,NF_GLOBAL,"hostname",&
               int(len_trim(hostname),kind=mpi_offset_kind),trim(hostname));call CHKERR(status);

      ! Add date and time as attributes
      if (myid==1) call date_and_time(date=datec,time=timec)
      call MPI_Bcast(datec,len(datec),MPI_CHARACTER,0,MPI_COMM_COMP,ierr)
      call MPI_Bcast(timec,len(timec),MPI_CHARACTER,0,MPI_COMM_COMP,ierr)
      status = nfmpi_put_att_text(f_id,NF_GLOBAL,"date",&
               int(len(datec),kind=mpi_offset_kind),datec);call CHKERR(status);
      status = nfmpi_put_att_text(f_id,NF_GLOBAL,"time",&
               int(len(timec),kind=mpi_offset_kind),timec);call CHKERR(status);

    end subroutine create_pnetcdf_file

    subroutine open_pnetcdf_file(fileloc,f_id)
      use amr_parameters,only:maxline
      use hints,only:hints_set

      integer,intent(out)::f_id
      character(LEN=*),intent(in)::fileloc

      integer::infos,ierr,mode

      mode = NF_NOWRITE
      call MPI_Info_create(infos,ierr)
      call hints_set(infos)
      status = nfmpi_open(MPI_COMM_COMP,fileloc,mode,infos,f_id)
      call CHKERR(status)
      call MPI_Info_free(infos,ierr)
    end subroutine open_pnetcdf_file

    subroutine pnetcdf_get_dim(f_id,name,dim)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      integer,intent(out)::dim

      integer::dim_id
      integer(kind=mpi_offset_kind)::dim8

      status = nfmpi_inq_dimid(f_id,name,dim_id)
      status = nfmpi_inq_dimlen(f_id,dim_id,dim8)

      dim = dim8 !Type conversion
    end subroutine pnetcdf_get_dim



    subroutine pnetcdf_get_var_real8(f_id,name,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      real(kind=8),intent(out)::var

      integer::var_id
      real(kind=8),dimension(1)::var1

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_var_double_all(f_id,var_id,var1)
      var = var1(1)
    end subroutine pnetcdf_get_var_real8

    subroutine pnetcdf_get_var_int4(f_id,name,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      integer(kind=4),intent(out)::var

      integer::var_id
      integer(kind=4),dimension(1)::var1

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_var_int_all(f_id,var_id,var1)
      var = var1(1)
    end subroutine pnetcdf_get_var_int4

    subroutine pnetcdf_get_var_real8_1D(f_id,name,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      real(kind=8),intent(out),dimension(:)::var

      integer::var_id

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_var_double_all(f_id,var_id,var)
    end subroutine pnetcdf_get_var_real8_1D

    subroutine pnetcdf_get_var_int4_2D(f_id,name,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      integer(kind=4),dimension(:,:),intent(out)::var

      integer::var_id

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_var_int_all(f_id,var_id,var)
    end subroutine pnetcdf_get_var_int4_2D

    subroutine pnetcdf_get_var_int8_2D(f_id,name,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      integer(kind=8),dimension(:,:),intent(out)::var

      integer::var_id
      real(kind=8),dimension(:,:),allocatable::var_dbl

      allocate(var_dbl(size(var,1),size(var,2)))

      status = nfmpi_inq_varid(f_id,name,var_id)
      !int*8 not supported by NetCDF => written as real*8 (bit by bit copy)
      status = nfmpi_get_var_double_all(f_id,var_id,var_dbl)

      var = reshape(transfer(var_dbl,1_8,size(var)),shape(var))

      deallocate(var_dbl)
    end subroutine pnetcdf_get_var_int8_2D

    subroutine pnetcdf_get_var_char(f_id,name,var)
      !TODO:secure the length of the string
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      character(len=*),intent(out)::var

      integer::var_id

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_var_text_all(f_id,var_id,var)
    end subroutine pnetcdf_get_var_char


    subroutine pnetcdf_get_vara_int4(f_id,name,offsets,sizes,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      integer(kind=4),intent(out)::var
      integer(kind=mpi_offset_kind),dimension(:),intent(in)::offsets,sizes

      integer::var_id
      integer(kind=4),dimension(1)::var1

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_vara_int_all(f_id,var_id,offsets,sizes,var1)
      var = var1(1)
    end subroutine pnetcdf_get_vara_int4

    subroutine pnetcdf_get_vara_real8(f_id,name,offsets,sizes,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      real(kind=8),intent(out)::var
      integer(kind=mpi_offset_kind),dimension(:),intent(in)::offsets,sizes

      integer::var_id
      real(kind=8),dimension(1)::var1

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_vara_double_all(f_id,var_id,offsets,sizes,var1)
      var = var1(1)
    end subroutine pnetcdf_get_vara_real8

    subroutine pnetcdf_get_vara_int4_1D(f_id,name,offsets,sizes,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      integer(kind=4),intent(out),dimension(:)::var
      integer(kind=mpi_offset_kind),dimension(:),intent(in)::offsets,sizes

      integer::var_id

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_vara_int_all(f_id,var_id,offsets,sizes,var)
    end subroutine pnetcdf_get_vara_int4_1D

    subroutine pnetcdf_get_vara_int4_2D(f_id,name,offsets,sizes,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      integer(kind=4),intent(out),dimension(:,:)::var
      integer(kind=mpi_offset_kind),dimension(:),intent(in)::offsets,sizes

      integer::var_id

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_vara_int_all(f_id,var_id,offsets,sizes,var)
    end subroutine pnetcdf_get_vara_int4_2D

    subroutine pnetcdf_get_vara_real8_1D(f_id,name,offsets,sizes,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      real(kind=8),intent(out),dimension(:)::var
      integer(kind=mpi_offset_kind),dimension(:),intent(in)::offsets,sizes

      integer::var_id

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_vara_double_all(f_id,var_id,offsets,sizes,var)
    end subroutine pnetcdf_get_vara_real8_1D

    subroutine pnetcdf_get_vara_real8_2D(f_id,name,offsets,sizes,var)
      integer,intent(in)::f_id
      character(len=*),intent(in)::name
      real(kind=8),intent(out),dimension(:,:)::var
      integer(kind=mpi_offset_kind),dimension(:),intent(in)::offsets,sizes

      integer::var_id

      status = nfmpi_inq_varid(f_id,name,var_id)
      status = nfmpi_get_vara_double_all(f_id,var_id,offsets,sizes,var)
    end subroutine pnetcdf_get_vara_real8_2D


    subroutine pnetcdf_put_var_int4(f_id,var_id,var)
      integer,intent(in)::f_id,var_id
      integer(kind=4),intent(in)::var

      status = nfmpi_put_var1_int(f_id,var_id,one1,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_var_int4

    subroutine pnetcdf_put_var_real4(f_id,var_id,var)
      integer,intent(in)::f_id,var_id
      real(kind=4),intent(in)::var

      status = nfmpi_put_var1_real(f_id,var_id,one1,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_var_real4

    subroutine pnetcdf_put_var_real8(f_id,var_id,var)
      integer,intent(in)::f_id,var_id
      real(kind=8),intent(in)::var

      status = nfmpi_put_var1_double(f_id,var_id,one1,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_var_real8

    subroutine pnetcdf_put_var_real8_1D(f_id,var_id,var,count)
      integer,intent(in)::f_id,var_id
      real(kind=8),dimension(:),intent(in)::var
      integer,optional,intent(in)::count

      integer(kind=mpi_offset_kind),dimension(1)::numelt

      if(present(count))then
         numelt(1) = count
      else
         numelt(1) = size(var,1)
      end if

      status = nfmpi_put_vara_double(f_id,var_id,one1,numelt,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_var_real8_1D

    subroutine pnetcdf_put_var_int4_2D(f_id,var_id,var,counts)
      integer,intent(in)::f_id,var_id
      integer(kind=4),dimension(:,:),intent(in)::var
      integer,dimension(2),optional,intent(in)::counts

      integer(kind=mpi_offset_kind),dimension(2)::numelt,offsets

      if(present(counts))then
         numelt(:) = counts(:)
      else
         numelt(:) = shape(var)
      end if
      offsets(:) = one

      status = nfmpi_put_vara_int(f_id,var_id,offsets,numelt,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_var_int4_2D

    subroutine pnetcdf_put_var_int8_2D(f_id,var_id,var,counts)
      integer,intent(in)::f_id,var_id
      integer(kind=8),dimension(:,:),intent(in)::var
      integer,dimension(2),optional,intent(in)::counts

      integer(kind=mpi_offset_kind),dimension(2)::numelt,offsets
      integer::tot_elt

      if(present(counts))then
         numelt(:) = counts(:)
      else
         numelt(:) = shape(var)
      end if
      offsets(:) = one
      tot_elt=numelt(1)*numelt(2)

      !int*8 not supported by NetCDF => written as real*8 (bit by bit copy)
      status = nfmpi_put_vara_double(f_id,var_id,offsets,numelt, &
                                reshape(transfer(var(:,:),1.0_8,tot_elt),numelt(:)))
      call CHKERR(status)
    end subroutine pnetcdf_put_var_int8_2D

    subroutine pnetcdf_put_var_char(f_id,var_id,var)
      integer,intent(in)::f_id,var_id
      character(len=*),intent(in)::var

      status = nfmpi_put_var_text(f_id,var_id,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_var_char



    subroutine pnetcdf_put_distrib_var_int4(f_id,var_id,var)
      use amr_commons,only:myid

      integer,intent(in)::f_id,var_id
      integer(kind=4),intent(in)::var

      integer(kind=mpi_offset_kind),dimension(1)::offset
      integer(kind=4),dimension(1)::var1

      offset(1) = myid
      var1(1) = var
      status = nfmpi_put_vara_int_all(f_id,var_id,offset,one1,var1)
      call CHKERR(status)
    end subroutine pnetcdf_put_distrib_var_int4

    subroutine pnetcdf_put_distrib_var_real4(f_id,var_id,var)
      use amr_commons,only:myid

      integer,intent(in)::f_id,var_id
      real(kind=4),intent(in)::var

      integer(kind=mpi_offset_kind),dimension(1)::offset
      real(kind=4),dimension(1)::var1

      offset(1) = myid
      var1(1) = var
      status = nfmpi_put_vara_real_all(f_id,var_id,offset,one1,var1)
      call CHKERR(status)
    end subroutine pnetcdf_put_distrib_var_real4

    subroutine pnetcdf_put_distrib_var_real8(f_id,var_id,var)
      use amr_commons,only:myid

      integer,intent(in)::f_id,var_id
      real(kind=8),intent(in)::var

      integer(kind=mpi_offset_kind),dimension(1)::offset
      real(kind=8),dimension(1)::var1

      offset(1) = myid
      var1(1) = var
      status = nfmpi_put_vara_double_all(f_id,var_id,offset,one1,var1)
      call CHKERR(status)
    end subroutine pnetcdf_put_distrib_var_real8

    subroutine pnetcdf_put_distrib_var_int4_1D(f_id,var_id,var)
      use amr_commons,only:myid

      integer,intent(in)::f_id,var_id
      integer(kind=4),dimension(:),intent(in)::var

      integer(kind=mpi_offset_kind),dimension(2)::offsets,sizes

      offsets(1) = 1
      offsets(2) = myid
      sizes(1) = size(var,1)
      sizes(2) = 1

      status = nfmpi_put_vara_int_all(f_id,var_id,offsets,sizes,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_distrib_var_int4_1D

    subroutine pnetcdf_put_distrib_var_int4_2D(f_id,var_id,var)
      use amr_commons,only:myid

      integer,intent(in)::f_id,var_id
      integer(kind=4),dimension(:,:),intent(in)::var

      integer(kind=mpi_offset_kind),dimension(3)::offsets,sizes

      offsets(1) = 1
      offsets(2) = 1
      offsets(3) = myid
      sizes(1) = size(var,1)
      sizes(2) = size(var,2)
      sizes(3) = 1

      status = nfmpi_put_vara_int_all(f_id,var_id,offsets,sizes,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_distrib_var_int4_2D

    subroutine pnetcdf_put_distrib_var_real8_1D(f_id,var_id,var)
      use amr_commons,only:myid

      integer,intent(in)::f_id,var_id
      real(kind=8),dimension(:),intent(in)::var

      integer(kind=mpi_offset_kind),dimension(2)::offsets,sizes

      offsets(1) = 1
      offsets(2) = myid
      sizes(1) = size(var,1)
      sizes(2) = 1

      status = nfmpi_put_vara_double_all(f_id,var_id,offsets,sizes,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_distrib_var_real8_1D

    subroutine pnetcdf_put_distrib_var_real8_2D(f_id,var_id,var)
      use amr_commons,only:myid

      integer,intent(in)::f_id,var_id
      real(kind=8),dimension(:,:),intent(in)::var

      integer(kind=mpi_offset_kind),dimension(3)::offsets,sizes

      offsets(1) = 1
      offsets(2) = 1
      offsets(3) = myid
      sizes(1) = size(var,1)
      sizes(2) = size(var,2)
      sizes(3) = 1

      status = nfmpi_put_vara_double_all(f_id,var_id,offsets,sizes,var)
      call CHKERR(status)
    end subroutine pnetcdf_put_distrib_var_real8_2D



    subroutine CHKERR(status)
      integer::status

      if ((status)/=NF_NOERR) then
         write(*,*) nfmpi_strerror(status)
         call clean_stop
      end if
    end subroutine CHKERR


    subroutine pnetcdf_write_multilevel_int_1D(f_id,var,var_id)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      integer,dimension(:),intent(in)::var
      integer,intent(in)::var_id

      integer(kind=mpi_offset_kind),dimension(1)::sz,offset
      integer::idx,igrid,ilevel,ibound,ncache,istart,i,levelstart,nlevel
      integer,dimension(ncache_max)::iig
!      integer,dimension(buffer_size)::iig
!      integer,dimension(nlevelmax)::offsets,sizes

#if 0
!A partir de pnetcdf 1.2.0 !!!
      integer,dimension(nlevelmax)::req,statuses

      nlevel=0
      levelstart=1
      do ilevel=1,nlevelmax
         idx=1
         sz(1)=0
         offset(1)=level_offsets(ilevel)+1
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
               sz(1)=sz(1)+ncache
               idx=idx+ncache
            end if
         end do
         nlevel=nlevel+1
         status = nfmpi_iput_vara_int(f_id,var_id,offset,sz,iig,req(ilevel))
         if(outputs(ilevel)) then
           status = nfmpi_wait_all(f_id,nlevel,req(levelstart:ilevel),statuses(levelstart:ilevel))
           nlevel=0
           levelstart=ilevel+1
         end if
      end do
#endif
#if 0
      offsets(:)=level_offsets(:)
      sizes(:)=size_loc_by_level(:)
      idx=1
      sz(1)=0
      nlevel=0
      levelstart=1
      do ilevel=1,nlevelmax
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
               sz(1)=sz(1)+ncache
               idx=idx+ncache
            end if
         end do

         nlevel=nlevel+1

         if(outputs(ilevel))then
!print *,myid,' : levst=',levelstart,' levend=',ilevel,' sizes=',sizes(levelstart:ilevel),' offsets=',offsets(levelstart:ilevel)
           call MPI_Type_indexed(nlevel,sizes(levelstart:ilevel),offsets(levelstart:ilevel),MPI_INTEGER,type_fs,status)
           call MPI_Type_commit(type_fs,status)
           status = nfmpi_put_vars_all(f_id,var_id,int(0,kind=mpi_offset_kind),int(1,kind=mpi_offset_kind),&
                                       int(1,kind=mpi_offset_kind),iig,type_fs)
           call MPI_Type_free(type_fs,status)
           idx=1
           sz(1)=0
           nlevel=0
           levelstart=ilevel+1
         endif
      end do
#endif
#if 1
      do ilevel=1,nlevelmax
         idx=1
         sz(1)=0
         offset(1)=level_offsets(ilevel)+1
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
               sz(1)=sz(1)+ncache
               idx=idx+ncache
            end if
         end do
         status = nfmpi_put_vara_int_all(f_id,var_id,offset,sz,iig)
      end do
#endif
    end subroutine pnetcdf_write_multilevel_int_1D

    subroutine pnetcdf_write_multilevel_int_2D(f_id,var,var_id,ndims)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      integer,dimension(:,:),intent(in)::var
      integer,intent(in)::var_id
      integer,intent(in)::ndims

      integer(kind=mpi_offset_kind)::sz
      integer::idx,igrid,ilevel,ibound,ncache,istart,i,ind
      integer,dimension(ncache_max)::iig
      integer(kind=mpi_offset_kind),dimension(2)::offsets,sizes

      sizes(2)=1
      do ind=1,ndims
         offsets(2)=ind
         do ilevel=1,nlevelmax
            idx=1
            sz=0
            offsets(1)=level_offsets(ilevel)+1
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
                     iig(i)=var(igrid,ind)
                     igrid=next(igrid)
                  end do
                  sz=sz+ncache
                  idx=idx+ncache
               end if
            end do
            sizes(1)=sz
            status = nfmpi_put_vara_int_all(f_id,var_id,offsets,sizes,iig)
         end do
      end do
    end subroutine pnetcdf_write_multilevel_int_2D

    subroutine pnetcdf_write_multilevel_real8_2D(f_id,var,var_id,ndims)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      real(kind=8),dimension(:,:),intent(in)::var
      integer,intent(in)::var_id
      integer,intent(in)::ndims

      integer(kind=mpi_offset_kind)::sz
      integer::idx,igrid,ilevel,ibound,ncache,istart,i,ind
      real(kind=8),dimension(ncache_max)::xxg
      integer(kind=mpi_offset_kind),dimension(2)::offsets,sizes

      sizes(2)=1
      do ind=1,ndims
         offsets(2)=ind
         do ilevel=1,nlevelmax
            idx=1
            sz=0
            offsets(1)=level_offsets(ilevel)+1
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
                     xxg(i)=var(igrid,ind)
                     igrid=next(igrid)
                  end do
                  sz=sz+ncache
                  idx=idx+ncache
               end if
            end do
            sizes(1)=sz
            status = nfmpi_put_vara_double_all(f_id,var_id,offsets,sizes,xxg)
         end do
      end do
    end subroutine pnetcdf_write_multilevel_real8_2D

    subroutine pnetcdf_write_multilevel_ncoarse_int_1D(f_id,var,var_id,ndims)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl,myid
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      integer,dimension(:),intent(in)::var
      integer,intent(in)::var_id
      integer,intent(in)::ndims

      integer(kind=mpi_offset_kind)::sz
      integer::idx,igrid,ilevel,ibound,ncache,istart,i,ind,iskip
      integer,dimension(ncache_max)::iig
      integer(kind=mpi_offset_kind),dimension(2)::offsets,sizes

      sizes(2)=1
      do ind=1,ndims
         offsets(2)=ind
         iskip=ncoarse+(ind-1)*ngridmax
         do ilevel=1,nlevelmax
            idx=1
            sz=0
            offsets(1)=level_offsets(ilevel)+1
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
            sizes(1)=sz
            status = nfmpi_put_vara_int_all(f_id,var_id,offsets,sizes,iig)
        end do
      end do
    end subroutine pnetcdf_write_multilevel_ncoarse_int_1D



    subroutine pnetcdf_read_multilevel_int_1D(f_id,name,var)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(out)::var

      integer(kind=mpi_offset_kind),dimension(1)::offsets,sizes
      integer::idx,igrid,ilevel,ibound,ncache,istart,i
      integer::var_id
      integer,dimension(ncache_max)::iig

      status = nfmpi_inq_varid(f_id,name,var_id)
      do ilevel=1,nlevelmax
        idx=1
        offsets(1)=level_offsets(ilevel)+1
        sizes(1)=size_loc_by_level(ilevel)
        status = nfmpi_get_vara_int_all(f_id,var_id,offsets,sizes,iig)
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
                  var(igrid)=iig(i)
                  igrid=next(igrid)
               end do
               idx=idx+ncache
            end if
         end do
      end do
    end subroutine pnetcdf_read_multilevel_int_1D

    subroutine pnetcdf_read_multilevel_int_2D(f_id,name,var,ndims,ngridmax2)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:,:),intent(out)::var
      integer,intent(in)::ndims
      integer,intent(in)::ngridmax2

      integer(kind=mpi_offset_kind),dimension(2)::sizes,offsets
      integer::idx,igrid,ilevel,ibound,ncache,istart,i,idim
      integer::var_id
      integer::tmp1,tmp2
      integer,dimension(ncache_max)::iig

      status = nfmpi_inq_varid(f_id,name,var_id)
      sizes(2)=1
      do idim=1,ndims
         offsets(2)=idim
         do ilevel=1,nlevelmax
            idx=1
            offsets(1)=level_offsets(ilevel)+1
            sizes(1)=size_loc_by_level(ilevel)
            status = nfmpi_get_vara_int_all(f_id,var_id,offsets,sizes,iig)
            do ibound=1,nboundary+ncpu
               if(ibound<=ncpu)then
                  ncache=numbl(ibound,ilevel)
                  istart=headl(ibound,ilevel)
               else
                  ncache=numbb(ibound-ncpu,ilevel)
                  istart=headb(ibound-ncpu,ilevel)
               end if
               if(ncache>0)then
                  if(ngridmax>ngridmax2.and.ilevel>1)then
                     do i=idx,idx-1+ncache
                        tmp1=(iig(i)-ncoarse-1)/ngridmax2
                        tmp2=iig(i)-ncoarse-tmp1*ngridmax2
                        iig(i)=ncoarse+tmp1*ngridmax+tmp2
                     end do
                  end if
                  igrid=istart
                  do i=idx,idx-1+ncache
                     var(igrid,idim)=iig(i)
                     igrid=next(igrid)
                  end do
                  idx=idx+ncache
               end if
            end do
         end do
      end do
    end subroutine pnetcdf_read_multilevel_int_2D

    subroutine pnetcdf_read_multilevel_real8_2D(f_id,name,var,ndims)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(out)::var
      integer,intent(in)::ndims

      integer(kind=mpi_offset_kind),dimension(2)::sizes,offsets
      integer::idx,igrid,ilevel,ibound,ncache,istart,i,idim
      integer::var_id
      real(kind=8),dimension(ncache_max)::xxg

      status = nfmpi_inq_varid(f_id,name,var_id)
      sizes(2)=1
      do idim=1,ndims
         offsets(2)=idim
         do ilevel=1,nlevelmax
            idx=1
            offsets(1)=level_offsets(ilevel)+1
            sizes(1)=size_loc_by_level(ilevel)
            status = nfmpi_get_vara_double_all(f_id,var_id,offsets,sizes,xxg)
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
                     var(igrid,idim)=xxg(i)
                     igrid=next(igrid)
                  end do
                  idx=idx+ncache
               end if
            end do
         end do
      end do
    end subroutine pnetcdf_read_multilevel_real8_2D

    subroutine pnetcdf_read_multilevel_ncoarse_int_1D(f_id,name,var,ndims)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(out)::var
      integer,intent(in)::ndims

      integer(kind=mpi_offset_kind),dimension(2)::offsets,sizes
      integer::idx,igrid,ilevel,ibound,ind,ncache,istart,i,iskip
      integer::var_id
      integer,dimension(ncache_max)::iig

      status = nfmpi_inq_varid(f_id,name,var_id)
      sizes(2)=1
      do ind=1,ndims
         iskip=ncoarse+(ind-1)*ngridmax
         offsets(2)=ind
         do ilevel=1,nlevelmax
            idx=1
            offsets(1)=level_offsets(ilevel)+1
            sizes(1)=size_loc_by_level(ilevel)
            status = nfmpi_get_vara_int_all(f_id,var_id,offsets,sizes,iig)
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
                     var(igrid+iskip)=iig(i)
                     igrid=next(igrid)
                  end do
                  idx=idx+ncache
               end if
            end do
         end do
      end do
    end subroutine pnetcdf_read_multilevel_ncoarse_int_1D

    subroutine pnetcdf_write_part_int_1D(f_id,var_id,var,npart,offset)
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer,intent(in)::f_id
      integer,intent(in)::var_id
      integer,dimension(:),intent(in)::var
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::offset

      integer(kind=mpi_offset_kind),dimension(1)::sizes,starts
      integer::ipart,i
      integer,dimension(npart)::ii

      sizes(1)=npart
      starts(1)=offset
      ipart=0
      do i=1,npartmax
        if(levelp(i)>0)then
          ipart=ipart+1
          ii(ipart)=var(i)
        end if
      end do
      status = nfmpi_put_vara_int_all(f_id,var_id,starts,sizes,ii)
    end subroutine pnetcdf_write_part_int_1D

    subroutine pnetcdf_write_part_real8_1D(f_id,var_id,var,npart,offset)
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer,intent(in)::f_id
      integer,intent(in)::var_id
      real(kind=8),dimension(:),intent(in)::var
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::offset

      integer(kind=mpi_offset_kind),dimension(1)::sizes,starts
      integer::ipart,i
      real(kind=8),dimension(npart)::xdp

      sizes(1)=npart
      starts(1)=offset
      ipart=0
      do i=1,npartmax
        if(levelp(i)>0)then
          ipart=ipart+1
          xdp(ipart)=var(i)
        end if
      end do
      status = nfmpi_put_vara_double_all(f_id,var_id,starts,sizes,xdp)
    end subroutine pnetcdf_write_part_real8_1D

    subroutine pnetcdf_write_part_real8_2D(f_id,var_id,var,npart,offset)
      use amr_parameters,only:ndim
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer,intent(in)::f_id
      integer,intent(in)::var_id
      real(kind=8),dimension(:,:),intent(in)::var
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::offset

      integer(kind=mpi_offset_kind),dimension(2)::sizes,starts
      integer::ipart,i,idim
      real(kind=8),dimension(npart)::xdp

      sizes(1)=npart
      sizes(2)=1
      starts(1)=offset
      do idim=1,ndim
        ipart=0
        do i=1,npartmax
          if(levelp(i)>0)then
            ipart=ipart+1
            xdp(ipart)=var(i,idim)
          end if
        end do
        starts(2)=idim
        status = nfmpi_put_vara_double_all(f_id,var_id,starts,sizes,xdp)
      end do
    end subroutine pnetcdf_write_part_real8_2D

    subroutine pnetcdf_read_part_int_1D(f_id,name,var,npart,offset)
      use amr_parameters,only:ndim

      integer,intent(in)::f_id
      character(LEN=*),intent(in)::name
      integer,dimension(:),intent(out)::var
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::offset

      integer::var_id
      integer(kind=mpi_offset_kind),dimension(1)::sizes,offsets

      status = nfmpi_inq_varid(f_id,name,var_id)
      sizes(1)=npart
      offsets(1)=offset
      status = nfmpi_get_vara_int_all(f_id,var_id,offsets,sizes,var)
    end subroutine pnetcdf_read_part_int_1D

    subroutine pnetcdf_read_part_real8_1D(f_id,name,var,npart,offset)
      use amr_parameters,only:ndim

      integer,intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:),intent(out)::var
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::offset

      integer::var_id
      integer(kind=mpi_offset_kind),dimension(1)::sizes,offsets

      status = nfmpi_inq_varid(f_id,name,var_id)
      sizes(1)=npart
      offsets(1)=offset
      status = nfmpi_get_vara_double_all(f_id,var_id,offsets,sizes,var)
    end subroutine pnetcdf_read_part_real8_1D

    subroutine pnetcdf_read_part_real8_2D(f_id,name,var,npart,offset)
      use amr_parameters,only:ndim

      integer,intent(in)::f_id
      character(LEN=*),intent(in)::name
      real(kind=8),dimension(:,:),intent(out)::var
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::offset

      integer::var_id
      integer(kind=mpi_offset_kind),dimension(2)::sizes,offsets
      integer::idim

      status = nfmpi_inq_varid(f_id,name,var_id)
      sizes(1)=npart
      sizes(2)=1
      offsets(1)=offset
      do idim=1,ndim
        offsets(2)=idim
        status = nfmpi_get_vara_double_all(f_id,var_id,offsets,sizes,var(:,idim))
      end do
    end subroutine pnetcdf_read_part_real8_2D
end module ramses_pnetcdf
#endif
