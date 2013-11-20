! Copyright (c) 2007-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

#ifdef MPIIO
module ramses_mpiio
  use amr_commons,only:MPI_COMM_COMP
  use mpi

  implicit none

  interface mpiio_write_distrib_var
     !These subroutines can be done with mpi_file_write_ordered
     !but it is slower on Blue Gene/P (20090813)
     module procedure mpiio_write_distrib_int,      mpiio_write_distrib_int_1D,  &
                      mpiio_write_distrib_int_2D,   mpiio_write_distrib_real8,   &
                      mpiio_write_distrib_real8_1D, mpiio_write_distrib_real8_2D
  end interface

  interface mpiio_read_distrib_var
     !These subroutines can be done with mpi_file_read_ordered
     !but it is slower on Blue Gene/P (20090813)
     module procedure mpiio_read_distrib_int,      mpiio_read_distrib_int_1D,  &
                      mpiio_read_distrib_int_2D,   mpiio_read_distrib_real8,   &
                      mpiio_read_distrib_real8_1D, mpiio_read_distrib_real8_2D
  end interface


  contains
    subroutine mpiio_create_file(fileloc,f_id)
      use hints,only:hints_set

      character(LEN=*),intent(in) :: fileloc
      integer,intent(out)         :: f_id

      integer :: ierr,infos

      call MPI_Info_create(infos,ierr)

      call hints_set(infos)

      call MPI_File_open(MPI_COMM_COMP,fileloc,MPI_MODE_WRONLY+MPI_MODE_CREATE+MPI_MODE_UNIQUE_OPEN,& 
                         infos,f_id,ierr) 
      call MPI_File_set_size(f_id,int(0,kind=mpi_offset_kind),ierr) !Empties file if existed before

      call mpi_info_free(infos,ierr)

    end subroutine mpiio_create_file


    subroutine mpiio_open_file(fileloc,f_id)
      use hints,only:hints_set

      character(LEN=*),intent(in) :: fileloc
      integer,intent(out)         :: f_id

      integer :: ierr,infos

      call MPI_Info_create(infos,ierr)

      call hints_set(infos)

      call MPI_File_open(MPI_COMM_COMP,fileloc,MPI_MODE_RDONLY+MPI_MODE_UNIQUE_OPEN,&
                        infos,f_id,ierr)
      call MPI_Info_free(infos,ierr)

    end subroutine mpiio_open_file


    subroutine mpiio_write_distrib_int(fid,var,disp)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      integer,intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      offset=disp+(myid-1)*sizeofinteger
      call MPI_File_write_at_all(fid,offset,var,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*sizeofinteger
    end subroutine mpiio_write_distrib_int

    subroutine mpiio_write_distrib_real8(fid,var,disp)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      real(kind=8),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      offset=disp+(myid-1)*sizeofdouble
      call MPI_File_write_at_all(fid,offset,var,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*sizeofdouble
    end subroutine mpiio_write_distrib_real8

    subroutine mpiio_write_distrib_int_1D(fid,var,disp,count)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      integer,dimension(:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,optional,intent(in)::count

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::dims
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      if(present(count))then
        dims=count
      else
        dims=size(var)
      end if

      offset=disp+(myid-1)*dims*sizeofinteger
      call MPI_File_write_at_all(fid,offset,var,dims,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*dims*sizeofinteger
    end subroutine mpiio_write_distrib_int_1D

    subroutine mpiio_write_distrib_int_2D(fid,var,disp)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      integer,dimension(:,:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      offset=disp+(myid-1)*size(var)*sizeofinteger
      call MPI_File_write_at_all(fid,offset,var,size(var),MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*size(var)*sizeofinteger
    end subroutine mpiio_write_distrib_int_2D

    subroutine mpiio_write_distrib_real8_1D(fid,var,disp,count)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      real(kind=8),dimension(:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,optional,intent(in)::count

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::dims
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      if(present(count))then
        dims=count
      else
        dims=size(var)
      end if
      offset=disp+(myid-1)*dims*sizeofdouble
      call MPI_File_write_at_all(fid,offset,var,dims,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*dims*sizeofdouble
    end subroutine mpiio_write_distrib_real8_1D

    subroutine mpiio_write_distrib_real8_2D(fid,var,disp,count)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      real(kind=8),dimension(:,:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,optional,intent(in)::count

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::dims
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      if(present(count))then
        dims=count
      else
        dims=size(var)
      end if
      offset=disp+(myid-1)*dims*sizeofdouble
      call MPI_File_write_at_all(fid,offset,var,dims,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*size(var)*sizeofdouble
    end subroutine mpiio_write_distrib_real8_2D


    subroutine mpiio_read_distrib_int(fid,var,disp)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      integer,intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      offset=disp+(myid-1)*sizeofinteger
      call MPI_File_read_at_all(fid,offset,var,1,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*sizeofinteger
    end subroutine mpiio_read_distrib_int

    subroutine mpiio_read_distrib_real8(fid,var,disp)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      real(kind=8),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      offset=disp+(myid-1)*sizeofdouble
      call MPI_File_read_at_all(fid,offset,var,1,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*sizeofdouble
    end subroutine mpiio_read_distrib_real8

    subroutine mpiio_read_distrib_int_1D(fid,var,disp,count)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      integer,dimension(:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,optional,intent(in)::count

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::dims
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      if(present(count))then
        dims=count
      else
        dims=size(var)
      end if

      offset=disp+(myid-1)*dims*sizeofinteger
      call MPI_File_read_at_all(fid,offset,var,dims,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*dims*sizeofinteger
    end subroutine mpiio_read_distrib_int_1D

    subroutine mpiio_read_distrib_int_2D(fid,var,disp)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      integer,dimension(:,:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      offset=disp+(myid-1)*size(var)*sizeofinteger
      call MPI_File_read_at_all(fid,offset,var,size(var),MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*size(var)*sizeofinteger
    end subroutine mpiio_read_distrib_int_2D

    subroutine mpiio_read_distrib_real8_1D(fid,var,disp,count)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      real(kind=8),dimension(:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,optional,intent(in)::count

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::dims
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      if(present(count))then
        dims=count
      else
        dims=size(var)
      end if
      offset=disp+(myid-1)*dims*sizeofdouble
      call MPI_File_read_at_all(fid,offset,var,dims,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*dims*sizeofdouble
    end subroutine mpiio_read_distrib_real8_1D

    subroutine mpiio_read_distrib_real8_2D(fid,var,disp)
      use amr_commons,only:myid,ncpu

      integer,intent(in)::fid
      real(kind=8),dimension(:,:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      offset=disp+(myid-1)*size(var)*sizeofdouble
      call MPI_File_read_at_all(fid,offset,var,size(var),MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

      disp=disp+ncpu*size(var)*sizeofdouble
    end subroutine mpiio_read_distrib_real8_2D


    subroutine mpiio_write_multilevel_int_1D(fid,var,disp)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      integer,dimension(:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer::idx,igrid,ilevel,ibound,ncache,istart,i,sz
      integer,dimension(ncache_max)::iig
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      do ilevel=1,nlevelmax
         idx=1
         sz=0
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
        offset=disp+level_offsets(ilevel)*sizeofinteger
        call MPI_File_write_at_all(fid,offset,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
        !call MPI_File_write_ordered(fid,iig,size,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
      end do
      disp=disp+size_max*sizeofinteger

    end subroutine mpiio_write_multilevel_int_1D

    subroutine mpiio_write_multilevel_int_2D(fid,var,disp,ndims)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      integer,dimension(:,:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::ndims

      integer::idx,igrid,ilevel,ind,ibound,ncache,istart,i,sz
      integer,dimension(ncache_max)::iig
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      do ind=1,ndims
         do ilevel=1,nlevelmax
            idx=1
            sz=0
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
            offset=disp+((ind-1)*size_max+level_offsets(ilevel))*sizeofinteger
            call MPI_File_write_at_all(fid,offset,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            !call MPI_File_write_ordered(fid,iig,size,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         end do
      end do

      disp=disp+ndims*size_max*sizeofinteger
    end subroutine mpiio_write_multilevel_int_2D

    subroutine mpiio_write_multilevel_real8_2D(fid,var,disp,ndims)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      real(kind=8),dimension(:,:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::ndims

      integer::idx,idim,igrid,ilevel,ibound,ncache,istart,i,sz
      real(kind=8),dimension(ncache_max)::xdp
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

     do idim=1,ndims
         do ilevel=1,nlevelmax
            idx=1
            sz=0
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
            !call MPI_File_write_ordered(fid,xdp,size,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            offset=disp+((idim-1)*size_max+level_offsets(ilevel))*sizeofdouble
            call MPI_File_write_at_all(fid,offset,xdp,sz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
         end do
      end do

      disp=disp+ndims*size_max*sizeofdouble

    end subroutine mpiio_write_multilevel_real8_2D

    subroutine mpiio_write_multilevel_ncoarse_int_1D(fid,var,disp,ndims)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      integer,dimension(:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::ndims

      integer::idx,igrid,ilevel,ibound,ind,ncache,istart,i,iskip,sz
      integer,dimension(ncache_max)::iig
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      do ind=1,ndims
         iskip=ncoarse+(ind-1)*ngridmax
         do ilevel=1,nlevelmax
            idx=1
            sz=0
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
            offset=disp+((ind-1)*size_max+level_offsets(ilevel))*sizeofinteger
            call MPI_File_write_at_all(fid,offset,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            !call MPI_File_write_ordered(fid,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
         end do
      end do
      disp=disp+ndims*size_max*sizeofinteger

    end subroutine mpiio_write_multilevel_ncoarse_int_1D


    subroutine mpiio_read_multilevel_int_1D(fid,var,disp)
      use amr_commons,only:ncpu,next,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      integer,dimension(:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

      integer::idx,igrid,ilevel,ibound,ncache,istart,i,sz
      integer,dimension(ncache_max)::iig
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::ierr,sizeofinteger

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      do ilevel=1,nlevelmax
         idx=1
         sz=size_loc_by_level(ilevel)
         offset=disp+level_offsets(ilevel)*sizeofinteger
         call MPI_File_read_at_all(fid,offset,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
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

      disp=disp+size_max*sizeofinteger

    end subroutine mpiio_read_multilevel_int_1D

    subroutine mpiio_read_multilevel_int_2D(fid,var,disp,ndims,ngridmax2)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      integer,dimension(:,:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::ndims
      integer,intent(in)::ngridmax2

      integer::idx,igrid,ilevel,ind,ibound,ncache,istart,i,sz
      integer,dimension(ncache_max)::iig
      integer::tmp1,tmp2
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      do ind=1,ndims
         do ilevel=1,nlevelmax
            idx=1
            sz=size_loc_by_level(ilevel)
            !call MPI_File_read_ordered(fid,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            offset=disp+((ind-1)*size_max+level_offsets(ilevel))*sizeofinteger
            call MPI_File_read_at_all(fid,offset,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
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
                     var(igrid,ind)=iig(i)
                     igrid=next(igrid)
                  end do
                  idx=idx+ncache
               end if
            end do
         end do
      end do

      disp=disp+ndims*size_max*sizeofinteger

    end subroutine mpiio_read_multilevel_int_2D

    subroutine mpiio_read_multilevel_real8_2D(fid,var,disp,ndims)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      real(kind=8),dimension(:,:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::ndims

      integer::idx,igrid,ilevel,ind,ibound,ncache,istart,i,sz
      real(kind=8),dimension(ncache_max)::xxg
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofdouble,ierr

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      do ind=1,ndims
         do ilevel=1,nlevelmax
            idx=1
            sz=size_loc_by_level(ilevel)
            !call MPI_File_read_ordered(fid,xxg,sz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
            offset=disp+((ind-1)*size_max+level_offsets(ilevel))*sizeofdouble
            call MPI_File_read_at_all(fid,offset,xxg,sz,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
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
                     var(igrid,ind)=xxg(i)
                     igrid=next(igrid)
                  end do
                  idx=idx+ncache
               end if
            end do
         end do
      end do

      disp=disp+ndims*size_max*sizeofdouble

    end subroutine mpiio_read_multilevel_real8_2D

    subroutine mpiio_read_multilevel_ncoarse_int_1D(fid,var,disp,ndims)
      use amr_commons,only:ncoarse,ncpu,next,ngridmax,headb,headl,numbb,numbl
      use amr_parameters,only:nboundary,nlevelmax
      use io_commons

      integer,intent(in)::fid
      integer,dimension(:),intent(out)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::ndims

      integer::idx,igrid,ilevel,ind,ibound,ncache,istart,i,iskip,sz
      integer,dimension(ncache_max)::iig
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      do ind=1,ndims
         iskip=ncoarse+(ind-1)*ngridmax
         do ilevel=1,nlevelmax
            idx=1
            sz=size_loc_by_level(ilevel)
            !call MPI_File_read_ordered(fid,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
            offset=disp+((ind-1)*size_max+level_offsets(ilevel))*sizeofinteger
            call MPI_File_read_at_all(fid,offset,iig,sz,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
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
      disp=disp+ndims*size_max*sizeofinteger

    end subroutine mpiio_read_multilevel_ncoarse_int_1D

    subroutine mpiio_write_part_int_1D(f_id,var,disp,npart,npart_tot,offset_part)
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer,intent(in)::f_id
      integer,dimension(:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::npart_tot
      integer(kind=mpi_offset_kind),intent(in)::offset_part

      integer::i,ipart
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofinteger,ierr
      integer,dimension(npart)::ii

      call MPI_Type_size(MPI_INTEGER,sizeofinteger,ierr)

      ipart=0
      do i=1,npartmax
         if(levelp(i)>0)then
            ipart=ipart+1
            ii(ipart)=var(i)
         end if
      end do
      offset=disp+offset_part*sizeofinteger
      call MPI_File_write_at_all(f_id,offset,ii,npart,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

      disp=disp+npart_tot*sizeofinteger
    end subroutine mpiio_write_part_int_1D

    subroutine mpiio_write_part_real8_1D(f_id,var,disp,npart,npart_tot,offset_part)
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer,intent(in)::f_id
      real(kind=8),dimension(:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::npart_tot
      integer(kind=mpi_offset_kind),intent(in)::offset_part

      integer::i,ipart
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofdouble,ierr
      real(kind=8),dimension(npart)::xdp

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      ipart=0
      do i=1,npartmax
         if(levelp(i)>0)then
            ipart=ipart+1
            xdp(ipart)=var(i)
         end if
      end do
      offset=disp+offset_part*sizeofdouble
      call MPI_File_write_at_all(f_id,offset,xdp,npart,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)

      disp=disp+npart_tot*sizeofdouble
    end subroutine mpiio_write_part_real8_1D

    subroutine mpiio_write_part_real8_2D(f_id,var,disp,npart,npart_tot,offset_part)
      use amr_parameters,only:ndim
      use pm_commons,only:levelp
      use pm_parameters,only:npartmax

      integer,intent(in)::f_id
      real(kind=8),dimension(:,:),intent(in)::var
      integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
      integer,intent(in)::npart
      integer(kind=mpi_offset_kind),intent(in)::npart_tot
      integer(kind=mpi_offset_kind),intent(in)::offset_part

      integer::i,idim,ipart
      integer(kind=MPI_OFFSET_KIND)::offset
      integer::sizeofdouble,ierr
      real(kind=8),dimension(npart)::xdp

      call MPI_Type_size(MPI_DOUBLE_PRECISION,sizeofdouble,ierr)

      do idim=1,ndim
         ipart=0
         do i=1,npartmax
            if(levelp(i)>0)then
               ipart=ipart+1
               xdp(ipart)=var(i,idim)
            end if
        end do
        offset=disp+((idim-1)*npart_tot+offset_part)*sizeofdouble
        call MPI_File_write_at_all(f_id,offset,xdp,npart,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
      end do

      disp=disp+ndim*npart_tot*sizeofdouble
    end subroutine mpiio_write_part_real8_2D

end module ramses_mpiio
#endif
