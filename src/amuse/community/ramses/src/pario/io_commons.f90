! Copyright (c) 2006-2011 IDRIS/CNRS
! Author: Philippe Wautelet (IDRIS/CNRS), wautelet@idris.fr
! Distributed under the CeCILL 2.0 license. For full terms see the file LICENSE.

module io_parameters
  integer,parameter::TAG_BAK_AMR=1000
  integer,parameter::TAG_BAK_HYD=1001
  integer,parameter::TAG_BAK_PAR=1002
  integer,parameter::TAG_BAK_POI=1003

  integer,parameter::TAG_OUT_INF=2000
  integer,parameter::TAG_OUT_AMR=2001
  integer,parameter::TAG_OUT_HYD=2002
  integer,parameter::TAG_OUT_COO=2003
  integer,parameter::TAG_OUT_PAR=2004
  integer,parameter::TAG_OUT_STA=2005

  integer,parameter::MAXFILES=11 ! max files/computational process

  integer,parameter::ILUN_IO=100111

end module io_parameters

module io_commons
  use io_parameters
  use amr_parameters,only:MAXLINE
  use mpi,only:mpi_offset_kind

  integer,allocatable,dimension(:,:,:)::numblio
  integer::count_out=0,count_bak=0
  integer::nbfiles=0
  character(LEN=MAXLINE)::scratchdir='.',permdir='.'
  character(LEN=MAXLINE),dimension(:),allocatable::filelist

  real(kind=8),parameter::multiplier=1.3 ! Size multiplier to determine buffer_size from ncache_max
  integer::ncache_max
  integer::buffer_size
  integer(kind=mpi_offset_kind)::size_max
  integer(kind=mpi_offset_kind),dimension(:),allocatable::level_offsets
  integer(kind=mpi_offset_kind),dimension(:),allocatable::size_loc_by_level
  logical,dimension(:),allocatable::outputs

end module io_commons

