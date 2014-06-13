subroutine FortranAssert(prepost,expression,filename,linenum)
  character(*), intent(in) :: prepost,expression, filename
  integer,      intent(in) :: linenum
  write(6,'(a,a,a,a,a,a,i5)') prepost," (", expression, ") failed in file ", filename," line ",linenum
  stop
end subroutine FortranAssert

subroutine FortranTrace(expression,filename,linenum)
  implicit none
  character(*), intent(in) :: expression, filename
  integer,      intent(in) :: linenum
  write(6,'(a,a,a,a,a,i5)') "Trace '", expression, "' in file ", filename," line ",linenum
end subroutine FortranTrace

subroutine FortranTraceInteger(iname,value,filename,linenum)
  implicit none
  character(*), intent(in) :: iname, filename
  integer,      intent(in) :: value
  integer,      intent(in) :: linenum
  write(6,'(a,a,a,i5,a,a,a,i5)') "Trace ", iname, " = ", value, " in file ", filename," line ",linenum
end subroutine FortranTraceInteger

subroutine FortranTraceReal(iname,value,filename,linenum)
  use real_kind
  implicit none
  character(*), intent(in) :: iname, filename
  real(double), intent(in) :: value
  integer,      intent(in) :: linenum
  write(6,'(1p,a,a,a,d24.16,a,a,a,i5)') "Trace ", iname, " = ", value, " in file ", filename," line ",linenum
end subroutine FortranTraceReal

subroutine FortranNanCheck(iname,value,filename, linenum)
   use ieee_arithmetic, only: isnan=>ieee_is_nan
   use real_kind
   implicit none
   character(*), intent(in) :: iname, filename
   real(double), intent(in) :: value
   integer,      intent(in) :: linenum
   if (isnan(value)) then
      write(6,'(1p,a,a,a,d24.16,a,a,a,i5)') "NaN check: ", iname, " = ", value, " in file ", filename," line ",linenum
      stop
   end if
end subroutine FortranNanCheck
