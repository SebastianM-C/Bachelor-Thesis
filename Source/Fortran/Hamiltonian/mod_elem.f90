!    Compute the matrix element on harmonic oscillator states
!             + k    l
!     < m | (a )  (a)  | n >
!

module mod_elem

  implicit none
  contains
    pure function elem(m, n, k, l)

      use, intrinsic :: iso_fortran_env
      implicit none

      integer(kind=int32), intent(in) :: m, n, k, l
      !real(kind=REAL128) :: elem
      real(kind=real64) :: elem

      integer(kind=int32) :: i

      elem = 0.

      if(m /= n - l + k) then
        return
      endif
      if(n < l) then
        return
      endif
      elem = 1.
      ! Apply creation operator
      do i = 0, l - 1
        !elem = elem * sqrt(real(n - i, 16))
        elem = elem * sqrt(real(n - i, 8))
      end do
      ! Apply anihilation operator
      do i = 1, k
        !elem = elem * sqrt(real(n - l + i, 16))
        elem = elem * sqrt(real(n - l + i, 8))
      end do

    end function elem

end module mod_elem
