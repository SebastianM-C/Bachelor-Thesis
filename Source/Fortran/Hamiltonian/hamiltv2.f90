! Compute the matrix elements of a boson Hamiltonian

program hamilt

  use, intrinsic :: iso_fortran_env
  use mod_elem
  implicit none

  integer(kind=int32) :: n1, n2, m1, m2, ios, err
  integer(kind=int64) :: n, nn, i, j, k, reclen
  !real(kind=REAL128) :: a, b, d
  !real(kind=REAL128), dimension(:, :), allocatable :: h
  real(kind=real64) :: a, b, d
  real(kind=real64), dimension(:, :), allocatable :: h
  integer(kind=int32), dimension(:, :), allocatable :: index

  ! Read Hamiltonian parameters
  open(unit=10, file='hamilt.inp', iostat=ios)
  if ( ios /= 0 ) stop "hamilt.inp: Error opening file "
  read(10, *) n, a, b, d
  close(unit=10, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 10"

  nn = n * (n + 1) / 2
  ! Allocate arrays
  allocate(h(nn, nn), stat=err)
  if (err /= 0) print *, "h: Allocation request denied"
  allocate(index(2, nn), stat=err)
  if (err /= 0) print *, "h: Allocation request denied"

  inquire(iolength=reclen) h(:, 1)
  !write(*, fmt=*) reclen
  open(unit=20, file='hamilt.bin', iostat=ios, form='unformatted', &
       access='direct', recl=reclen)
  if ( ios /= 0 ) stop "hamilt.bin: Error opening file "

  ! Compute index
  k = 0
  do i = 0, n - 1
    do j = 0, n - i - 1
      k = k + 1
      index(1, k) = i
      index(2, k) = j
      !write(*,"(3i3x)") i, j, k
    enddo
  end do
  !write(*,*) "------"
  !write(*, "(2i3x)") (index(1,i), index(2, i), i=1, nn)
  !write(*,*) "------"

  do i = 1, nn
    do j = 1, nn
      m1 = index(1, i)
      m2 = index(2, i)
      n1 = index(1, j)
      n2 = index(2, j)

      h(i, j) = a * (elem(m1, n1, 1, 1) * elem(m2, n2, 0, 0) +          &
                     elem(m1, n1, 0, 0) * elem(m2, n2, 1, 1))           &
              + 0.25 * b * (3 * elem(m1, n1, 1, 0) * elem(m2, n2, 2, 0) &
              + 3 * elem(m1, n1, 0, 1) * elem(m2, n2, 0, 2)             &
              - elem(m1, n1, 3, 0) * elem(m2, n2, 0, 0)                 &
              - elem(m1, n1, 0, 3) * elem(m2, n2, 0, 0))                &
              + 0.75 * b * (elem(m1, n1, 0, 1) * elem(m2, n2, 2, 0)     &
              + elem(m1, n1, 1, 0) * elem(m2, n2, 0, 2)                 &
              - elem(m1, n1, 1, 2) * elem(m2, n2, 0, 0)                 &
              - elem(m1, n1, 2, 1) * elem(m2, n2, 0, 0)                 &
              + 2 * elem(m1, n1, 0, 1) * elem(m2, n2, 1, 1)             &
              + 2 * elem(m1, n1, 1, 0) * elem(m2, n2, 1, 1))            &
              + 0.375 * d * (elem(m1, n1, 2, 2) * elem(m2, n2, 0, 0)    &
              + elem(m1, n1, 0, 0) * elem(m2, n2, 2, 2))                &
              + 0.125 * d * (elem(m1, n1, 2, 0) * elem(m2, n2, 0, 2)    &
              + elem(m1, n1, 0, 2) * elem(m2, n2, 2, 0))                &
              + 0.500 * d * elem(m1, n1, 1, 1) * elem(m2, n2, 1, 1)     &
              + 0.250 * d * (elem(m1, n1, 1, 3) * elem(m2, n2, 0, 0)    &
              + elem(m1, n1, 3, 1) * elem(m2, n2, 0, 0)                 &
              + elem(m1, n1, 0, 0) * elem(m2, n2, 1, 3)                 &
              + elem(m1, n1, 0, 0) * elem(m2, n2, 3, 1)                 &
              + elem(m1, n1, 0, 2) * elem(m2, n2, 1, 1)                 &
              + elem(m1, n1, 2, 0) * elem(m2, n2, 1, 1)                 &
              + elem(m1, n1, 1, 1) * elem(m2, n2, 0, 2)                 &
              + elem(m1, n1, 1, 1) * elem(m2, n2, 2, 0))                &
              + 0.0625 * d * (elem(m1, n1, 4, 0) * elem(m2, n2, 0, 0)   &
              + elem(m1, n1, 0, 4) * elem(m2, n2, 0, 0)                 &
              + elem(m1, n1, 0, 0) * elem(m2, n2, 4, 0)                 &
              + elem(m1, n1, 0, 0) * elem(m2, n2, 0, 4)                 &
              + 2 * elem(m1, n1, 2, 0) * elem(m2, n2, 2, 0)             &
              + 2 * elem(m1, n1, 0, 2) * elem(m2, n2, 0, 2))
    enddo
  enddo

  do i = 1, nn
  !  write(*, '(1000(F35.32))') (h(i,j), j=1,nn)
     write(20, rec=i) (h(i,j), j=1,nn)
  enddo

  deallocate(h)
  deallocate(index)
  close(unit= 20, iostat=ios)
  if ( ios /= 0 ) stop "Error closing file unit 20"

end program hamilt
