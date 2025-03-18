subroutine CIS(nBas, nO, ERI_MO, ex_en, e)

  ! CIS module

  implicit none

  ! Input variables
  integer, intent(in)                 :: nBas
  integer, intent(in)                 :: nO
  double precision, intent(in)        :: ERI_MO(nBas, nBas, nBas, nBas)
  double precision, intent(in)        :: e(nBas)

  ! Local variables
  integer                             :: p, q, r, s, i, j, a, b, delta, size, vir, i_tmp, j_tmp
  double precision, allocatable       :: matr(:,:), tmp(:)
  integer                             :: kronecker_delta

  ! Output variables
  double precision, intent(out) :: ex_en

  ! Hello world
  write(*,*)
  write(*,*) '**************************************'
  write(*,*) '|          CIS calculation           |'
  write(*,*) '**************************************'
  write(*,*)

  
  ! Variables indicating size of CIS matrix
  vir = (nBas - nO)
  size = vir * nO

  ! Allocation of CIS matrix
  allocate(matr(size, size))
  allocate(tmp(size))

  ! Initiallisation of CIS matrix
  matr(:,:) = 0.0d0

  ! Construct CIS matrix
  do i = 1, size
    a = mod(i-1, vir) + nO + 1 ! Virtual part
    i_tmp = (i-1)/vir + 1 ! Occupied part

    do j = 1, size
      b = mod(j-1, vir) + nO + 1 ! Virtual part
      j_tmp = (j-1)/vir + 1 ! occupied part

      if (i == j) then ! Diagonal elements
        matr(i, j) = (e(a) - e(i_tmp))
      end if

      matr(i, j) = matr(i, j) + 2 * ERI_MO(i_tmp, b, a, j_tmp) - ERI_MO(i_tmp, b, j_tmp, a) ! Overlap elements

    end do
  end do

  ! Diagonalize matrix
  call diagonalize_matrix(size, matr, tmp)

  ! Store excitation energies
  !ex_en = tmp ! This didn't work for some reason

  ! Print excitation energies
  write(*,*) "Excitation energies:"
  do i = 1, size
     write(*, '(I6, F16.8)') i, tmp(i)
  end do


end subroutine CIS
