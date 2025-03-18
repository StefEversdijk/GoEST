subroutine TDHF(nBas, nO, ERI_MO, ex_en, e)

  ! TDHF module

  implicit none

  ! Input variables
  integer, intent(in)                 :: nBas
  integer, intent(in)                 :: nO
  double precision, intent(in)        :: ERI_MO(nBas, nBas, nBas, nBas)
  double precision, intent(in)        :: e(nBas)

  ! Local variables
  integer                             :: p, q, r, s, i, j, a, b, delta, size, vir, i_tmp, j_tmp
  double precision, allocatable       :: matr_a(:,:), matr_b(:,:), matr_dif(:,:), matr_sum(:,:), matr_dif_eigen(:),tmp_matr_dif(:,:), tmp_matr_dif_eigen(:,:), matr_c(:,:), tmp(:)
  integer                             :: kronecker_delta

  ! Output variables
  double precision, intent(out) :: ex_en

  ! Hello world
  write(*,*)
  write(*,*) '**************************************'
  write(*,*) '|          TDHF calculation          |'
  write(*,*) '**************************************'
  write(*,*)

  ! Variables indicating the size of A and B matrix
  vir = (nBas - nO)
  size = vir * nO

  ! Allocation of TDHF matrices and temporary intermediates
  allocate(matr_a(size, size))
  allocate(matr_b(size, size))
  allocate(matr_dif(size,size))
  allocate(matr_sum(size,size))
  allocate(matr_dif_eigen(size))
  allocate(tmp_matr_dif(size,size))
  allocate(tmp_matr_dif_eigen(size,size))
  allocate(matr_c(size,size))
  allocate(tmp(size))

  ! Initialisation of matrices
  matr_a(:,:) = 0.0d0
  matr_b(:,:) = 0.0d0
  tmp_matr_dif_eigen(:,:) = 0.0d0

  ! Construct A and B matrix
  do i = 1, size
    a = mod(i-1, vir) + nO + 1 ! Virtual part
    i_tmp = (i-1)/vir + 1 ! Occupied part

    do j = 1, size
      b = mod(j-1, vir) + nO + 1 ! Virtual part
      j_tmp = (j-1)/vir + 1 ! Occupied part

      if (i == j) then ! Diagonal elements
        matr_a(i, j) = (e(a) - e(i_tmp))
      end if

      matr_a(i, j) = matr_a(i, j) + 2 * ERI_MO(i_tmp, b, a, j_tmp) - ERI_MO(i_tmp, b, j_tmp, a) ! Overlap elements for A
      matr_b(i, j) = matr_b(i, j) + 2 * ERI_MO(i_tmp, j_tmp, a, b) - ERI_MO(i_tmp, j_tmp, b, a) ! Overlap elements for B

    end do
  end do

 matr_dif(:,:) = matr_a(:,:) - matr_b(:,:) ! A-B matrix
 matr_sum(:,:) = matr_a(:,:) + matr_b(:,:) ! A+B Matrix

 tmp_matr_dif = matr_dif

 call diagonalize_matrix(size, tmp_matr_dif, matr_dif_eigen) ! Diagonalisation of the A-B matrix

 do i = 1, size
   tmp_matr_dif_eigen(i,i) = sqrt(matr_dif_eigen(i)) ! square root of the eigenvalue
 end do

 matr_dif = matmul(tmp_matr_dif,matmul(tmp_matr_dif_eigen,transpose(tmp_matr_dif))) ! Obtain the (A-B)^(-1/2) matrix

 matr_c = matmul(matr_dif,matmul(matr_sum,matr_dif)) ! Obtain the Hermitian formulation

 ! Diagonalize matrix
 call diagonalize_matrix(size, matr_c, tmp)

 tmp = sqrt(tmp)
  
 !Print excitation energies
 write(*,*) "Excitation energies:"
 do i = 1, size
   write(*, '(I6, F16.8)') i, tmp(i)
 end do

end subroutine TDHF
