subroutine AO_to_MO(nBas,c,ERI_AO,ERI_MO)

  implicit none

! Input variables
  integer, intent(in) :: nBas
  double precision, intent(in) :: c(nBas, nBas)
  double precision, intent(in) :: ERI_AO(nBas, nBas, nBas, nBas)

! Local variables
  integer :: mu, nu, la, si
  integer :: p, q, r, s
  double precision, allocatable :: scr(:,:,:,:), tmp1(:,:,:,:), tmp2(:,:,:,:), tmp3(:,:,:,:)
  integer :: alloc_status

! Output variables
  double precision, intent(out) :: ERI_MO(nBas, nBas, nBas, nBas)

! Memory allocation with error checking
  allocate(scr(nBas, nBas, nBas, nBas), stat=alloc_status)
  if (alloc_status /= 0) stop "Error: Memory allocation failed for scr"
  allocate(tmp1(nBas, nBas, nBas, nBas), stat=alloc_status)
  if (alloc_status /= 0) stop "Error: Memory allocation failed for tmp1"
  allocate(tmp2(nBas, nBas, nBas, nBas), stat=alloc_status)
  if (alloc_status /= 0) stop "Error: Memory allocation failed for tmp2"
  allocate(tmp3(nBas, nBas, nBas, nBas), stat=alloc_status)
  if (alloc_status /= 0) stop "Error: Memory allocation failed for tmp3"

! initiallisation
  scr = 0.0d0
  tmp1 = 0.0d0
  tmp2 = 0.0d0
  tmp3 = 0.0d0

! AO to MO transformation
  do p = 1, nBas
    do mu = 1, nBas
      tmp1(p,:,:,:) = tmp1(p,:,:,:) + c(mu,p)*ERI_AO(mu,:,:,:)
    end do
    do q = 1, nBas
      do nu = 1, nBas
        tmp2(p,q,:,:) = tmp2(p,q,:,:) + c(nu,q)*tmp1(p,nu,:,:)
      end do
      do r = 1, nBas
        do la = 1, nBas
          tmp3(p,q,r,:) = tmp3(p,q,r,:) + c(la,r)*tmp2(p,q,la,:)
        end do
        do s = 1, nBas
          do si = 1, nBas
            scr(p,q,r,s) = scr(p,q,r,s) + c(si,s)*tmp3(p,q,r,si)
          end do
        end do
      end do
    end do
  end do

  ERI_MO = scr

! Deallocate memory
  deallocate(scr, tmp1, tmp2, tmp3)

end subroutine AO_to_MO
