!!! wien2wannier/SRC_w2w/abc.f

module     abc_m; contains
subroutine abc(l,jatom,pei,pi12lo,pe12lo,jlo,lapw)
  use param,  only: unit_out, Nrf
  use struct, only: RMT
  use loabc,  only: alo
  use atspdt, only: P, DP
  use const,  only: R8

  implicit none

  integer,  intent(in) :: l, jatom, jlo
  real(R8), intent(in) :: pei, pi12lo, pe12lo
  logical,  intent(in) :: lapw

  integer  :: irf, jrf
  real(R8) :: xac, xbc, clo, alonorm
  !---------------------------------------------------------------------
  !
  do irf=1,nrf
     alo(l,jlo,irf)=0.d0
  end do
  if (lapw) then
     irf=2+jlo
     xac=p(l,irf)*dp(l,2)-dp(l,irf)*p(l,2)
     xac=xac*rmt(jatom)*rmt(jatom)
     xbc=p(l,irf)*dp(l,1)-dp(l,irf)*p(l,1)
     xbc=-xbc*rmt(jatom)*rmt(jatom)
     clo=xac*(xac+2.0D0*pi12lo)+xbc* &
          (xbc*pei+2.0D0*pe12lo)+1.0D0
     clo=1.0D0/sqrt(clo)
     write(unit_out,*)clo
     if (clo.gt.2.0D2) clo=2.0d2
     alo(l,jlo,1)=clo*xac
     alo(l,jlo,2)=clo*xbc
     alo(l,jlo,irf)=clo
     write(unit_out,10) l, alo(l,jlo,1), alo(l,jlo,2), alo(l,jlo,irf)
  else
     if (jlo.eq.1) then
        alonorm=sqrt(1.d0+(P(l,1)/P(l,2))**2 * PEI)
        alo(l,jlo,1)=1.d0/alonorm
        alo(l,jlo,2)=-P(l,1)/P(l,2)/alonorm
     else
        xbc=-P(l,1)/P(l,1+jlo)
        xac=sqrt(1+xbc**2+2*xbc*PI12LO)
        alo(l,jlo,1)=1.d0/xac
        alo(l,jlo,1+jlo)=xbc/xac
     end if
     write (unit_out,10)l,(alo(l,jlo,jrf),jrf=1,nrf)
  end if
10 FORMAT ('LO COEFFICIENT: l,A,B,C  ',i2,5X,6F12.5)
  return
end subroutine abc
end module     abc_m


!!/---
!! Local Variables:
!! mode: f90
!! End:
!!\---
!!
!! Time-stamp: <2016-07-18 15:10:05 assman@faepop71.tu-graz.ac.at>
