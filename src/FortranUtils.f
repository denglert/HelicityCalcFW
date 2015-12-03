***********************************************************************

      real*8 function LorentzScalarProd(p1,p2)
      real*8 p1(0:3),p2(0:3)

      LorentzScalarProd=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)

      RETURN
      END

***********************************************************************

      complex*16 function ComplexLorentzScalarProd(p1,p2)

      complex*16 p1(0:3)
      real*8 p2(0:3)

      ComplexLorentzScalarProd=p1(0)*p2(0)-p1(1)*p2(1)
     & -p1(2)*p2(2)-p1(3)*p2(3)

      RETURN 
      END

***********************************************************************

      real*8 function CalcP3Mag(p)

      real*8 p(0:3)

      CalcP3Mag=p(1)*p(1)+p(2)*p(2)+p(3)*p(3)
      CalcP3Mag=sqrt(CalcP3Mag)

      RETURN 
      END


***********************************************************************

      real*8 function CalcP4MagSq(p)

      real*8 p(0:3)

      CalcP4MagSq=p(0)*p(0)-p(1)*p(1)-p(2)*p(2)-p(3)*p(3)

      RETURN 
      END


***********************************************************************

      subroutine Constructk0(k0,p)
      real*8 k0(0:3), p(0:3)
      real*8 pmag
      pmag = CalcP3Mag(p)

      k0(0) = 1
      do i=1,3
      k0(i) = -p(i)/pmag
      enddo

      do nu=0,3
      k0(nu) = k0(nu)/(p(0)+pmag)
      enddo

      END


***********************************************************************

      subroutine ConstructSpinPolarization(s,p,k0)
      real*8 s(0:3), p(0:3), k0(0:3)
      real*8 LorentzScalarProd

      real*8 m
      real*8 pk0


      m=LorentzScalarProd(p,p)
      m=sqrt(m)
      pk0=LorentzScalarProd(p,k0)

      do mu=0,3
      s(mu)=p(mu)/m-(m/pk0)*k0(mu)
      enddo

      END


***********************************************************************

      subroutine PrintLorentzVector(p)
      real*8 p(0:3)
      do nu=0,3
      print*, nu, p(nu)
      enddo
      END
