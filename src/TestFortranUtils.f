      program test_fortranutils
      real*8 LorentzScalarProd
      real*8 CalcP3Mag
      complex*16 ComplexLorentzScalarProd

      real*8 p(0:3),p1(0:3),p2(0:3)
      real*8 k0(0:3),k0_1(0:3),k0_2(0:3)

      real*8 prod1
      real*8 p3mag1
      real*8 p4magsq1
      real*8 sp

      real*8 realvec(0:3)
      complex*16 cvec(0:3)
      complex*16 complexprod1

      print*, 'Test FortanUtils suite'

      p1(0)=6
      p1(1)=0
      p1(2)=3
      p1(3)=4

      p2(0)=0
      p2(1)=1
      p2(2)=1
      p2(3)=3

      print*, 'p1'
      call PrintLorentzVector(p1)

      print*, 'p2'
      call PrintLorentzVector(p2)

      prod1 = LorentzScalarProd(p1,p2)
      print*, 'LorentzScalarProd', prod1

      p3mag1 = CalcP3Mag(p1)
      print*, 'p3mag1', p3mag1

      call Constructk0(k0_1,p1)
      call Constructk0(k0_2,p2)

      
      print*, 'k0_1'
      call PrintLorentzVector(k0_1)

      print*, 'k0_2'
      call PrintLorentzVector(k0_2)

      p4magsq1=CalcP4MagSq(k0_1)
      print*, 'k0_1^2:', p4magsq1


      print*, 'Testing ComplexLorentzScalarProd function:'

      realvec(0)=1.
      realvec(1)=0.
      realvec(2)=1.
      realvec(3)=0.5

      cvec(0)=(1.0,2.0)
      cvec(1)=(0.0,0.0)
      cvec(2)=(-0.5,1.0)
      cvec(3)=(1.0,-2.0)

      complexprod1=ComplexLorentzScalarProd(cvec,realvec)

      print*, 'ComplexLorentzScalarProd', complexprod1

      print*, ''
      print*, 'Testing ConstructSpinPolarization'
      print*, ''

      p(0)=5.0
      p(1)=4.0
      p(2)=0.0
      p(3)=0.0


      call Constructk0(k0,p)
      call ConstructSpinPolarization(s,p,k0)

      print*, 'p'
      call PrintLorentzVector(p)
      print*, 'k0'
      call PrintLorentzVector(k0)
      print*, 's'
      call PrintLorentzVector(s)

      sp=LorentzScalarProd(p,s)
      print*, 'sp', sp

      end program test_fortranutils
