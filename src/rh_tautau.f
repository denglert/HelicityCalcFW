c Function rh_tatau
c Input 'p1' and 'p2'
c Returns 8 byte real value 
      real*8 function rh_tautau(p1,p2,cth)

c Variable a-b, d-h o-z are all real
c Variable 'c' is complex, and i,j,k,l,m,n are not defined
      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)

c p1 and p2 are
      dimension p1(0:3),p2(0:3)

      dimension cth(2,2)

C Define a structure type 'tc' wicth member a(2,2), b(2,2) and c(2,2), 2x2 matrices
      structure /tc/
        double complex a(2,2),b(2,2),c(2,2)
      end structure
C Define a variable called 'tsc' of type 'tc'
      record/tc/tsc

c i1 and i2 define the polarizations 
      integer i1, i2

c Make 'rmtau' a common variable
      real*8 rmtau
      COMMON/masses/rmtau

c Make 'czero' and 'cim' two component variables wicth initialization
      PARAMETER (czero=(0.d0,0.d0),cim=(0.d0,1.d0))

c     Check if rmtau is cthe tau mass
       print *,'Inside rh_tautau FORTRAN func., value of rmtau:', rmtau

c     Fortran momenta
       do i=0,3                          
          print *,'FORTRAN p1(i):',  p1(i), i
       enddo

       do i=0,3                          
          print *,'FORTRAN p2(i):',  p2(i), i
       enddo

C ******************************************************************
c The product of p1 and k0 and p2 and k0
c where k0 = (1,1,0,0)

      factor = 1
      p1k0=factor*(p1(0)-p1(1))
      p2k0=factor*(p2(0)-p2(1))

c subroutine TH(qu,qd,abc)
c computes the coefficients a,b,c (d=0) for the insertion of a scalar particle
c into a fermion line.
* TSC -- qu=p1,qd=p2,a=tsc.a,b=tsc.b,c=tsc.c,cr=1.d0,cl=1.d0                    
c subroutine TSC(qu,qd,abc,cr,cl,den)
c computes the coefficients a,b,c (d=0) for the insertion of a pseudoscalar
c particle into a fermion line.
      auxa=-p1k0*p2(2)+p2k0*p1(2)
      cauxa=-cim*(p2(3)*p1k0-p1(3)*p2k0)
      tsc.a(1,2)=(auxa+cauxa)
      tsc.a(2,1)=(-auxa+cauxa)
      tsc.b(1,1)=p2k0
      tsc.b(2,2)=p2k0
      tsc.c(1,1)=p1k0
      tsc.c(2,2)=p1k0
* mline -- res=cth(&1,&2),abcd=tsc.,m1=rmtau,m2=(-rmtau),den=0,nsum=0            
c subroutine mline(mline,abcd,m1,m2,den,nsum)
c computes the final massive fermion line A+m1*B+m2*C+m1*m2*D.
c The name of the final result mline must contain the fermion indices specified
c in the abcd coefficients. For example, if abcd=name(i1,i2) so that the
c a-coefficient is given by name(i1,i2).a(i3,i4), then the mline must be named as
c mline(i1?,i2?,&1,&2). In the output file, mline will appear as
c mline(i1,i2,i3,i4) within a do loop over all 4 indices.
c If one wants to change the order of the indices in mline, e.g. mline(i1,i4,i2,i3),
c the input must be the following: mline(i1?,&2,i2?,&1).
      do iut=1,2                          
      do jut=1,2
      cth(iut,jut)=tsc.a(iut,jut)+rmtau*tsc.b(iut,jut)+(-rmtau)*t
     & sc.c(iut,jut)
      enddo
      enddo

      res=0.d0

       do i1=1,2
         do i2=1,2
           cres=cth(i1,i2)                  
           res=res+dreal(cres)**2+dimag(cres)**2
         enddo
       enddo


      print *,'FORTRAN res:', res
c Just one polarization

c     cres=cth(i1,i2)                  
c     res=res+dreal(cres)**2+dimag(cres)**2

      rh_tautau=res/p1k0/p2k0

      print *,'FORTRAN rh_tautau:', rh_tautau

      RETURN
      END
