c Function rh_tatau
c Input 'p1' and 'p2'
c Returns 8 byte real value 
      real*8 function rh_tautau(p1,p2)

c Variable a-b, d-h o-z are all real
c Variable 'c' is complex, and i,j,k,l,m,n are not defined
      implicit real*8 (a-b,d-h,o-s,u-z)
      implicit double complex (c)

c p1 and p2 are
      dimension p1(0:3),p2(0:3)

      dimension th(2,2)

C Define a type type 'tc' with member a(2,2), b(2,2) and c(2,2), 2x2 matrices
      type tc
        double complex a(2,2),b(2,2),c(2,2)
      end type
C Define a variable called 'tsc' of type 'tc'
      type(tc)tsc

c Make 'rmtau' a common variable
      real rmtau
      COMMON/masses/rmtau

c Make 'czero' and 'cim' two component variables with initialization
      PARAMETER (czero=(0.d0,0.d0),cim=(0.d0,1.d0))

c     Check if rmtau is the tau mass
      print *,'rmtau', rmtau

C ******************************************************************
c The product of p1 and k0 and p2 and k0
c where k0 = (1,1,0,0)
      p1k0=p1(0)-p1(1)
      p2k0=p2(0)-p2(1)

* TSC -- qu=p1,qd=p2,a=tsc.a,b=tsc.b,c=tsc.c,cr=1.d0,cl=1.d0                    
      auxa=-p1k0*p2(2)+p2k0*p1(2)
      cauxa=-cim*(p2(3)*p1k0-p1(3)*p2k0)
      tsc%a(1,2)=(auxa+cauxa)
      tsc%a(2,1)=(-auxa+cauxa)
      tsc%b(1,1)=p2k0
      tsc%b(2,2)=p2k0
      tsc%c(1,1)=p1k0
      tsc%c(2,2)=p1k0
* mline -- res=th(&1,&2),abcd=tsc.,m1=rmtau,m2=(-rmtau),den=0,nsum=0            
      do iut=1,2                          
      do jut=1,2
      th(iut,jut)=tsc%a(iut,jut)+rmtau*tsc%b(iut,jut)+(-rmtau)*t
     & sc%c(iut,jut)
      enddo
      enddo

      res=0.d0

      do i1=1,2
        do i2=1,2
          cres=th(i1,i2)                  
          res=res+dreal(cres)**2+dimag(cres)**2
        enddo
      enddo
      rh_tautau=res/p1k0/p2k0

      RETURN
      END
