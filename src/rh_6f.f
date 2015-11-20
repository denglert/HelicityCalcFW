* Function rh_6f
c Input 'p3', 'p4', 'p5', 'p6', 'p7' and 'p8'
* Process and momenta convention:
* H(p) -> e-(p3) vebar(p4) vmu(p5) mu+(p6) vtau(p7) vtaubar(p8)                 
      real*8 function rh_6f(p3,p4,p5,p6,p7,p8)

      implicit real*8 (a-b,d-h,o-z)
      implicit double complex (c)

      dimension p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3),p8(0:3)
      dimension p734(0:3),p568(0:3)

      dimension cres(2,2),cdec_taum(2,2),cdec_taup(2,2),ch_tautau(2,2)
      dimension cres_test(2,2)

      structure/pol/
        double complex e(0:3),ek0
      end structure
      record/pol/cw34,cw56

      structure/tc/
        double complex a(2,2),b(2,2),c(2,2),d(2,2)
      end structure
      record/tc/tw7_734,th734_568,tw568_8,t7_568,t568_8

      COMMON/masses/rmtau
      COMMON/couplings/wcl,gh_tautau
      COMMON/amplitudes/taup,taum

      PARAMETER (czero=(0.d0,0.d0),cim=(0.d0,1.d0))

*     print*,'p3(e-)', p3
*     print*,'p4(vebar)', p4
*     print*,'p5(vmu)', p5
*     print*,'p6(mu+)', p6
*     print*,'p7(vtau)', p7
*     print*,'p8(vtaubar)', p8

      do mu=0,3
      p734(mu)=p7(mu)+p3(mu)+p4(mu)
      p568(mu)=p5(mu)+p6(mu)+p8(mu)
      enddo

*pk0                                                                            
      p3k0=p3(0)-p3(1)
      p4k0=p4(0)-p4(1)
      p5k0=p5(0)-p5(1)
      p6k0=p6(0)-p6(1)
      p7k0=p7(0)-p7(1)
      p8k0=p8(0)-p8(1)
      p734k0=p734(0)-p734(1)
      p568k0=p568(0)-p568(1)

      cden34=1.d0
* quqd -- p=p3,q=p4                                                             
      quqd=p3(0)*p4(0)-p3(1)*p4(1)-p3(2)*p4(2)-p3(3)*p4(3)
      ccl=wcl/cden34
* TW10 -- qu=p3,qd=p4,v=0,a=cw34.e(0),cl=ccl,nsum=0                             
      eps_0=-p3(2)*p4(3)+p4(2)*p3(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p3k0*p4(0)+p4k0*p3(0)
      cw34.e(0)=ccl*(auxa-ceps_0)
* TW10 -- qu=p3,qd=p4,v=1,a=cw34.e(1),cl=ccl,nsum=0                             
      auxa=-quqd+p3k0*p4(1)+p4k0*p3(1)
      cw34.e(1)=ccl*(auxa-ceps_0)
* TW10 -- qu=p3,qd=p4,v=2,a=cw34.e(2),cl=ccl,nsum=0                             
      eps_0=-p3k0*p4(3)+p4k0*p3(3)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(2)+p4k0*p3(2)
      cw34.e(2)=ccl*(auxa-ceps_0)
* TW10 -- qu=p3,qd=p4,v=3,a=cw34.e(3),cl=ccl,nsum=0                             
      eps_0=p3k0*p4(2)-p4k0*p3(2)
      ceps_0=eps_0*cim
      auxa=p3k0*p4(3)+p4k0*p3(3)
      cw34.e(3)=ccl*(auxa-ceps_0)
      cw34.ek0=cw34.e(0)-cw34.e(1)

      cden56=1.d0
* quqd -- p=p5,q=p6                                                             
      quqd=p5(0)*p6(0)-p5(1)*p6(1)-p5(2)*p6(2)-p5(3)*p6(3)
      ccl=wcl/cden56
* TW10 -- qu=p5,qd=p6,v=0,a=cw56.e(0),cl=ccl,nsum=0                             
      eps_0=-p5(2)*p6(3)+p6(2)*p5(3)
      ceps_0=eps_0*cim
      auxa=-quqd+p5k0*p6(0)+p6k0*p5(0)
      cw56.e(0)=ccl*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=1,a=cw56.e(1),cl=ccl,nsum=0                             
      auxa=-quqd+p5k0*p6(1)+p6k0*p5(1)
      cw56.e(1)=ccl*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=2,a=cw56.e(2),cl=ccl,nsum=0                             
      eps_0=-p5k0*p6(3)+p6k0*p5(3)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(2)+p6k0*p5(2)
      cw56.e(2)=ccl*(auxa-ceps_0)
* TW10 -- qu=p5,qd=p6,v=3,a=cw56.e(3),cl=ccl,nsum=0                             
      eps_0=p5k0*p6(2)-p6k0*p5(2)
      ceps_0=eps_0*cim
      auxa=p5k0*p6(3)+p6k0*p5(3)
      cw56.e(3)=ccl*(auxa-ceps_0)
      cw56.ek0=cw56.e(0)-cw56.e(1)

      cden734=1.d0
* quqd -- p=p7,q=p734                                                           
      quqd=p7(0)*p734(0)-p7(1)*p734(1)-p7(2)*p734(2)-p7(3)*p734(
     & 3)
      ccl=wcl/cden734
* TW -- qu=p7,qd=p734,v=cw34.e,a=tw7_734.a,b=tw7_734.b,c=tw7_734.c,d=tw7_73     
* 4.d,cl=ccl,nsum=0                                                             
      ceps_0=-cw34.ek0*(p7(2)*p734(3)-p734(2)*p7(3))+p7k0*(cw34.
     & e(2)*p734(3)-p734(2)*cw34.e(3))-p734k0*(cw34.e(2)*p7(3)-p
     & 7(2)*cw34.e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw34.e(3)*p7k0+p7(3)*cw34.ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw34.e(3)*p734k0+p734(3)*cw34.ek0
      ceps_2=ceps_2*cim
      cvqu=cw34.e(0)*p7(0)-cw34.e(1)*p7(1)-cw34.e(2)*p7(2)-cw34.
     & e(3)*p7(3)
      cvqd=cw34.e(0)*p734(0)-cw34.e(1)*p734(1)-cw34.e(2)*p734(2)
     & -cw34.e(3)*p734(3)
      cauxa=-cw34.ek0*quqd+p7k0*cvqd+p734k0*cvqu
      cauxb=-cw34.ek0*p734(2)+p734k0*cw34.e(2)
      cauxc=+cw34.ek0*p7(2)-p7k0*cw34.e(2)
      tw7_734.a(2,2)=ccl*(cauxa-ceps_0)
      tw7_734.b(1,2)=ccl*(cauxb-ceps_2)
      tw7_734.c(2,1)=ccl*(-cauxc+ceps_1)
      tw7_734.d(1,1)=ccl*cw34.ek0

      cden568=1.d0/gh_tautau
      ccr=1.d0/cden568
      ccl=1.d0/cden568
* TSC -- qu=p734,qd=p568,a=th734_568.a,b=th734_568.b,c=th734_568.c,cr=ccr,c     
* l=ccl                                                                         
      auxa=-p734k0*p568(2)+p568k0*p734(2)
      cauxa=-cim*(p568(3)*p734k0-p734(3)*p568k0)
      th734_568.a(1,2)=ccl*(auxa+cauxa)
      th734_568.a(2,1)=ccr*(-auxa+cauxa)
      th734_568.b(1,1)=ccr*p568k0
      th734_568.b(2,2)=ccl*p568k0
      th734_568.c(1,1)=ccl*p734k0
      th734_568.c(2,2)=ccr*p734k0
* quqd -- p=p568,q=p8                                                           
      quqd=p568(0)*p8(0)-p568(1)*p8(1)-p568(2)*p8(2)-p568(3)*p8(
     & 3)
* TW -- qu=p568,qd=p8,v=cw56.e,a=tw568_8.a,b=tw568_8.b,c=tw568_8.c,d=tw568_     
* 8.d,cl=wcl,nsum=0                                                             
      ceps_0=-cw56.ek0*(p568(2)*p8(3)-p8(2)*p568(3))+p568k0*(cw5
     & 6.e(2)*p8(3)-p8(2)*cw56.e(3))-p8k0*(cw56.e(2)*p568(3)-p56
     & 8(2)*cw56.e(3))
      ceps_0=ceps_0*cim
      ceps_1=-cw56.e(3)*p568k0+p568(3)*cw56.ek0
      ceps_1=ceps_1*cim
      ceps_2=-cw56.e(3)*p8k0+p8(3)*cw56.ek0
      ceps_2=ceps_2*cim
      cvqu=cw56.e(0)*p568(0)-cw56.e(1)*p568(1)-cw56.e(2)*p568(2)
     & -cw56.e(3)*p568(3)
      cvqd=cw56.e(0)*p8(0)-cw56.e(1)*p8(1)-cw56.e(2)*p8(2)-cw56.
     & e(3)*p8(3)
      cauxa=-cw56.ek0*quqd+p568k0*cvqd+p8k0*cvqu
      cauxb=-cw56.ek0*p8(2)+p8k0*cw56.e(2)
      cauxc=+cw56.ek0*p568(2)-p568k0*cw56.e(2)
      tw568_8.a(2,2)=wcl*(cauxa-ceps_0)
      tw568_8.b(1,2)=wcl*(cauxb-ceps_2)
      tw568_8.c(2,1)=wcl*(-cauxc+ceps_1)
      tw568_8.d(1,1)=wcl*cw56.ek0
* TWTSC -- aa=t7_568.a,bb=t7_568.b,cc=t7_568.c,dd=t7_568.d,a1=tw7_734.a,b1=     
* tw7_734.b,c1=tw7_734.c,d1=tw7_734.d,a2=th734_568.a,b2=th734_568.b,c2=th734    * _568.c,prq=p734q,m=rmtau,nsum=0                                               
      t7_568.b(1,2)=rmtau*(tw7_734.d(1,1)*th734_568.a(1,2)+tw7_7
     & 34.b(1,2)*th734_568.b(2,2))/p734k0
      t7_568.d(1,2)=tw7_734.b(1,2)*th734_568.c(2,2)/p734k0
      t7_568.b(1,1)=(tw7_734.d(1,1)*p734q*th734_568.b(1,1)+tw7_73
     & 4.b(1,2)*th734_568.a(2,1))/p734k0
      t7_568.d(1,1)=rmtau*tw7_734.d(1,1)*th734_568.c(1,1)/p734k0
      t7_568.a(2,2)=rmtau*(tw7_734.c(2,1)*th734_568.a(1,2)+tw7_7
     & 34.a(2,2)*th734_568.b(2,2))/p734k0
      t7_568.c(2,2)=tw7_734.a(2,2)*th734_568.c(2,2)/p734k0
      t7_568.a(2,1)=(tw7_734.c(2,1)*p734q*th734_568.b(1,1)+tw7_73
     & 4.a(2,2)*th734_568.a(2,1))/p734k0
      t7_568.c(2,1)=rmtau*tw7_734.c(2,1)*th734_568.c(1,1)/p734k0
* TTW -- aa=t568_8.a,bb=t568_8.b,cc=t568_8.c,dd=t568_8.d,a1=t7_568.a,b1=t7_     
* 568.b,c1=t7_568.c,d1=t7_568.d,a2=tw568_8.a,b2=tw568_8.b,c2=tw568_8.c,d2=tw    
* 568_8.d,prq=p568q,m=rmtau,nsum=0                                              
      t568_8.c(1,1)=rmtau*(t7_568.a(1,1)*tw568_8.d(1,1)+t7_568.c
     & (1,2)*tw568_8.c(2,1))/p568k0
      t568_8.d(1,1)=(t7_568.d(1,1)*p568q*tw568_8.d(1,1)+t7_568.b(
     & 1,2)*tw568_8.c(2,1))/p568k0
      t568_8.a(1,2)=rmtau*(t7_568.a(1,1)*tw568_8.b(1,2)+t7_568.c
     & (1,2)*tw568_8.a(2,2))/p568k0
      t568_8.b(1,2)=(t7_568.d(1,1)*p568q*tw568_8.b(1,2)+t7_568.b(
     & 1,2)*tw568_8.a(2,2))/p568k0
      t568_8.c(2,1)=(t7_568.c(2,1)*p568q*tw568_8.d(1,1)+t7_568.a(
     & 2,2)*tw568_8.c(2,1))/p568k0
      t568_8.d(2,1)=rmtau*(t7_568.b(2,1)*tw568_8.d(1,1)+t7_568.d
     & (2,2)*tw568_8.c(2,1))/p568k0
      t568_8.a(2,2)=(t7_568.c(2,1)*p568q*tw568_8.b(1,2)+t7_568.a(
     & 2,2)*tw568_8.a(2,2))/p568k0
      t568_8.b(2,2)=rmtau*(t7_568.b(2,1)*tw568_8.b(1,2)+t7_568.d
     & (2,2)*tw568_8.a(2,2))/p568k0
* mline -- res=cres(&1,&2),abcd=t568_8.,m1=0,m2=0,den=0,nsum=0                  
      do iut=1,2
      do jut=1,2
      cres(iut,jut)=t568_8.a(iut,jut)
      enddo
      enddo

* mline -- res=cdec_taum(&1,&2),abcd=tw7_734.,m1=0,m2=(-rmtau),den=0,           
* nsum=0                                                                        
      do iut=1,2
      do jut=1,2
      cdec_taum(iut,jut)=tw7_734.a(iut,jut)-rmtau*tw7_734.c(iut,jut)
      enddo
      enddo

* mline -- res=cdec_taup(&1,&2),abcd=tw568_8.,m1=rmtau,m2=0,den=0,nsum=0        
      do iut=1,2
      do jut=1,2
      cdec_taup(iut,jut)=tw568_8.a(iut,jut)+rmtau*tw568_8.b(iut,jut)
      enddo
      enddo

* mline -- res=ch_tautau(&1,&2),abcd=th743_568.,m1=rmtau,m2=(-rmtau),den=0,     
* nsum=0                                                                        
      do iut=1,2
      do jut=1,2
      ch_tautau(iut,jut)=th734_568.a(iut,jut)+rmtau*th734_568.b(iut,jut)
     &                   -rmtau*th734_568.c(iut,jut)
      enddo
      enddo


      res=0.d0
      do i=1,2
      do j=1,2
      res=res+cres(i,j)*conjg(cres(i,j))
      enddo
      enddo
      res=res/p3k0/p4k0/p5k0/p6k0/p7k0/p8k0
      rh_6f=res

* H->tau+tau- amplitude
      res_htautau=0.d0
      do i=1,2
      do j=1,2
      res_htautau=res_htautau+ch_tautau(i,j)*conjg(ch_tautau(i,j))
      enddo
      enddo
      res_htautau=res_htautau/p734k0/p568k0

c     do i=1,2
c     do j=1,2
c     print*,'i ', i, 'j', j
c     print*,'cdec_taup(i.j)', cdec_taup(i,j)
c     print*,'||^2', conjg(cdec_taup(i,j))*cdec_taup(i,j)
c     enddo
c     enddo

c     print*,'' 

* tau- amplitude
      res_taum=0.d0
      do i=1,2
      do j=1,2
      res_taum=res_taum+cdec_taum(i,j)*conjg(cdec_taum(i,j))
      enddo
      enddo
      res_taum=res_taum/p3k0/p4k0/p7k0/p734k0

* tau+ amplitude
      res_taup=0.d0
      do i=1,2
      do j=1,2
      res_taup=res_taup+cdec_taup(i,j)*conjg(cdec_taup(i,j))     
      enddo
      enddo
      res_taup=res_taup/p5k0/p6k0/p8k0/p568k0

      res_nwa=res_htautau*res_taum*res_taup
      print*,'res', rea
      print*,'res_nwa', res_nwa
      print*,'res_htautau', res_htautau
      print*,'res_taum', res_taum
      print*,'res_taup', res_taup
      write(*,"(A, F10.5)") ' res_taum/res_taup', res_taum/res_taup
* Ratio between the complete chain (with full spin correlation) and
* narrow width approximation (nwa)
* Both production and decay are summed over all polarizations
c      print*,'ratio = res/res_nwa',res/res_nwa
c      print*,''

*test                                                                           

* Ratio between the amplitude of the complete result and the same result
* computed via the mulitplication of mlines in 4 sets of polarizations
      do i7=1,2
      do i8=1,2
      cres_test(i7,i8)=czero
      do i734=1,2
      do i568=1,2
      cres_test(i7,i8)=cres_test(i7,i8)+
     &                 cdec_taum(i7,i734)*ch_tautau(i734,i568)*
     &                 cdec_taup(i568,i8)/p734k0/p568k0


      enddo
      enddo
c      print*,'i7i8=',i7,i8
c      print*,'cres_test(i7,i8)',cres_test(i7,i8)
c      print*,'cres(i7,i8)',cres(i7,i8)
c      print*,'ratio_i7j8 = cres_test/cres',cres_test(i7,i8)/cres(i7,i8)
      enddo
      enddo

      RETURN
      END
