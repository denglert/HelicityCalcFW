      SUBROUTINE vegas(region,ndim,fxn,init,ncall,itmx,nprn,tgral,sd,
     * chi2a,acc,xi,it,ndo,si,swgt,schi)
C     VEGAS - An adaptive multi-dimensional integration program
c     arguments:
c     region - input, real*8 region(2*ndim), integral region, i=1 xmin, i=2 xmax, ...  
c     ndim - input, integer, number of dimensions
c     fxn - input, external, function subprogram which computes the integrand f(x)
c     init - input, integer,
c     ncall - input, integer, the approx. number of integran evaluations per iteration
c     itmx - input, integer, the maximum number of iterations
c     nprn - input, integer, printing options
c     tgral - input, real*8,
c     sd -  input?, real*8, the standard deviation of Sbar from I
c     chi2a - output, real*8, chi2/n.o.d
c     acc - input, real*8, relative accuracy limit, algorithm stops here
c     xi - output, real*8 xi(ndmx,mxdim)
c     it - output, integer,
c     ndo - output, integer,
c     si - output, double precision,
c     swgt - output, double precision
c     schi - output, double precision,
Cc ho aggiunto acc
      INTEGER init,itmx,ncall,ndim,nprn,ndmx,mxdim,ncall_eff
      REAL*8 tgral,chi2a,sd,acc,region(2*ndim),fxn,alph,tiny
      PARAMETER (alph=1.5,ndmx=50,mxdim=10,tiny=1.e-30)
C      PARAMETER (ALPH=1.5,NDMX=100,MXDIM=10,TINY=1.e-30)
      EXTERNAL fxn
CU    USES fxn,ran2,rebin
      INTEGER i,idum,it,j,k,mds,nd,ndo,ng,npg,ia(mxdim),kg(mxdim)
      REAL*8 calls,dv2g,dxg,f,f2,f2b,fb,rc,ti,tsi,wgt,xjac,xn,xnd,xo,
     *d(ndmx,mxdim),di(ndmx,mxdim),dt(mxdim),dx(mxdim),r(ndmx),x(mxdim),
     *xi(ndmx,mxdim),xin(ndmx),ran2
      double precision schi,si,swgt,resl,standdevl
      COMMON/abresl/resl(10),standdevl(10)
      COMMON /abrann/ idum
      COMMON/abchia/calls
      COMMON/abstat/ncall_eff
      COMMON/abfla2/irepeat,nevent,nflevts
      DATA mds/1/

      SAVE
      IF(init.LE.0)THEN
        mds=1
        ndo=1
C**
        it=1
C**
        DO 11 j=1,ndim
          xi(1,j)=1.
   11   CONTINUE
      ENDIF
      IF (init.LE.1)THEN
        si=0.
        swgt=0.
        schi=0.
C**
        it=1
C**
      ENDIF
      print *,'Inside vegas.f'
      IF (init.LE.2)THEN
        nd=ndmx
        ng=1
        IF(mds.NE.0)THEN
          ng=(ncall/2.+0.25)**(1./ndim)
          mds=1
          IF((2*ng-ndmx).ge.0)then
            mds=-1
            npg=ng/ndmx+1
            nd=ng/npg
            ng=npg*nd
          ENDIF
        ENDIF
        k=ng**ndim
        npg=max(ncall/k,2)
        calls=npg*k
        dxg=1./ng
        dv2g=(calls*dxg**ndim)**2/npg/npg/(npg-1.)
        xnd=nd
        dxg=dxg*xnd
        xjac=1./calls
        DO 12 j=1,ndim
          dx(j)=region(j+ndim)-region(j)
          xjac=xjac*dx(j)
   12   CONTINUE
        IF(nd.NE.ndo)THEN
          DO 13 i=1,nd
            r(i)=1.
   13     CONTINUE
          DO 14 j=1,ndim
            CALL rebin(ndo/xnd,nd,r,xin,xi(1,j))
   14     CONTINUE
          ndo=nd
        ENDIF
        IF(nprn.GE.0) WRITE(*,200) ndim,calls,it,itmx
      ENDIF
      DO 28 it=it,itmx
C**
        IF(it.GE.2.AND.acc*abs(tgral).ge.sd) RETURN
C**
        ti=0.
        tsi=0.
        DO 16 j=1,ndim
          kg(j)=1
          DO 15 i=1,nd
            d(i,j)=0.
            di(i,j)=0.
   15     CONTINUE
   16   CONTINUE
   10   CONTINUE


c aggiunta per far scegliere a caso cella in cui valutare funzione
c  per generazione piatta con numero determinato di eventi nflevts
        IF (irepeat.EQ.2) THEN
          DO WHILE (nevent.LT.nflevts)
            wgt=xjac
            DO  j=1,ndim

              kg(j)=ran2(idum)*ng+1

              xn=(kg(j)-ran2(idum))*dxg+1.
              ia(j)=max(min(int(xn),ndmx),1)
              IF(ia(j).GT.1)THEN
                xo=xi(ia(j),j)-xi(ia(j)-1,j)
                rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
              ELSE
                xo=xi(ia(j),j)
                rc=(xn-ia(j))*xo
              ENDIF
              x(j)=region(j)+rc*dx(j)
              wgt=wgt*xo*xnd
            END DO
            f=wgt*fxn(x,wgt)
          END DO

          RETURN
        ENDIF



        fb=0.
        f2b=0.
        DO 19 k=1,npg
          wgt=xjac
          DO 17 j=1,ndim
            xn=(kg(j)-ran2(idum))*dxg+1.
            ia(j)=max(min(int(xn),ndmx),1)
            IF(ia(j).GT.1)THEN
              xo=xi(ia(j),j)-xi(ia(j)-1,j)
              rc=xi(ia(j)-1,j)+(xn-ia(j))*xo
            ELSE
              xo=xi(ia(j),j)
              rc=(xn-ia(j))*xo
            ENDIF
            x(j)=region(j)+rc*dx(j)
            wgt=wgt*xo*xnd
   17     CONTINUE
          f=wgt*fxn(x,wgt)
          f2=f*f
          fb=fb+f
          f2b=f2b+f2
          DO 18 j=1,ndim
            di(ia(j),j)=di(ia(j),j)+f
            IF(mds.GE.0) d(ia(j),j)=d(ia(j),j)+f2
   18     CONTINUE
   19   CONTINUE
        f2b=sqrt(f2b*npg)
        f2b=(f2b-fb)*(f2b+fb)
        IF (f2b.LE.0.) f2b=tiny
        ti=ti+fb
        tsi=tsi+f2b
        IF(mds.LT.0)THEN
          DO 21 j=1,ndim
            d(ia(j),j)=d(ia(j),j)+f2b
   21     CONTINUE
        ENDIF
        DO 22 k=ndim,1,-1
          kg(k)=mod(kg(k),ng)+1
          IF(kg(k).NE.1) GOTO 10
   22   CONTINUE
        tsi=tsi*dv2g
        wgt=1./tsi
        si=si+dble(wgt)*dble(ti)
        schi=schi+dble(wgt)*dble(ti)**2
        swgt=swgt+dble(wgt)
        tgral=si/swgt
        chi2a=max((schi-si*tgral)/(it-.99d0),0.d0)
        sd=sqrt(1./swgt)
        tsi=sqrt(tsi)
        IF(nprn.GE.0)THEN
C**  aggiunta di ncall_eff e sua inizializzazione dopo ogni iterazione
C**  ho modificato anche il FORMAT 201
          WRITE(*,201) it,ncall_eff,it,ti,tsi,tgral,sd,chi2a
          ncall_eff=0
          resl(it)=ti
          standdevl(it)=tsi
c          IF(nprn.NE.0)THEN
c            DO 23 j=1,ndim
c              WRITE(*,202) j,(xi(i,j),di(i,j),i=1+nprn/2,nd,nprn)
c   23       CONTINUE
c          ENDIF
        ENDIF
c** era qui l'aggiunta
        DO 25 j=1,ndim
          xo=d(1,j)
          xn=d(2,j)
          d(1,j)=(xo+xn)/2.
          dt(j)=d(1,j)
          DO 24 i=2,nd-1
            rc=xo+xn
            xo=xn
            xn=d(i+1,j)
            d(i,j)=(rc+xn)/3.
            dt(j)=dt(j)+d(i,j)
   24     CONTINUE
          d(nd,j)=(xo+xn)/2.
          dt(j)=dt(j)+d(nd,j)
   25   CONTINUE
        DO 27 j=1,ndim
          rc=0.
          DO 26 i=1,nd
            IF(d(i,j).lt.tiny) d(i,j)=tiny
            r(i)=((1.-d(i,j)/dt(j))/(log(dt(j))-log(d(i,j))))**alph
            rc=rc+r(i)
   26     CONTINUE
          CALL rebin(rc/xnd,nd,r,xin,xi(1,j))
   27   CONTINUE
   28 CONTINUE
      RETURN
  200 FORMAT(/' input parameters for vegas:  ndim=',i3,'  ncall=',
     *f12.0/28x,'  it=',i5,'  itmx=',i5)
  201 FORMAT(/' iteration no.',i3,':',12x,'effective ncall=',i11/
     *' iteration no.',i3,': ','integral =',g14.7,'+/- ',g9.2/
     *' all iterations:   integral =',g14.7,'+/- ',g9.3,' chi**2/it'
     *'n =',g9.2)
  202 FORMAT(/' data for axis ',i2/'    X       delta i       ',
     *'   x       delta i       ','    x       delta i       ',/(1x,
     *f7.5,1x,g11.4,5x,f7.5,1x,g11.4,5x,f7.5,1x,g11.4))
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,1')5c).

      FUNCTION ran2(idum)
      INTEGER idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      REAL*8 ran2,am,eps,rnmx
      PARAMETER (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     *ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
      INTEGER idum2,j,k,iv(ntab),iy
      COMMON/absalv/iv,iy,idum2
      DATA idum2/123456789/, iv/ntab*0/, iy/0/
      IF (idum.LE.0) THEN
        idum=max(-idum,1)
        idum2=idum
        DO 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          IF (idum.LT.0) idum=idum+im1
          IF (j.LE.ntab) iv(j)=idum
   11   CONTINUE
        iy=iv(1)
      ENDIF
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      IF (idum.LT.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      IF (idum2.LT.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      IF(iy.LT.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,1')5c).

      SUBROUTINE rebin(rc,nd,r,xin,xi)
      INTEGER nd
      REAL*8 rc,r(*),xi(*),xin(*)
      INTEGER i,k
      REAL*8 dr,xn,xo
      k=0
      xn=0.
      dr=0.
      DO 11 i=1,nd-1
    1   IF(rc.GT.dr)THEN
          k=k+1
          dr=dr+r(k)
          xo=xn
          xn=xi(k)
          GOTO 1
        ENDIF
        dr=dr-rc
        xin(i)=xn-(xn-xo)*dr/r(k)
   11 CONTINUE
      DO 12 i=1,nd-1
        xi(i)=xin(i)
   12 CONTINUE
      xi(nd)=1.
      RETURN
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,1')5c).

