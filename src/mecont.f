C Copyright (C) 2009 
C Shu-Kay Ng, School of Medecine, Logan Campus, Griffith University Meadowbrook, Australia
C Kim-Anh Lê Cao, ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
C and Queensland Facility for Advanced Bioinformatics, The University of Queensland, Australia
C
C This program is free software; you can redistribute it and/or
C modify it under the terms of the GNU General Public License
C as published by the Free Software Foundation; either version 2
C of the License, or (at your option) any later version.
C
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this program; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


C fortran functions to be called by R
C independence model + binary classification with continous variables only


c KA: 6/04/2009 modified as in ngme2 [see HME/currentprog]


c ---------------- SUBROUTINE E-STEP --------------------------
cka      subroutine estep(n,nv,ng,iount,x,sw,swx,swxx,xlogl,tau,u,w)
cka     subroutine estep(n,nv,ng,iount,x,sw,swx,swc,tau,u,w)

cka removed computations of sw and swx

      subroutine econt(n,nv,ng, x,
     +  tau, w, det, inv, indclass, mean, t, loglik)

      INTEGER n,nv,ng, indclass
cka      double precision sw(ng),gum, loglik,swx(ng,nv)
      double precision gum, loglik
      double precision twopi,al(ng)
      double precision x(n,nv),tau(n,ng)
      double precision u(n,ng), w(ng,nv),sumga(ng)
      double precision mean(nv,ng),t(ng),det(ng)
      double precision inv(nv,nv,ng)
      intrinsic log, sqrt, exp

c   obtain the sufficient stat
      twopi=2.*3.141592653
      loglik = 0

c ----------LOOP on the samples without loo, here n = num-1, here we
c work directly on the learning data sets

      do 800 j=1,n
        gum=0.
	do 400 k=1,ng

          sumga(k)=0.
	  u(j,k) = 0.
          al(k)=0.
	  do 300 l=1,nv
	    do 250 m=1,nv
	      al(k)=al(k)+(x(j,l)-mean(l,k))*inv(l,m,k)*
     +	      (x(j,m)-mean(m,k))
250         continue
300       continue

	  if (al(k).gt.200.) then
c	    print *, 'al(k) gt 200, k', al(k), k
	    al(k)=0.
	    goto 400
	  endif

	  al(k)=exp(-.5*al(k))/(sqrt(det(k))*(twopi)**(nv*.5))

350       if (t(k).lt.1.E-30.or.al(k).lt.1.E-30) then
c		print *, 't(k).lt.1.E-40.or.al(k).lt.1.E-40', t(k), al(k)
		goto 400
	  endif
c	  sumga(k)=0.  
	  do 360 l=1,nv
	    sumga(k)=sumga(k)+w(k,l)*x(j,l)
360       continue
	  u(j,k)=1./(1.+exp(sumga(k)))

c KA for those who develop metastases, 
	  if (j.gt.indclass) u(j,k)=u(j,k)*exp(sumga(k))
	  gum=gum+t(k)*al(k)*u(j,k)
400     continue
	if (gum.eq.0.) then
	  do 420 k=1,ng
	    tau(j,k)=0.
420       continue
c	  print *,'zero mixture density for entity in e-step',j
cka          do 200 k=1,ng
cka	  sw(k)=0.
cka200       continue
	  goto 800
        endif

	do 550 k=1,ng
	  tau(j,k)=t(k)*al(k)*u(j,k)/gum
cka	  sw(k)=sw(k)+tau(j,k)
c KA: compute loglikelihood
	  loglik = loglik + tau(j,k)* log(t(k)*al(k)*u(j,k))

cka	  do 500 l=1,nv
cka	    swx(k,l)=swx(k,l)+tau(j,k)*x(j,l)
cka500       continue
550     continue
800   continue
      return 
      end



c -------------- SUBROUTINE M-STEP ---------------------------
cka      subroutine mstep(ng,sw,swx,w)

      subroutine mcont(n,nv,ng, x, sw, swx, 
     +  tau, w, var, mean, indclass, t)

      INTEGER ng,indclass
      double precision sw(ng), sum, w(ng,nv),swx(ng,nv)
      double precision x(n,nv)
      double precision mean(nv,ng),t(ng)
      double precision var(nv,nv,ng)
      integer n,nv,gp
      double precision tau(n,ng)
c external function:
      external evalfwcont, HYBRD1cont
      integer sizehy,ifail, iflag
c KA : be careful to allocate right dimension for work()
      double precision tem(nv),xtol, taugp(n)
      double precision fvec(nv)
      intrinsic exp
      double precision work(nv*(3*nv +13)/2 + 10)

      ifail = 0

      sum=0.
      do 100 k=1,ng
	sum=sum+sw(k)
100   continue

c     Update t(k) and the mean and variance
      do 500 k=1,ng
	t(k)=sw(k)/sum
	do 400 l=1,nv
	  mean(l,k)=swx(k,l)/sw(k)
	  do 300 m=1,l
	    var(l,m,k)=0.
	  do 295 j=1,n 
	    var(l,m,k)=var(l,m,k)+(x(j,l)-mean(l,k))
     +	      *(x(j,m)-mean(m,k))*tau(j,k)
295	continue
300       continue
400     continue

     	do 450 l=1,nv
	do 460 m=1,l
		var(l,m,k) = var(l,m,k)/sw(k)
		var(m,l,k)=var(l,m,k)
460     continue
450   	continue
500   continue
cc      call matinv(nv,ng,var,inv,det)   to inverse later !!

c  Update w
      do 600 gp=1,ng
	do 510 l=1,nv
	  tem(l)=w(gp,l)
510     continue

	do 520 j=1,n
	  taugp(j)=tau(j,gp)
520     continue

c ---- call HYBRD1
	xtol=0.000000001
cka	sizehy=1000
	sizehy=nv*(3*nv+13)/2 + 10
cka	call HYBRD1(evalfwcont,nv,tem,fvec,xtol,ifail,work,sizehy)
	call  HYBRD1cont(nv,tem,fvec,xtol,ifail,work,sizehy,
     +  n, x, taugp, indclass, nv)
c	if (ifail.ne.1) print *, 'ifail= ',ifail,iount,gp
cka	print *, (fvec(l),l=1,nv)
	do 530 l=1,nv
	  w(gp,l)=tem(l)
530     continue
600   continue
      return 
      end

c-----------SUBROUTINE EVALFW --------------------------
cka      subroutine evalfwcont(ne,xw,fvecw,iflag)
      subroutine evalfwcont(nv, xw, fvec, iflag,
     + n, x, taugp, indclass, nvx)

      INTEGER nv,iflag, indclass, n, nvx
      double precision xw(nv),fvec(nv),tmp(nvx),sumxi,temp
      double precision x(n,nvx),taugp(n)

      intrinsic exp

      do 200 l=1,nv
	tmp(l)=xw(l)
	fvec(l)=0.
200   continue
      do 500 j=1,n
        sumxi=0.
	do 400 l=1,nv
          sumxi=sumxi+tmp(l)*x(j,l)
400     continue
	temp=exp(sumxi)
	sumxi=temp/(1.+temp)
	temp=0.
	if (j.gt.indclass) temp=1.
	do 450 l=1,nv
          fvec(l)=fvec(l)+taugp(j)*(temp-sumxi)*x(j,l)
450     continue
500   continue
      return
      end
      
c ----------------- SUBROUTINE CHECK ----------------------------
cka      subroutine check(n,nv,ng,ind,x,w)
cka  removed pred

      subroutine checkcont(n,nv,ng, ind,xj,
     +  w, mean,inv,det,t,
     +  gum)

      INTEGER n,nv,ng,ind,out
      double precision gum(2),twopi,al(ng)
      double precision u(n,ng),w(ng,nv),sumga(ng)
      double precision mean(nv, ng),t(ng),det(ng)
      double precision inv(nv,nv, ng)
      double precision xj(n,nv)

      intrinsic log,sqrt

      twopi=2.*3.141592653
      j=ind

c here deals only with 2 classes
      do 800 jj=1,2
        out=jj-1
        gum(jj)=0.
	do 400 k=1,ng
c should be xj
          al(k)=0.
	  do 300 l=1,nv
	    do 250 m=1,nv
	      al(k)=al(k)+(xj(j,l)-mean(l,k))*inv(l,m,k)*
     +	      (xj(j,m)-mean(m,k))
250         continue
300       continue
	  if (al(k).gt.200.) then
c	    print *, 'al(k) gt 200, k in check', al(k), k
	    al(k)=0.
	    goto 400
	  endif
	  al(k)=exp(-.5*al(k))/(sqrt(det(k))*(twopi)**(nv*.5))
350       if (t(k).lt.1.E-30.or.al(k).lt.1.E-30) then
c	        print *, 't(k)ltE-40.or.al(k)ltE-40 in check',t(k),al(k)		
		goto 400
	  endif
c ka: should be xj
	  sumga(k)=0.  
	  do 360 l=1,nv
	    sumga(k)=sumga(k)+w(k,l)*xj(j,l)
360       continue
	  u(j,k)=1./(1.+exp(sumga(k)))
	  if (out.eq.1) u(j,k)=u(j,k)*exp(sumga(k))
	  gum(jj)=gum(jj)+t(k)*al(k)*u(j,k)
400     continue
	if (gum(jj).eq.0.) then
	  print *,'zero mixture density for entity ',j
	  goto 800
        endif
800   continue
cka      pred=0
cka      if (gum(2).gt.gum(1)) pred=1
c      error=0
c      if (j.le.indclass.and.pred.eq.1) error=1
c      if (j.gt.indclass.and.pred.eq.0) error=1
c      write (28,*) j,gum(1),gum(2),out,error,(t(k),k=1,ng)
c      total=total+error
      return 
      end
