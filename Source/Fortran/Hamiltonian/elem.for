      function elem(m,n,k,l)
c
c compute the matrix element on h.o. states
c             + k    l
c     < m | (a )  (a)  | n >
c
c     implicit double precision (a-h,o-z)
      elem=0.
      if(n-l.ne.m-k) go to 90
      if(m.le.k-1) go to 90
      if(n.le.l-1) go to 90
      elem=1.
   10 if(l.eq.0) go to 20
      do i=1,l
        x=n-i+1
        elem=elem*sqrt(x)
      end do
   20 if(k.eq.0) go to 90
      do j=1,k
        y=n-l+j
        elem=elem*sqrt(y)
      end do
   90 continue
      return
      end
