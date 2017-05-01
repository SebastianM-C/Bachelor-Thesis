c compute matrix elements of a boson hamiltonian
c and diagonalise it
c
c     implicit double precision (a-h,o-z)
      parameter(nd1=2000000000,nd2=100000)
      dimension aa(nd1), xx(nd1),ee(nd2)
      dimension index(nd2)
      common/param/a,b,d
      common/pinteg/m1,m2,m3,iii
      common/ener/dd(nd2)
c      real(kind=8), allocatable, dimension(:) :: aa
c      integer :: reclen
c      allocate(aa(nd1))
c      inquire(iolength=reclen) aa
 1000 format(/' nn =',i5
     *       /' i  =',i5
     *       /' j  =',i5)
 2000 format(/30x,'Matrix elements'/)
 2001 format(/32x,'  Eigenvalues'/)
 2002 format(/25x,'Eigenvector =',i8,'  E =',f12.5/)
 2003 format(i6,f12.5)
 2004 format(/30x,'   Structure'/)
 2010 format(10f8.3)
 2020 format(10i8)
 3000 format(' *** k =',i10)
 3001 format(' i j i1 j1 i2 j2 a =',6i4,f8.3)
      open(unit=10,file='hamilt.inp',status='old')
      open(unit=20,file='hamilt.bin',status='new',form='unformatted')
c     *     access='direct',action='write',recl=reclen)
      open(unit=30,file='hamilt.dat',status='new')
      open(unit=40,file='rebde.dat',status='new')
      open(unit=50,file='reuns.dat',status='new')
      open(unit=60,file='reuna.dat',status='new')
      open(unit=70,file='eigenvectors.out',status='new')
      open(unit=80,file='index.out',status='new')

      read(10,*) n,a,b,d
c     nn=n*n
      nn=n*(n+1)/2
      i=0
      do i1=1,n
        do i2=1,n-i1+1
          i=i+1
          index(i)=i1*100+i2
          write(80, *) i1 - 1, " ", i2 - 1
c         i=n*(i1-1)+i2
          j=0
          do j1=1,n
            do j2=1,n-j1+1
              j=j+1
c             j=n*(j1-1)+j2
              k=nn*(i-1)+j
              aa(k)=ham(i1-1,i2-1,j1-1,j2-1)
c              if(k.gt.nd1) write(*,3000)k
c             if(aa(k).ne.0)
c    *        write(20,3001)i,j,i1,j1,i2,j2,aa(k)
            end do
          end do
        end do
      end do
c     nn=i
      write(* ,1000)nn,i,j
c199      write(20,1000)nn,i,j
c     write(20,2000)
c     call wmat(20,nn,1,nn,1,aa)
c      is=1
c      call sdiag(nn,aa,dd,xx,ee,is)
c199      write(20,2001)
c199      write(20,2010)(dd(i),i=1,nn)
c199      write(20,2004)
c199      write(20,2020)(index(i),i=1,nn)
c      write(*, *) char(10), " index:"
c      write(*, *) (index(i), i = 1, nn)
c      write(80, *) (index(i), i = 1, nn)
C      write(*, *) "H transposed:", char(10)
c      do i = 1, nn
c        do j = 1, nn
c          write(20, '(f12.8,$)') aa((i - 1) * nn + j)
c        end do
c        write(20, *) char(10)
c      end do
c      write(*, *) "eigenvectors: ", char(10)
c      do i = 1, nn
c        write(*, '(a, i5)') "v", i
c        do j = 1, nn
c          write(*, '(f12.8,$)') xx((i - 1) * nn + j)
c        end do
c        write(*, *) char(10)
c      end do
c      write(*,*) "eigenvectors: ", char(10), (xx(i),i=1,nn*nn)
      do j=1,nn
c199        write(20,2002)j,dd(j)
c        do k=1, nn
c          write(*, *) ((j-1)*nn+k), k, xx((j-1)*nn+k), dd(j)
c        enddo
c        write(70, *) (xx((j - 1) * nn + k), k = 1, nn)
        write(20) (aa((j - 1) * nn + k), k = 1, nn)
c        write(30,2003)j,dd(j)


c
c        if(j.eq.1) then
c          m2=m2+1
c          write(50,2003) m2,dd(j)
c          go to 99
c          endif
c          dif1=dd(j+1)-dd(j)
c          dif=dd(j)-dd(j-1)
c          if(dif.le.0.0005) then
c           m1=m1+1
c           write(40,2003) m1,dd(j)
c          endif
c          if(dif.gt.0.0005.and.dif1.gt.0.0005) then
c           call select(nn,j,j,2,xx)
c            if(iii.eq.1) then
c             m2=m2+1
c             write(50,2003) m2,dd(j)
c            endif
c            if(iii.eq.2) then
c             m3=m3+1
c             write(60,2003) m3,dd(j)
c            endif
c          endif
c 99        continue
c        call wmat(20,nn,j,j,2,xx)
      end do
c      deallocate(aa)
      close(unit=10)
      close(unit=20)
      close(unit=30)
      close(unit=40)
      close(unit=50)
      close(unit=60)
      stop
      end

      function ham(m1,m2,n1,n2)
c
c compute matrix elements of the Hamiltonian
c
c     implicit double precision (a-h,o-z)
      common/param/a,b,d
      ham=     a*(elem(m1,n1,1,1)*elem(m2,n2,0,0)
     +           +elem(m1,n1,0,0)*elem(m2,n2,1,1))
      ham=ham
     + +0.25*b*(3*elem(m1,n1,1,0)*elem(m2,n2,2,0)
     +         +3*elem(m1,n1,0,1)*elem(m2,n2,0,2)
     -           -elem(m1,n1,3,0)*elem(m2,n2,0,0)
     -           -elem(m1,n1,0,3)*elem(m2,n2,0,0))
     +   +0.75*b*(elem(m1,n1,0,1)*elem(m2,n2,2,0)
     +           +elem(m1,n1,1,0)*elem(m2,n2,0,2)
     -           -elem(m1,n1,1,2)*elem(m2,n2,0,0)
     -           -elem(m1,n1,2,1)*elem(m2,n2,0,0)
     +         +2*elem(m1,n1,0,1)*elem(m2,n2,1,1)
     +         +2*elem(m1,n1,1,0)*elem(m2,n2,1,1))
      ham=ham
     +  +0.375*d*(elem(m1,n1,2,2)*elem(m2,n2,0,0)
     +           +elem(m1,n1,0,0)*elem(m2,n2,2,2))
     +  +0.125*d*(elem(m1,n1,2,0)*elem(m2,n2,0,2)
     +           +elem(m1,n1,0,2)*elem(m2,n2,2,0))
     +  +0.500*d* elem(m1,n1,1,1)*elem(m2,n2,1,1)
     +  +0.250*d*(elem(m1,n1,1,3)*elem(m2,n2,0,0)
     +           +elem(m1,n1,3,1)*elem(m2,n2,0,0)
     +           +elem(m1,n1,0,0)*elem(m2,n2,1,3)
     +           +elem(m1,n1,0,0)*elem(m2,n2,3,1)
     +           +elem(m1,n1,0,2)*elem(m2,n2,1,1)
     +           +elem(m1,n1,2,0)*elem(m2,n2,1,1)
     +           +elem(m1,n1,1,1)*elem(m2,n2,0,2)
     +           +elem(m1,n1,1,1)*elem(m2,n2,2,0))
     + +0.0625*d*(elem(m1,n1,4,0)*elem(m2,n2,0,0)
     +           +elem(m1,n1,0,4)*elem(m2,n2,0,0)
     +           +elem(m1,n1,0,0)*elem(m2,n2,4,0)
     +           +elem(m1,n1,0,0)*elem(m2,n2,0,4)
     +         +2*elem(m1,n1,2,0)*elem(m2,n2,2,0)
     +         +2*elem(m1,n1,0,2)*elem(m2,n2,0,2))
      return
      end

      subroutine wmat(nf,n,n1,n2,ind,a)
c
c write matrix a(n,n) in the file nf
c if nf = 0 write on display
c
c ind = 1 : write from  line  n1 to  line  n2
c       2 : write from column n1 to column n2
c
c     implicit double precision (a-h,o-z)
      dimension a(n,n)
 2000 format(10f8.3)
 4003 format(i5,f8.3)
c199      do i=n1,n2
c199        go to (10,20) ind
c199   10   if(nf.ne.0) write(20,2000)(a(i,j),j=1,n)
c199        if(nf.eq.0) write(* ,2000)(a(i,j),j=1,n)
c199        go to 30
c199   20   if(nf.ne.0) write(20,2000)(a(j,i),j=1,n)
c199        if(nf.eq.0) write(* ,2000)(a(j,i),j=1,n)
   30   continue
c199       enddo
      return
      end

      subroutine select(n,n1,n2,ind,bb)
c
c read matrix a(n,n)
c if nf = 0 write on display
c
c ind = 1 : write from  line  n1 to  line  n2
c       2 : write from column n1 to column n2
c
c     implicit double precision (a-h,o-z)
      parameter(nd1=200000000,nd2=100000)
      dimension bb(nd1)
      common/pinteg/m1,m2,m3,iii
      common/ener/dd(nd2)
 2000 format(10f8.3)
 4003 format(i5,f8.3)
       i=n1
        go to (10,20) ind
 10      continue
c   10   if(nf.ne.0) write(20,2000)(a(i,j),j=1,n)
c        if(nf.eq.0) write(* ,2000)(a(i,j),j=1,n)
c        go to 30
c   20   if(nf.ne.0) write(20,2000)(a(j,i),j=1,n)
c        if(nf.eq.0) write(* ,2000)(a(j,i),j=1,n)
  20     continue
             j1=(i-1)*n+1
             j2=i*n
           do j=j1,j2
             nr=j-j1+1
             xxx=abs(bb(j))
            if(xxx.gt.0.0001.and.nr.ne.1) then
c               write(*,*) bb(j)
               k1=(nr-1)/2
               k2=2*k1
              if(k2.eq.nr-1) then
               iii=1
              endif
              if(k2.ne.nr-1) then
               iii=2
              endif
           go to 70
          endif
       enddo
70     continue
c       enddo
      return
      end

      subroutine descom(n,i,j,i1,i2,j1,j2)
c
c find i1,i2,j1,j2 from:
c
c i = n * (i1-1) + i2
c j = n * (j1-1) + j2
c
      i1=i/n+1
      j1=j/n+1
      j1=i-n*(i1-1)
      j2=i-n*(i2-1)
      return
      end
