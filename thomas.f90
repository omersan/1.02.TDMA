!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Thomas Algorithm to solve tridiagonal system of equations
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Aug. 18, 2015
!-----------------------------------------------------------------------------!

program tdma_test
implicit none
integer::i,n
real*8, dimension(:),allocatable::a,b,c,r,x

n = 4

allocate(a(n),b(n),c(n),r(n),x(n))
!---------------------------------------------!
!Given system
!---------------------------------------------!
do i=1,n
a(i) = 1.0d0
b(i) = 4.0d0
c(i) = 2.0d0 
end do

r(1) = 6.0d0
r(2) = 7.0d0
r(3) = 7.0d0 
r(n) = 5.0d0 

!r(1) = 8.0d0
!r(2) = 15.0d0
!r(3) = 22.0d0 
!r(n) = 19.0d0 

!Thomas Algorithm
call tdma(a,b,c,r,x,1,n)

!Print solution 
do i=1,n
write(*,*) i, x(i)
end do

       
end


!------------------------------------------------------------------!
!Tridiagonal matrix algorithm (TDMA)
!Thomas algorithm
!solution tridiagonal systems
!a: lower diagonal
!b: main diagonal
!c: upper diagonal
!r: source vector
!x: solution vector
!   for indices s(start) to e(end)
!   i: s,s+1,s+2, ....,i,....,e 
!
!Note: a(s) and c(e) are dummy coefficients, not used.
!------------------------------------------------------------------!

subroutine tdma(a,b,c,r,x,s,e)
implicit none
integer s,e,i
real*8, dimension(s:e) ::a,b,c,r,x    

! forward elimination phase
do i=s+1,e
b(i) = b(i) - a(i)/b(i-1)*c(i-1)
r(i) = r(i) - a(i)/b(i-1)*r(i-1)
end do
! backward substitution phase 
x(e) = r(e)/b(e)
do i=e-1,s,-1
x(i) = (r(i)-c(i)*x(i+1))/b(i)
end do

return
end
