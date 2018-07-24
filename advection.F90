program advection1d_driver

implicit none

integer :: nx, nt
real, allocatable, dimension(:) :: x, f, u
real, allocatable, dimension(:) :: fp, up
real :: dx, dt, cour, s, tltest, tol

real, parameter :: pi = 3.14159265358979323846

integer :: j

real :: dp(2)
real, allocatable, dimension(:) :: f0, f1, f2
real, allocatable, dimension(:) :: u0, u1, u2
real, allocatable, dimension(:) :: fp11, up11
real, allocatable, dimension(:) :: fp12, up12
real, allocatable, dimension(:) :: fp21, up21
real, allocatable, dimension(:) :: fp22, up22

real :: ttest = 1.0e-2

!Number of grid points
!---------------------
nx = 100

!Allocate memory
!---------------
allocate(x(nx))
allocate(f(nx))
allocate(fp(nx))
allocate(u(nx))
allocate(up(nx))
allocate(f0(nx))
allocate(f1(nx))
allocate(f2(nx))
allocate(u0(nx))
allocate(u1(nx))
allocate(u2(nx))
allocate(fp11(nx))
allocate(up11(nx))
allocate(fp12(nx))
allocate(up12(nx))
allocate(fp21(nx))
allocate(up21(nx))
allocate(fp22(nx))
allocate(up22(nx))


!Grid
!----
s = 100.0
x = (/((j*2*pi/(nx-1)),j=0,nx-1)/)
dx = x(2) - x(1)
dt = dx/s
nt = int(nx*s)


!Initial condition
!-----------------
f = sin(x)
u = -1.0

f0 = f   !Save a copy
u0 = u


!Stability check
!---------------
cour = maxval(u*dt/dx)
if (cour > 1.0) then
  print*, 'Unstable time step, exiting.'
  return
endif


!Call nonlinear advection scheme
!-------------------------------
call advection1d(nx,nt,dx,dt,x,f,u)


!Write forecast initial and final state
!--------------------------------------
call writef(nx,f0,'q_init.txt')
call writef(nx,f ,'q_final.txt')


!Tangent test
!------------
print*, ' '
print*, 'Tangent linear test'
print*, '-------------------'

tltest = 10.0
do j = 1,10

  f1 = f
  u1 = u
  
  call random_number(f)
  call random_number(u)
  
  fp = ttest*f
  up = ttest*u
  
  f2 = f0 + fp
  u2 = u0 + up
  
  call advection1d(nx,nt,dx,dt,x,f2,u2)
  
  f = f0
  u = u0
  
  call advection1d_tl(nx,nt,dx,dt,x,f,fp,u,up)
  
  tltest = min(tltest,maxval((f2-f1)/fp))

  print*, 'Pert scaled by, tangent test result:', ttest, ', ', maxval((f2-f1)/fp)
  
  ttest = ttest/10

enddo

!Tolerance for passing/failing tangent linear test
tol = 1.0e-4

print*, ' '
print*, 'Best: ', tltest, tol
print*, ' '
if (abs(tltest-1) < tol) then
  print*, 'TL Test Passed with tolerance', tol
else
  print*, 'TL Test Failed with tolerance', tol
endif
print*, ' '


!Dot product test
!----------------

call random_number(fp11)
call random_number(up11)
call random_number(fp22)
call random_number(up22)

fp12 = fp11
up12 = up11
fp21 = fp22
up21 = up22

f = f0
u = u0
call advection1d_tl(nx,nt,dx,dt,x,f,fp12,u,up12)

f = f0
u = u0
call advection1d_ad(nx,nt,dx,dt,x,f,fp21,u,up21)

dp = 0.0
do j = 1,nx
  dp(1) = dp(1) + fp11(j)*fp21(j)
  dp(1) = dp(1) + up11(j)*up21(j)
  dp(2) = dp(2) + fp12(j)*fp22(j)
  dp(2) = dp(2) + up12(j)*up22(j)
enddo

print*, 'Dot product test'
print*, '----------------'
print*, '<x,xhat>, <y,yhat>, relative difference:', dp(1), dp(2), (dp(2) - dp(1))/dp(1)
print*, ' '


!Run adjoint with a 1 somewhere
!------------------------------




!Deallocate memory
!-----------------
deallocate(x)
deallocate(f)
deallocate(fp)
deallocate(u)
deallocate(up)
deallocate(f0)
deallocate(f1)
deallocate(f2)
deallocate(u0)
deallocate(u1)
deallocate(u2)
deallocate(fp11)
deallocate(up11)
deallocate(fp12)
deallocate(up12)
deallocate(fp21)
deallocate(up21)
deallocate(fp22)
deallocate(up22)


!------------------------------------
contains
!------------------------------------


subroutine advection1d(nx,nt,dx,dt,x,f,u)

implicit none
integer, intent(in)    :: nx, nt
real,    intent(in)    :: dx, dt, x(nx)
real,    intent(inout) :: f(nx)
real,    intent(in)    :: u(nx)

integer :: n,j
real :: fold(nx)

do n = 1,nt

  fold = f

  do j = 1,nx

     if (j == 1) then
       f(j) = fold(j) + dt*u(j)/dx * (fold(j) - fold(nx))
     else
       f(j) = fold(j) + dt*u(j)/dx * (fold(j) - fold(j-1))
     endif

  enddo

enddo

end subroutine advection1d

!------------------------------------

subroutine advection1d_tl(nx,nt,dx,dt,x,f,fd,u,ud)

implicit none
integer, intent(in)    :: nx, nt
real,    intent(in)    :: dx, dt, x(nx)
real,    intent(inout) :: f(nx), fd(nx)
real,    intent(in)    :: u(nx), ud(nx)

integer :: n,j
real :: fold(nx), foldd(nx)


end subroutine advection1d_tl

!------------------------------------

subroutine advection1d_ad(nx,nt,dx,dt,x,f,fb,u,ub)

implicit none
integer, intent(in)    :: nx, nt
real,    intent(in)    :: dx, dt, x(nx)
real,    intent(inout) :: f(nx), fb(nx)
real,    intent(inout) :: u(nx), ub(nx)

integer :: n,j,branch
real :: fold(nx), foldb(nx)
real :: tempb, tempb0


end subroutine advection1d_ad

!------------------------------------

subroutine writef(nx,f,fname)

implicit none
integer, intent(in) :: nx
real,    intent(in) :: f(nx)
character(len=*)    :: fname

integer :: unit

!Write f to file for reading by plotting software
open(unit = 101, file = trim(fname))
do j = 1,nx
  write(101,*) f(j)
enddo
close(101)

end subroutine writef

!------------------------------------

end program advection1d_driver
