! Program of solving the point reactor neutron kinetics equations.
! First order Taylor polynomials in integral of one step is used.
! author : ikheu

program main
implicit none

integer::n						      ! number of delayed neutron
real*8,allocatable::lambda(:)	   ! decay constant
real*8,allocatable::beta(:)		! effective delayed neutron fraction
real*8::q						      ! neutron source
real*8::betaAll					   ! effective delayed neutron total fraction
real*8::age						      ! prompt neutron generation time
real*8::reacn,reaco				   ! reactivity
real*8,external::reac			   ! function of reactivity
real*8::nt						      ! neutron flux
real*8,allocatable::Cit(:)		   ! delayed neutron precursors
real*8::nt1						      ! derivative of neutron flux with time
real*8::dt						      ! time step
real*8::totaltime				      ! total problem time
real*8::time					      ! current time
real*8::editt					      ! edit time
real*8,allocatable::G1i(:)		   !-------------
real*8,allocatable::G2i(:)		   !
real*8::F1						      !
real*8::F2						      ! intermediate variable
real*8::tave					      !	
real*8::tn						      !
real*8::m1,m2,m3,m4,m5,m6,m7,m8,m9
integer::j1,j2,j3				      !--------------

! --- step1: parameters setting ---
n=6
allocate(lambda(n),beta(n),Cit(n),G1i(n),G2i(n))
lambda=(/0.0127d0,0.0317d0,0.115d0,0.311d0,1.40d0,3.87d0/)
beta=(/0.000266d0,0.001491d0,0.001316d0,0.002849d0,0.000896d0,0.000182d0/)
!betaAll=0.007
age=2.0d-5
q=0.0d0
nt=1.0d0
time=0.0d0
dt=0.01d0
totaltime=10.0d0
editt=1.0d0
betaAll=sum(beta)
Cit=nt*beta/age/lambda
j1=nint(editt/dt)
j2=0
tn=time
open(10,file="nt.txt")
open(20,file="Cit.txt")
write(10,*)"  time/s     n(t)/cm^-3"

! --- step2: begin to caculate ---
100 continue
reaco=reac(time)
time=time+dt
j2=j2+1
reacn=reac(time)
tave=time-dt/2.0
G1i=exp(-lambda*dt)*(-1+exp(lambda*dt))/lambda
G2i=exp(-lambda*dt)*(-1)*(-1+exp(lambda*dt)-lambda*dt)/lambda**2
!integral term
F1=0.5d0*(reacn+reaco)/age*dt
F2=-0.5d0*F1*dt            !F2=-0.25d0*(reacn+reaco)/age*dt**2
! intermediate variable
m1=sum(Cit)
m2=sum(exp(-lambda*dt)*Cit)
m3=sum(beta*G2i)/age
m4=sum(lambda*exp(-lambda*dt)*Cit)
m5=1-sum(lambda*beta*G2i)/age
m6=sum(lambda*beta*G1i)/age
m7=sum(beta*G1i)/age
! target variable
nt=(nt+q*dt+m1-m2+(F2-m3)*(m4+q)/m5)/(1-F1+m7-(F2-m3)*((reacn-betaAll)/age+m6)/m5)
nt1=(((reacn-betaAll)/age+m6)*nt+m4+q)/m5
Cit=exp(-lambda*dt)*Cit+beta/age*(G1i*nt+G2i*nt1)

! --- step3: output ---
if(mod(j2,j1)==0)then
   write(10,"(f10.6,es20.6)")time,nt
end if
tn=time
if (time < totaltime) goto 100

stop
end