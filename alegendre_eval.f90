! This code is modified from ALegendreEval by James Bremer .
! 
! This file contains code for evaluating the scaled, normalized 
! associated Legendre functions of the first kind.
! The code for the computation of the scaled, normalized 
! associated Legendre functions of the second kind has been 
! excluded!
! 
! It runs in time independent of the degree dnu and order dmu
! and can be applied when 0 <= dnu <= 1,000,000 and 
! 0 <= dmu <= dnu.  
! 
! Reference:    James Bremer, "An algorithm for the numerical 
! evaluation of the associated Legendre functions in time 
! independent of degree and order."  arXiv:1707.03287


module alegendreeval

use, intrinsic :: iso_c_binding
use besseleval

implicit double precision (a-h,o-z)

!
!  The following structure stores one of the precomputed expansions.
!
type       alegendre_expansion_data

integer                          :: ifbetas,ifover,ifsmalldmu,ifinverse
double precision                 :: dnu1,dnu2
integer                          :: ncoefs,k,nintsab,nintscd,nintsef
double precision, allocatable    :: ab(:,:),cd(:,:),ef(:,:)

!
!  Phase function and its derivative
!

integer                          :: ncoefsalpha,ncoefsalphap
double precision, allocatable    :: coefsalpha(:),coefsalphap(:)
integer, allocatable             :: iptrsalpha(:,:,:),iptrsalphap(:,:,:)

!
!  Second derivative of the phase function at the turning point
!

integer                          :: ncoefsalphapp
double precision, allocatable    :: coefsalphapp(:)
integer, allocatable             :: iptrsalphapp(:,:)

!
!  Logarithms
!

integer                          :: ncoefsbeta1,ncoefsbeta2,nintsabb
double precision, allocatable    :: coefsbeta1(:),coefsbeta2(:),abb(:,:)
integer, allocatable             :: iptrsbeta1(:,:,:),iptrsbeta2(:,:,:)

!
!  For the inverse
!

integer                          :: ncoefsalphainv,nintsinv
double precision, allocatable    :: coefsalphainv(:),abinv(:,:)
integer, allocatable             :: iptrsalphainv(:,:,:)

!
!  Some meta data regarding the precomputation
!

double precision                 :: epsrequired,epsphase,epsdisc,epscomp
double precision                 :: dmemory,dtime

end type   alegendre_expansion_data


type(alegendre_expansion_data), private :: expdata1,expdata2,expdata3,expdata4
type(alegendre_expansion_data), private :: expdata5,expdata6
type(alegendre_expansion_data), private :: expdata7,expdata8
integer, private                        :: ifloaded,maxdegree
double precision, private               :: pi

data pi            / 3.14159265358979323846264338327950288d0  /
data ifloaded      / 0         /
data maxdegree     / 100000000 /

contains


! subroutine alegendre_eval(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq)
subroutine alegendre_eval(dnu,dmu,t,vallogp,valp)
implicit double precision (a-h,o-z)

!
!  Evaluate the *scaled*, *normalized* associated Legendre functions of the first 
!  and second kinds (1) and (2) of degrees nu and orders -dmu at the point cos(t), 
!  where  0 < t <= pi/2 and (dnu,dmu) is in the set (3).
!
!  When (dnu,dmu,t) is in the oscillatory region (4), also return the values of a 
!  nonoscillatory phase function alpha for the associated Legendre differential 
!  equation and its derivative.
!
!  When (dnu,dmu,t) is in the nonoscillatory region (5), also return the values of
!  the logarithms of the functions (1) and (2).
!
!  Input parameters:
!    dnu - the degree of the associated Legendre functions to evaluate
!    dmu - the order of the associated Legendre functions to evaluate
!    t - the argument at which to evaluate them
!
!  Output parameters:
!    vallogp - the value of the logarithm of (1) if (dnu,dmu,t) is in (5)
!    valp - the value of (1)
!

! alpha     = 0
! alphader  = 0
vallogp   = 0
! vallogq   = 0
valp      = 0
! valq      = 0

!
!  Check to see that the precomputed expansions have been loaded
!  

if (ifloaded .eq. 0) then
print *,"alegendre_eval:  alegendre_eval_init must be called before alegendre_eval"
stop
endif

!
!  Perform range checking. 
!

if (dnu .lt. 0 .OR. dnu .gt. maxdegree   .OR. &
    dmu .lt. 0 .OR. dmu .gt. dnu         .OR. &
    t .le. 0                             .OR. &
    t .gt. pi) then
print *,"alegendre_eval: parameters out of range"
print *,"dnu = ",dnu
print *,"dmu = ",dmu
print *,"t   = ",t
stop
endif

dlambda        = dnu+0.5d0
ifoscillatory  = 1

if (dmu .gt. 0.5d0) then
deta    = sqrt(dmu**2-0.25d0)
a       = asin(deta/dlambda)
if (t .lt. a) ifoscillatory = 0
endif

ifsmalldmu = 0
if (dmu .lt. 1) ifsmalldmu = 1
call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)


!
!  For small dnu, always use series expansions
!

if (dnu .lt. 2) then
! call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
call alegendre_taylor(dnu,dmu,t,vallogp,valp,ifoscillatory)
return
endif

!
!  Handle the case of 2 <= dnu < 10
!


if (dnu .lt. 10) then


if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
! call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
call alegendre_taylor(dnu,dmu,t,vallogp,valp,ifoscillatory)
else
! call alegendre_expeval(expdata1,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata1,dnu,dmu,t,vallogp,valp,a,b,c)
endif

else

if (t .lt. a) then
! call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
call alegendre_taylor(dnu,dmu,t,vallogp,valp,ifoscillatory)
else
! call alegendre_expeval(expdata2,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata2,dnu,dmu,t,vallogp,valp,a,b,c)
endif

endif

return
endif


!
!  Now 10 <= dnu <= 10,000 
!

if (dnu .le. 10000) then

if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
! call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
call alegendre_taylor(dnu,dmu,t,vallogp,valp,ifoscillatory)
else
! call alegendre_expeval(expdata3,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata3,dnu,dmu,t,vallogp,valp,a,b,c)
endif

else

if (t .lt. c) then
! call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
call alegendre_taylor(dnu,dmu,t,vallogp,valp,ifoscillatory)
else
! call alegendre_expeval(expdata4,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata4,dnu,dmu,t,vallogp,valp,a,b,c)
endif
endif

return
endif


!`
!  Now 10,000 < dnu <= 1,000,000
!

if (dnu .le. 1000000) then


if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
! call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
call alegendre_macdonald(dnu,dmu,t,vallogp,valp)
else
! call alegendre_expeval(expdata5,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata5,dnu,dmu,t,vallogp,valp,a,b,c)
endif

else

if (t .lt. c) then
!call alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,ifoscillatory)
! call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
call alegendre_macdonald(dnu,dmu,t,vallogp,valp)
else
! call alegendre_expeval(expdata6,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata6,dnu,dmu,t,vallogp,valp,a,b,c)
endif
endif

return
endif



!
!  Now 1,000,000 < dnu <= 100,000,000
!


if (ifsmalldmu .eq. 1) then

if (t .lt. a) then
! call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
call alegendre_macdonald(dnu,dmu,t,vallogp,valp)
else
! call alegendre_expeval(expdata7,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata7,dnu,dmu,t,vallogp,valp,a,b,c)
endif

else

if (t .lt. c) then
! call alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
call alegendre_macdonald(dnu,dmu,t,vallogp,valp)
else
! call alegendre_expeval(expdata8,dnu,dmu,t,alpha,alphader,vallogp,vallogq,valp,valq,a,b,c)
call alegendre_expeval(expdata8,dnu,dmu,t,vallogp,valp,a,b,c)
endif
endif

return



end subroutine




subroutine alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
implicit double precision (a-h,o-z)


if (ifsmalldmu .eq. 0) then
deta    = sqrt(dmu**2-0.25d0)
dlambda = dnu+0.5d0
a       = asin(deta/dlambda) 
b       = pi/2
c       = asin(deta/dlambda)/100.0d0
else
a       = 1.0d0/dnu**1.1d0
b       = pi/2
endif

end subroutine


! subroutine alegendre_taylor(dnu,dmu,t,alpha,alphader,vallogp,vallogq, &
!   valp,valq,ifoscillatory)
subroutine alegendre_taylor(dnu,dmu,t,vallogp,valp,ifoscillatory)
implicit double precision (a-h,o-z)

!
!  Evaluate the SCALED, NORMALIZED associated Legendre functions of the first
!  and second kinds (1) and (2) of degree dnu and order -dmu at the point cos(t)
!  using series expansions.
!

!
!  In the nonoscillatory regime, compute the logarithms.
!

if (ifoscillatory .eq. 0) then
call alegendre_ptaylor_log(dnu,dmu,t,vallogp,dsignp)
! call alegendre_qtaylor_log(dnu,dmu,t,vallogq)
valp    = exp(vallogp)
! valq    = exp(vallogq)
return
endif

!
!  In the oscillatory regime, use the series expansions and compute
!  alpha' and alpha as needed.  This should only be done for small
!  values of dnu, dmu and t which are less than the first zeros
!  of P and Q.
!


call alegendre_ptaylor(dnu,dmu,t,valp)
! call alegendre_qtaylor(dnu,dmu,t,valq)

! alphader = (2*dnu+1.0d0)/pi*1/(valp**2 + valq**2)
! alpha    = atan2(-valq,valp) + 2*pi
vallogp  = 0
! vallogq  = 0

end subroutine


subroutine alegendre_ptaylor(dnu,dmu,t,val)
implicit double precision (a-h,o-z)

!
!  Evaluate the SCALED, NORMALIZED associated Legendre function (1)
!  using a Taylor expansion. 
!
!  This routine should only be used in the oscillatory regime, and then
!  only for small values of dnu and dmu.
!
!  Input parameters:
!    dnu - the degree of the *normalized* associated Legende function to evaluate
!    dmu - the order of the *normalize* associated Legendre function to evaluate
!    t - the point on the interval (0,\pi/2) at which to evaluate the function
!  
!  Output parameters:
!    val - the value of P_\nu^{-\mu} (\cos(t))
!

eps0     = epsilon(0.0d0)/10.0d0
i1       = 0
i2       = 49
if (floor(dnu) .eq. dnu) then
i2 = dnu
endif

if (dmu .lt. 0 .AND. ( floor(dmu) .eq. dmu )) then
i1 = -dmu
endif


dsum = 0 
dd   = 1.0d0/alegendre_gamma(dmu+i1+1.0d0) * 1/alegendre_gamma(i1+1.0d0) * sin(t/2)**(2*i1) * (-1)**(i1)


if (i1 .ne. 0) then
dd   = dd * alegendre_gamma(dnu+i1+1.0d0)/alegendre_gamma(dnu-i1+1.0d0)
endif

xx   = sin(t/2)**2

do i=i1,i2
dsum = dsum  + dd

dd   = -dd * (dnu+i+1.0d0)*(dnu-i) * 1.0d0/(i+1.0d0) * 1.0d0/(dmu+i+1.0d0) * xx 
if (abs(dd) .lt. eps0*abs(dsum)) exit
end do

val = tan(t/2)**dmu * dsum

!
!  Now apply the normalization and scaling factors
!

call alegendre_gamma_ratio2(dnu,dmu,dd)
val = val * sqrt( (dnu+0.5d0) * exp(dd)) * sqrt(sin(t))

end subroutine


subroutine alegendre_ptaylor_log(dnu,dmu,t,vallog,dsign)
implicit double precision (a-h,o-z)

!
!  Calculate the value of the NORMALIZED, SCALED, associated legendre
!  function of the first kind (1).
!
!  Input parameters:
!    dnu - the degree of the *normalized* associated Legende function to evaluate
!    dmu - the order of the *normalize* associated Legendre function to evaluate
!    t - the point on the interval (0,\pi/2) at which to evaluate the function
!  
!  Output parameters:
!    dsign - the sign of P_nu^{-mu}(cos(t))
!    vallog - logarithm of the absolute value of P_dnu^{-dmu} ( cos (t ))
!

eps0     = epsilon(0.0d0)/10.0d0
i1       = 0
i2       = 49

if (floor(dnu) .eq. dnu) then
i2 = dnu
endif

if (dmu .lt. 0 .AND. ( floor(dmu) .eq. dmu )) then
i1 = -dmu
endif

dsign = (-1.0d0)**i1
dd   = -log_gamma(dmu+i1+1.0d0)-log_gamma(i1+1.0d0) + 2*i1 * log(sin(t/2)) 


if (i1 .ne. 0) then
dd   = dd + log_gamma(dnu+i1+1.0d0) - log_gamma(dnu-i1+1.0d0)
endif


ddd = alegendre_gamma(dmu+i1+1.0d0)

if (ddd .lt. 0)  dsign = -dsign

dd2   = 1.0d0
dsum  = 0 
xx    = sin(t/2)**2


do i=i1,i2
dsum = dsum  + dd2
dd2   = -dd2 * (dnu+i+1.0d0)*(dnu-i) * 1.0d0/(i+1.0d0) * 1.0d0/(dmu+i+1.0d0) * xx 
if (abs(dd2) .le. eps0*abs(dsum)) exit
end do

if (dsum .lt. 0) then
dsum  = -dsum
dsign  = -dsign
endif


vallog = dd + log(dsum) + dmu * log(tan(t/2))


!
!  Now the normalization factor
!

call alegendre_gamma_ratio2(dnu,dmu,dd)


vallog = vallog + 0.5d0 * (log(dnu+0.5d0) + dd + log(sin(t)) )
val    = dsign*exp(vallog)


end subroutine


subroutine alegendre_gamma_ratio2(dnu,dmu,val)
implicit double precision (a-h,o-z)

!
!  Evaluate log ( Gamma(dnu+dmu+1) / Gamma(dnu-dmu+1) )
!


if (dnu .lt. 10 .OR. abs(dmu) .gt. dnu/1000.0d0) then
val = log_gamma(dnu+dmu+1) - log_gamma(dnu-dmu+1)
else

x = dnu+1
y = dmu

dsum = 1-(0.1d1*y)/x-(y*(1-0.3d1*y+0.2d1*y**2))/(0.6d1*x**2)+  &
(y**2*(1-0.3d1*y+0.2d1*y**2))/(0.6d1*x**3)+(y*(6+0.5d1*y-  &
0.9d2*y**2+0.155d3*y**3-0.96d2*y**4+  &
0.2d2*y**5))/(0.36d3*x**4)-(y**2*(6+0.5d1*y-0.9d2*y**2+  &
0.155d3*y**3-0.96d2*y**4+0.2d2*y**5))/(0.36d3*x**5)-(y*(360+  &
0.126d3*y-0.2863d4*y**2-0.1323d4*y**3+0.14385d5*y**4-  &
0.18711d5*y**5+0.10518d5*y**6-0.2772d4*y**7+  &
0.28d3*y**8))/(0.4536d5*x**6)+(y**2*(360+0.126d3*y-  &
0.2863d4*y**2-0.1323d4*y**3+0.14385d5*y**4-0.18711d5*y**5+  &
0.10518d5*y**6-0.2772d4*y**7+0.28d3*y**8))/(0.4536d5*x**7)+  &
(y*(45360+0.7956d4*y-0.32274d6*y**2-0.58505d5*y**3+  &
0.8106d6*y**4+0.207158d6*y**5-0.205002d7*y**6+  &
0.2238855d7*y**7-0.115776d7*y**8+0.323336d6*y**9-  &
0.4704d5*y**10+0.28d4*y**11))/(0.54432d7*x**8)-(y**2*(45360+  &
0.7956d4*y-0.32274d6*y**2-0.58505d5*y**3+0.8106d6*y**4+  &
0.207158d6*y**5-0.205002d7*y**6+0.2238855d7*y**7-  &
0.115776d7*y**8+0.323336d6*y**9-0.4704d5*y**10+  &
0.28d4*y**11))/(0.54432d7*x**9)-(y*(5443200+0.54648d6*y-  &
0.37374084d8*y**2-0.3394248d7*y**3+0.82588957d8*y**4+  &
0.6408765d7*y**5-0.101524412d9*y**6-0.12413874d8*y**7+  &
0.154651321d9*y**8-0.150057435d9*y**9+0.72418826d8*y**10-  &
0.20401128d8*y**11+0.3409472d7*y**12-0.31416d6*y**13+  &
0.1232d5*y**14))/(0.3592512d9*x**10)



val  = 2*y*log(x) + log(dsum)

endif

end subroutine


! subroutine alegendre_expeval(expdata,dnu,dmu,t,aval,apval,bval1,bval2,valp,valq,a,b,c)
subroutine alegendre_expeval(expdata,dnu,dmu,t,bval1,valp,a,b,c)
implicit double precision (a-h,o-z)

double precision, intent(in)       :: dnu,t
type(alegendre_expansion_data)     :: expdata
double precision, intent(out)      :: bval1

double precision :: aval,apval,bval2

aval   = 0
apval  = 0
bval   = 0


k           = expdata%k
ifsmalldmu  = expdata%ifsmalldmu
ifover      = expdata%ifover
ifbetas     = expdata%ifbetas


if (ifover .eq. 1) then
dnu0   = 1/dnu
else
dnu0   = dnu
endif



!call alegendre_evalabc(ifsmalldmu,dnu,dmu,a,b,c)
call compute_dmu0(ifsmalldmu,dnu,dmu,dmu0)



if (t .ge. a) then

u = (t-a)/(b-a)
call alegendre_findint(expdata%nintsef,expdata%ef,dnu0,intef,e0,f0)
call alegendre_findint(expdata%nintscd,expdata%cd,dmu0,intcd,c0,d0)
call alegendre_findint(expdata%nintsab,expdata%ab,u,intab,a0,b0)


iptr1  = expdata%iptrsalpha(intab,intcd,intef)
iptr2  = expdata%iptrsalphap(intab,intcd,intef)

call alegendre_tensor_eval2(expdata%ncoefsalpha,expdata%coefsalpha,expdata%ncoefsalphap, &
  expdata%coefsalphap,iptr1,iptr2,a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,aval,apval)


if (ifover .eq. 1) then
aval  = aval  * dnu
apval = apval * dnu
endif

dd   =  sqrt((1+2*dnu)/(pi*apval))
valp =  cos(aval ) * dd
! valq = -sin(aval ) * dd


else

u = (t-c)/(a-c)


call alegendre_findint(expdata%nintsef,expdata%ef,dnu0,intef,e0,f0)
call alegendre_findint(expdata%nintscd,expdata%cd,dmu0,intcd,c0,d0)
call alegendre_findint(expdata%nintsabb,expdata%abb,u,intabb,a0,b0)


iptr1  = expdata%iptrsbeta1(intabb,intcd,intef)
iptr2  = expdata%iptrsbeta2(intabb,intcd,intef)

call alegendre_tensor_eval2(expdata%ncoefsbeta1,expdata%coefsbeta1, &
  expdata%ncoefsbeta2,expdata%coefsbeta2,iptr1,iptr2, &
  a0,b0,c0,d0,e0,f0,u,dmu0,dnu0,bval1,bval2)

if (ifover .eq. 1) then
bval1  = bval1  * dnu
! bval2  = bval2  * dnu
endif

bval1 = bval1+dnu
! bval2 = bval2-dnu

valp = exp(bval1)
! valq = exp(bval2)

endif

end subroutine


subroutine compute_dmu0(ifsmalldmu,dnu,dmu,dmu0)
implicit double precision (a-h,o-z)


if (ifsmalldmu .eq. 1 ) then
dmu0 = dmu
else
dmu0 = (dmu-1)/(dnu-1)
endif

end subroutine


subroutine alegendre_findint(nints,ab,x,int,a,b)
implicit double precision (a-h,o-z)

integer             :: int,nints
double precision    :: ab(2,nints),x,a,b

integer             :: int0

!eps0 = epsilon(0.0d0)

! intl   = 1
! intr   = nints

do int = 1,nints-1
b = ab(2,int)
if (x .le. b) exit
end do

int0 = int
a = ab(1,int)
b = ab(2,int)

end subroutine


subroutine alegendre_tensor_eval2(ncoefs1,coefs1,ncoefs2,coefs2,iptr1,iptr2, &
  a,b,c,d,e,f,x,y,z,val1,val2)
implicit double precision (a-h,o-z)

double precision :: coefs1(ncoefs1),coefs2(ncoefs2)
double precision :: polsx(0:100),polsy(0:100),polsz(0:100)

nz1     = coefs1(iptr1)
ny1     = coefs1(iptr1+1)
nx1     = coefs1(iptr1+2)

nz2     = coefs2(iptr2)
ny2     = coefs2(iptr2+1)
nx2     = coefs2(iptr2+2)

nx      = max(nx1,nx2)
ny      = max(ny1,ny2)
nz      = max(nz1,nz2)

xx      = (x-(b+a)/2)*2/(b-a)
yy      = (y-(d+c)/2)*2/(d-c)
zz      = (z-(f+e)/2)*2/(f-e)


call alegendre_chebs(xx,nx,polsx)
call alegendre_chebs(yy,ny,polsy)
call alegendre_chebs(zz,nz,polsz)

iptr      = iptr1+3
val1      = 0

do k=0,nz1
n     = coefs1(iptr)
iptr  = iptr+1
do j=0,n
m     = coefs1(iptr)
iptr  = iptr+1
do i=0,m
val1  = val1 + coefs1(iptr)*polsx(i)*polsy(j)*polsz(k)
iptr  = iptr+1
end do
end do
end do

iptr      = iptr2+3
val2      = 0

do k=0,nz2
n     = coefs2(iptr)
iptr  = iptr+1
do j=0,n
m     = coefs2(iptr)
iptr  = iptr+1
do i=0,m
val2  = val2 + coefs2(iptr)*polsx(i)*polsy(j)*polsz(k)
iptr  = iptr+1
end do
end do
end do

end subroutine


subroutine alegendre_chebs(x,n,pols)
implicit double precision (a-h,o-z)

integer          :: n 
double precision :: pols(0:n),x

!
!  Evaluate the Chebyshev polynomials of degree 0 through n at a specified point
!  using the standard 3-term recurrence relation.
!
!  Input parameters:
!
!    x - point at which to evaluate the polynomials
!    n - the order of the polynomials to evaluate
!
!  Output parameters:
!
!    pols - this user-supplied and allocated array of length n+1 will
!      contain the values of the polynomials of order 0 through n upon
!      return
!

if (x .eq. 1.0d0) then
do i=0,n
pols(i) = 1.0d0
end do
return
endif

if (x .eq. -1.0d0) then
pols(0) = 1.0d0
do i=1,n
pols(i) = -pols(i-1)
end do
return
endif

pols(0) = 1.0d0
if (n .eq. 0) return

pols(1) = x
if (n .eq. 1) return

xx1 = 1.0d0
xx2 = x

do i=1,n-1
xx        = 2*x*xx2 - xx1
pols(i+1) = xx
xx1       = xx2
xx2       = xx
end do

end subroutine


! subroutine alegendre_macdonald(dnu,dmu,t,vallogp,vallogq,valp,valq)
! implicit double precision (a-h,o-z)
subroutine alegendre_macdonald(dnu,dmu,t,vallogp,valp)
implicit double precision (a-h,o-z)

! double precision :: valsj(0:15),valslogj(0:15),valsratj(0:15),dsignsj(0:15)

!
!  Evaluate the logarithms of the SCALED, NORMALIZED associated Legendre
!  functions (1) and (2) when dnu is positive and large and t is small using 
!  Macdonald's asymptotic expansions (see Proc. Royal Soc. of London, 1914 
!  pgs. 221-222).  
!
!  These expansions are valid for all 0 <= dmu <= dnu and small t.
!
!  Care has been taken to reduce the potential for numerical overflow 
!  and underflow.
!
!  Input parameters:
!    dnu - the degree of the Legendre function to evaluate; dnu >=0
!    dmu - the order of the Legendre function to evaluate; 0<= dmu <= dnu.
!    t - the argument; t must be small
!
!  Output parameters:
!    vallogp - the logarithm of the function (1)
!    vallogq - the logarithm of the function (2)
!

! dimension xs(30), vals(30)

! data xs  / -0.100000000000000000000000000000000000D+01,  &
!            -0.959492973614497389890368057066327693D+00,  &
!            -0.841253532831181168861811648919367640D+00,  &
!            -0.654860733945285064056925072466293503D+00,  &
!            -0.415415013001886425529274149229623161D+00,  &
!            -0.142314838273285140443792668616369629D+00,  &
!             0.142314838273285140443792668616369701D+00,  &
!             0.415415013001886425529274149229623209D+00,  &
!             0.654860733945285064056925072466293599D+00,  &
!             0.841253532831181168861811648919367736D+00,  &
!             0.959492973614497389890368057066327693D+00,  &
!             0.100000000000000000000000000000000000D+01,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00,  &
!             0.000000000000000000000000000000000000D+00 /


!
!  Evaluate P
! 

call alegendre_macdonald_logp0(dnu,dmu,t,vallogp)


!
!  Now evaluate Q_dnu^{-dmu} using the connection formula and interpolation
!

! diff    = 1-2*(dmu - nint(dmu))
! dmu0    = nint(2*dmu)/2.0d0
! diff    = abs(dmu-dmu0)


! if (diff .gt. 0.005d0) then
! call alegendre_macdonald_logq0(dnu,dmu,t,vallogq0,dsignq0)
!dd      = 1.0d0/cos(dmu*pi)* dsignq0  + pi/2 * tan(dmu*pi) * exp(vallogp-vallogq0)
! dd      = 1.0d0/cos(dmu*pi)* dsignq0  + tan(dmu*pi) * exp(vallogp-vallogq0)
! if (dd .lt. 0) dd = -dd
! vallogq = vallogq0 + log(dd)
! else
! a = dmu0-0.1d0
! b = dmu0+0.1d0
! n = 12

! if (n .eq. 1) then
! xs(1) = 0.0d0
! else
! h = pi/(n-1)
! do i=1,n
! xs(n-i+1) = cos(h*(i-1))
! end do
! endif


! do i=1,n
! x0    = xs(i) *(b-a)/2 + (b+a)/2
! call alegendre_macdonald_logq0(dnu,x0,t,vallogq0,dsignq0)
! call alegendre_macdonald_logp0(dnu,x0,t,vallogp0)
!dd      = 1.0d0/cos(x0*pi)* dsignq0  + pi/2 * tan(x0*pi) * exp(vallogp0-vallogq0)
! dd      = 1.0d0/cos(x0*pi)* dsignq0  + tan(x0*pi) * exp(vallogp0-vallogq0)
! if (dd .lt. 0) dd = -dd

! vallogq = vallogq0 + log(dd)
! vals(i) = vallogq
! end do

!
!  Interpolate using the barycentric formula
!

! xx   = (2*dmu - (b+a) ) /(b-a)

! sum1 = 0
! sum2 = 0
! dd1  = 1.0d0
! do i=1,n
! dd=1.0d0
! if (i .eq. 1 .OR. i .eq. n) dd = 0.5d0
! diff = xx-xs(i)

!
!  Handle the case in which the target node coincides with one of
!  of the Chebyshev nodes.
!
! if(abs(diff) .le. eps0) then
! val = vals(i)
! exit
! endif

!
!  Otherwise, calculate the sums.
!

! dd   = (dd1*dd)/diff
! dd1  = - dd1
! sum1 = sum1+dd*vals(i)
! sum2 = sum2+dd
! dd   = - dd
! end do

! vallogq    = sum1/sum2

! endif

!
!  Calculate the values of the functions
!

valp = exp(vallogp) 
! valq = exp(vallogq)

end subroutine


subroutine alegendre_macdonald_logp0(dnu,dmu,t,vallogp)
implicit double precision (a-h,o-z)

!
!  Evaluate the logarithm of the NORMALIZED, SCALED associated Legendre function 
!  of the first kind (1) using Macdonald's asymptoptic expansion. 
!
!  Care has been taken to minimize the potential for numerical overflow 
!  and underflow.
!


double precision :: valsj(0:15),valslogj(0:15),valsratj(0:15),dsignsj(0:15)

x  = (2*dnu+1)   * sin(t/2)

!
!  Evaluate the ratios J_{\mu+k}(x) / J_{\mu}(x)
!


do i=0,6
call alegendre_log_besselj(dmu+i,x,valslogj(i),dsignsj(i),valsj(i))
end do


do i=0,6
valsratj(i) = dsignsj(i)*exp(valslogj(i)-valslogj(0))
end do


!
!  Compute the logarithm of P 
!

term0  = valsratj(0)


term1  = valsratj(1)/(0.2d1*x)-0.1d1*valsratj(2)+(x*valsratj(3))/0.6d1
dterm1 = (0.3d1*x*valsratj(0)-0.6d1*(1+x**2)*valsratj(1)+x*((-3+  &
x**2)*valsratj(2)+x*(0.8d1*valsratj(3)-  &
0.1d1*x*valsratj(4))))/(0.12d2*x**2)
term1  = sin(t/2)**2 * term1


term2 = (9*valsratj(2))/(0.8d1*x**2)-(29*valsratj(3))/(0.6d1*x)+  &
(31*valsratj(4))/0.12d2-(11*x*valsratj(5))/0.3d2+  &
(x**2*valsratj(6))/0.72d2
term2  = sin(t/2)**4 * term2


! term3 = (75*valsratj(3))/(0.16d2*x**3)-(751*valsratj(4))/(0.24d2*x**2)+  &
! (1381*valsratj(5))/(0.48d2*x)-(1513*valsratj(6))/0.18d3+  &
! (4943*x*valsratj(7))/0.504d4-(17*x**2*valsratj(8))/0.36d3+  &
! (x**3*valsratj(9))/0.1296d4
! term3  = sin(t/2)**6 * term3


! term4 = (3675*valsratj(4))/(0.128d3*x**4)-(20877*valsratj(5))/(0.8d2*x**3)  &
! +(98683*valsratj(6))/(0.288d3*x**2)-  &
! (110177*valsratj(7))/(0.72d3*x)+(610843*valsratj(8))/0.2016d5-  &
! (44887*x*valsratj(9))/0.1512d5+(67117*x**2*valsratj(10))/0.4536d6-  &
! (23*x**3*valsratj(11))/0.648d4+(x**4*valsratj(12))/0.31104d5
! term4  = sin(t/2)**8 * term4


! term5 = (59535*valsratj(5))/(0.256d3*x**5)-  &
! (1714703*valsratj(6))/(0.64d3*x**4)+  &
! (52827053*valsratj(7))/(0.1152d5*x**3)-  &
! (3985691*valsratj(8))/(0.144d4*x**2)+  &
! (94053343*valsratj(9))/(0.12096d6*x)-  &
! (7002137*valsratj(10))/0.6048d5+  &
! (64804643*x*valsratj(11))/0.66528d7-  &
! (643931*x**2*valsratj(12))/0.13608d7+  &
! (141703*x**3*valsratj(13))/0.108864d8-  &
! (29*x**4*valsratj(14))/0.15552d6+(x**5*valsratj(15))/0.93312d6
! term5  = sin(t/2)**10 * term5

dsum     = term0 + term1 + term2
vallogp = -dmu * log(cos(t/2)* ( dnu+0.5d0) ) + valslogj(0) + log(dsum)


!
! Apply the normalizing factor and scaling factors
!

call alegendre_gamma_ratio2(dnu,dmu,dd)
vallogp = vallogp + 0.5d0 * (log(dnu+0.5d0) + dd + log(sin(t)) )


end subroutine


subroutine alegendre_log_besselj(dmu,t,vallogj,dsignj,valj)
implicit double precision (a-h,o-z)
!
!  Evaluate log( | J_\dmu(t) | ) and  J_\dmu(t).
! 

if (dmu .gt. 0) then
call bessel_eval(dmu,t,alpha,alphader,vallogj,vallogy,valj,valy)
dsignj  =  1.0d0

if (vallogj .eq. 0) then
dsignj  = 1.0d0
vallogj = log(abs(valj))
if (valj .lt. 0) dsignj = -1.0d0
endif

return
endif

call bessel_eval(-dmu,t,alpha,alphader,vallogj,vallogy,valj,valy)

if (vallogj .eq. 0) then
valj   = cos(pi*dmu)*valj-sin(pi*dmu)*valy
dsignj = 1.0d0
if (valj .lt. 0) dsignj = -1.0d0
vallogj = log(abs(valj))
return
endif

dd     = cos(pi*dmu)*exp(vallogj-vallogy) - sin(pi*dmu)

if (dd .lt. 0) then
dd     = -dd
dsignj = -1.0d0
endif

vallogj = vallogyj + log ( dd ) 
valj    = dsignj*exp(vallogj)

end subroutine


function alegendre_gamma(x)
implicit double precision (a-h,o-z)
double precision :: alegendre_gamma

!
!  This routine exists to deal with a bug in certain versions of the gfortran
!  library.  The bug causes the sign of the gamma function of negative values 
!  to be calculated incorrectly in certain cases.
!

dgamma = gamma(x)
if (x .lt. 0) then
nn1 = floor(x)
if ((nn1/2)*2 .eq. nn1 .AND. dgamma .lt. 0) then
dgamma = -dgamma
endif
endif

alegendre_gamma = dgamma

end function


subroutine alegendre_eval_init(dsize)
implicit double precision (a-h,o-z)

!
!  Read the precomputed expansions from the binary file "alegendre_data.bin"
!
!  Input parameters:  
!    None
!
!  Output parameters:
!    dsize - a (pretty good) estimate of the  memory occupied by the table read
!     from the disk, in megabytes
!

!
!  Do nothing and return if the table is already loaded into memory.
!

if (ifloaded .eq. 1) return

!
!  Initialize bessel_eval, which is used in alegendre_macdonald.
!

call bessel_eval_init(dsize0)

!
! Read the table into memory, set the ifloaded flag and return the
! size of the table in megabytes.
!


iw = 200
open (iw, FILE = 'alegendre_data.bin1', form = 'UNFORMATTED', status = 'OLD', &
  access = 'stream', err = 1000)
call alegendre_read_expansion(iw,expdata1)
call alegendre_read_expansion(iw,expdata2)
call alegendre_read_expansion(iw,expdata3)
call alegendre_read_expansion(iw,expdata4)
close (iw)

dsize    = expdata1%dmemory + expdata2%dmemory + expdata3%dmemory + expdata4%dmemory 

iw = 201
open (iw, FILE = 'alegendre_data.bin2', form = 'UNFORMATTED', status = 'OLD', &
  access = 'stream', err = 2000)
call alegendre_read_expansion(iw,expdata5)
call alegendre_read_expansion(iw,expdata6)
close (iw)

dsize    = dsize + expdata5%dmemory + expdata6%dmemory 

!
!  Uncomment this to allow for degrees between 1,000,000 and 100,000,000
!

iw = 202
open (iw, FILE = 'alegendre_data.bin3', form = 'UNFORMATTED', status = 'OLD', &
  access = 'stream', err = 2000)
call alegendre_read_expansion(iw,expdata7)
call alegendre_read_expansion(iw,expdata8)
close (iw)
dsize    = dsize + expdata7%dmemory + expdata8%dmemory 


ifloaded  = 1

return

1000 continue

print *,"alegendre_eval_init: unable to open and/or read alegendre_data.bin1"
stop

2000 continue

print *,"alegendre_eval_init: unable to open and/or read alegendre_data.bin2"
stop

end subroutine


subroutine alegendre_read_expansion(iw,expdata)
implicit double precision (a-h,o-z)
type(alegendre_expansion_data), intent(out) :: expdata


call alegendre_read_integer_binary(iw,expdata%ifbetas)
call alegendre_read_integer_binary(iw,expdata%ifover)
call alegendre_read_integer_binary(iw,expdata%ifsmalldmu)
call alegendre_read_integer_binary(iw,expdata%ifinverse)

call alegendre_read_double_binary(iw,expdata%dnu1)
call alegendre_read_double_binary(iw,expdata%dnu2)
call alegendre_read_double_binary(iw,expdata%dtime)
call alegendre_read_double_binary(iw,expdata%dmemory)
call alegendre_read_double_binary(iw,expdata%epsdisc)
call alegendre_read_double_binary(iw,expdata%epscomp)
call alegendre_read_double_binary(iw,expdata%epsphase)
call alegendre_read_double_binary(iw,expdata%epsrequired)


call alegendre_read_integer_binary(iw,expdata%ncoefs)
call alegendre_read_integer_binary(iw,expdata%nintsab)
call alegendre_read_integer_binary(iw,expdata%nintscd)
call alegendre_read_integer_binary(iw,expdata%nintsef)
call alegendre_read_integer_binary(iw,expdata%ncoefsalpha)
call alegendre_read_integer_binary(iw,expdata%ncoefsalphap)
call alegendre_read_integer_binary(iw,expdata%ncoefsalphapp)

allocate(expdata%ab(2,expdata%nintsab))
allocate(expdata%cd(2,expdata%nintscd))
allocate(expdata%ef(2,expdata%nintsef))
allocate(expdata%coefsalpha(expdata%ncoefsalpha))
allocate(expdata%coefsalphap(expdata%ncoefsalphap))
allocate(expdata%coefsalphapp(expdata%ncoefsalphapp))

allocate(expdata%iptrsalpha(expdata%nintsab,expdata%nintscd,expdata%nintsef))
allocate(expdata%iptrsalphap(expdata%nintsab,expdata%nintscd,expdata%nintsef))
allocate(expdata%iptrsalphapp(expdata%nintscd,expdata%nintsef))

nn = expdata%nintsab*2
call alegendre_read_double_array_binary(iw,nn,expdata%ab)
nn = expdata%nintscd*2
call alegendre_read_double_array_binary(iw,nn,expdata%cd)
nn = expdata%nintsef*2
call alegendre_read_double_array_binary(iw,nn,expdata%ef)

call alegendre_read_double_array_binary(iw,expdata%ncoefsalpha,expdata%coefsalpha)
call alegendre_read_double_array_binary(iw,expdata%ncoefsalphap,expdata%coefsalphap)
call alegendre_read_double_array_binary(iw,expdata%ncoefsalphapp,expdata%coefsalphapp)

nn = expdata%nintsab * expdata%nintscd * expdata%nintsef
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalpha)
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalphap)
nn = expdata%nintscd * expdata%nintsef
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalphapp)

if (expdata%ifinverse .eq. 1) then
call alegendre_read_integer_binary(iw,expdata%nintsinv)
call alegendre_read_integer_binary(iw,expdata%ncoefsalphainv)

allocate(expdata%abinv(2,expdata%nintsinv))
allocate(expdata%coefsalphainv(expdata%ncoefsalphainv))
allocate(expdata%iptrsalphainv(expdata%nintsinv,expdata%nintscd,expdata%nintsef))

nn = expdata%nintsinv*2
call alegendre_read_double_array_binary(iw,nn,expdata%abinv)

call alegendre_read_double_array_binary(iw,expdata%ncoefsalphainv,expdata%coefsalphainv)
nn = expdata%nintsinv * expdata%nintscd * expdata%nintsef
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsalphainv)
endif

if (expdata%ifbetas .eq. 1) then
call alegendre_read_integer_binary(iw,expdata%nintsabb)
call alegendre_read_integer_binary(iw,expdata%ncoefsbeta1)
call alegendre_read_integer_binary(iw,expdata%ncoefsbeta2)

allocate(expdata%abb(2,expdata%nintsabb))
allocate(expdata%coefsbeta1(expdata%ncoefsbeta1))
allocate(expdata%coefsbeta2(expdata%ncoefsbeta2))
allocate(expdata%iptrsbeta1(expdata%nintsabb,expdata%nintscd,expdata%nintsef))
allocate(expdata%iptrsbeta2(expdata%nintsabb,expdata%nintscd,expdata%nintsef))

nn = expdata%nintsabb*2
call alegendre_read_double_array_binary(iw,nn,expdata%abb)
call alegendre_read_double_array_binary(iw,expdata%ncoefsbeta1,expdata%coefsbeta1)
call alegendre_read_double_array_binary(iw,expdata%ncoefsbeta2,expdata%coefsbeta2)
nn = expdata%nintsef *  expdata%nintscd * expdata%nintsabb
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsbeta1)
call alegendre_read_integer_array_binary(iw,nn,expdata%iptrsbeta2)
endif

end subroutine


subroutine alegendre_read_integer_array_binary(iw,n,idata)
implicit double precision (a-h,o-z)
integer :: idata(n)
integer*8 ix
do i=1,n
read(iw) ix
idata(i) = ix
end do
end subroutine


subroutine alegendre_read_integer_binary(iw,idata)
implicit double precision (a-h,o-z)
integer*8 i 
read (iw) i
idata = i
end subroutine


subroutine alegendre_read_double_binary(iw,data)
implicit double precision (a-h,o-z)
real*8 :: x
read (iw) x
data = x
end subroutine


subroutine alegendre_read_double_array_binary(iw,n,data)
implicit double precision (a-h,o-z)
double precision :: data(n)
real*8 x 
read (iw) data
end subroutine

end module