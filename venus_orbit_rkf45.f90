program venus_orbit_rkf45
implicit none
 
! Veenuse orbiidi simulatsioon (2D) RKF45 meetodiga.
! ODE süsteem: y = [x, y, vx, vy]
!   x'  = vx
!   y'  = vy
!   vx' = -mu*x/r^3
!   vy' = -mu*y/r^3
! kus r = sqrt(x^2 + y^2), mu = G*M.
! Integratsioon tehakse teegialamprogrammiga rkf45 (srkf45.for).
 
! Tehakse 2 simulatsiooni:
!   1) v = 1.0 * V1  -> (peaaegu) ringorbiit
!   2) v = 1.2 * V1  -> elliptilisem orbiit
! Tulemused salvestatakse eraldi failidesse:
!   orbit_v1.dat  ja  orbit_v12.dat

! ODE dimensioon
integer, parameter :: neqn = 4

! Olek ja RKF45 töömassivid (rkf45 nõuab work >= 3+6*neqn ja iwork >= 5)
real    :: y(neqn), work(3+6*neqn)
integer :: iwork(5), iflag

! Aeg, väljundi samm ja tolerantsid
integer :: i, nt
real    :: t, tout, dt, relerr, abserr

! Füüsikalised konstandid ja algtingimuste parameetrid
real :: G, Mv, Rv, mu, h0, v1, vfactor

! Abimuutujad (diagnostika/energia)
real :: r, vabs, Ekin, Epot, Etot

! Kaks juhtumit (V1 ja 1.2*V1)
real, dimension(2) :: vfactors
character(len=13), dimension(2) :: fnames
integer :: icase

! Jagame gravitatsiooniparameetri mu alamprogrammiga f common kaudu
common /params/ mu

! Ütle kompilaatorile, et f on alamprogramm (external), mida rkf45 kutsub
external f

! 1) Konstandid (SI ühikud)
G  = 6.67430e-11     ! gravitatsioonikonstant, m^3/(kg*s^2)
Mv = 4.8675e24       ! Veenuse mass, kg
Rv = 6051800.        ! Veenuse raadius, m
mu = G*Mv            ! mu = G*M, m^3/s^2

! 2) Algtingimused: h0 = 200 km kõrgus
! Positsioon: (x0, y0) = (Rv+h0, 0)
! Kiirus: tangentsiaalne (vx0=0, vy0=vfactor*V1)
h0 = 200000.

! 3) Integreerimise seaded
! dt = väljundi samm
dt = 1.
nt = 15000           ! koguaeg ~15000 s
relerr = 1.e-8
abserr = 1.e-6

! Kahe juhtumi definitsioon
vfactors = (/ 1.0, 1.2 /)
fnames = (/ "orbit_v1.dat ", "orbit_v12.dat" /)

! 4) Jooksuta kaks simulatsiooni
do icase = 1, 2
vfactor = vfactors(icase)

! reset algtingimused iga case jaoks
t     = 0.
iflag = 1

y(1) = Rv + h0        ! x0
y(2) = 0.             ! y0

v1   = sqrt(mu / y(1)) ! V1 kõrgusel r0 = Rv+h0
y(3) = 0.              ! vx0
y(4) = vfactor * v1    ! vy0

! Ava väljundfail (iga case eraldi fail)
open(20, file=trim(fnames(icase)))
write(20,'(A)') "# t  x  y  vx  vy  r  v  Ekin  Epot  Etot (energiad per 1 kg)"
write(20,'(A,F6.2,A,ES12.4)') "# vfactor = ", vfactor, " ; V1 (m/s) = ", v1

! Integreerimine: iga iteratsioon viib t -> tout=t+dt
do i = 1, nt
tout = t + dt

call rkf45(f, neqn, y, t, tout, relerr, abserr, iflag, work, iwork)

! Kui rkf45 tagastab muu lipu, jätkame standardrežiimis
if (iflag /= 2) iflag = 2

! Kindluse mõttes (rkf45 viib t tavaliselt tout-ni)
t = tout

! Diagnostika: kaugus ja kiirus
r    = sqrt(y(1)*y(1) + y(2)*y(2))
vabs = sqrt(y(3)*y(3) + y(4)*y(4))

! Energia per 1 kg: Ekin = 0.5*v^2, Epot = -mu/r
Ekin = 0.5 * vabs*vabs
Epot = -mu / r
Etot = Ekin + Epot

! Salvesta rida faili
write(20,'(10ES15.6)') t, y(1), y(2), y(3), y(4), r, vabs, Ekin, Epot, Etot
end do

close(20)
end do

print *, "Valmis: orbit_v1.dat ja orbit_v12.dat"
end program venus_orbit_rkf45


! Derivaatide alamprogramm RKF45 jaoks
! Sisend:  t, y=[x,y,vx,vy]
! Väljund: yp = dy/dt
subroutine f(t, y, yp)
implicit none
real :: t, y(4), yp(4)
real :: mu, r

common /params/ mu

! r = kaugus Veenuse keskpunktist
r = sqrt(y(1)*y(1) + y(2)*y(2))

! ODE:
yp(1) = y(3)                 ! x'  = vx
yp(2) = y(4)                 ! y'  = vy
yp(3) = -mu * y(1) / (r**3)  ! vx' = -mu*x/r^3
yp(4) = -mu * y(2) / (r**3)  ! vy' = -mu*y/r^3
end subroutine f



