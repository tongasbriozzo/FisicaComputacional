!*******************************************************************************
!
!  Programa Tiempo de Auto Correlación
!
!  Este programa muestrea la temperatura cinética, las energías cinética,
!  potencial y total y la presión instantaneas en un sistema de Dinámica
!  Browneana con potencial de Lennard-Jones durante la termalización y el
!  posterior equilibrio termodinámico. Obtiene ademas valores medios y varianzas
!  además de unh histograma de velocidades.
!
!  Gaston Briozzo, Física Computacional, FaMAFyC, UNC, 07/07/2021
!
!*******************************************************************************

program run
use precision, only: pr
use mtmod,     only: grnd
use dbmod
implicit none

  integer                      :: i, nh, ns
  integer, parameter           :: nl=8, N=nl*nl*nl, t_eq=2500, nsim=100000, nhis=50
  real(pr)                     ::  Ti,  Ec,  Ep,  Et,  P
  real(pr)                     ::  Ti0, Ec0, Ep0, Et0, P0
  real(pr)                     ::  Ti2, Ec2, Ep2, Et2, P2
  real(pr), dimension(1:nhis)  ::  Tit, Ect, Ept, Ett, Pt
  real(pr), dimension(1:nhis)  :: HTi, HEc, HEp, HEt, HP
  real(pr), dimension(1:nhis)  :: CTi, CEc, CEp, CEt, CP
  real(pr), dimension(1:N,1:3) :: x, v, xnew
  real(pr), parameter          :: pi=acos(-1._pr)
  real(pr), parameter          :: dt=0.002_pr
  real(pr), parameter          :: r_cut = 2.5_pr, r2cut=r_cut*r_cut, r_inv=1._pr/r_cut, r_i3=r_inv*r_inv*r_inv
  real(pr), parameter          :: T=2.5_pr, rho=0.8_pr, eta=2.87_pr
  real(pr), parameter          :: Vol=real(N,pr)/rho, L=(Vol)**(1._pr/3._pr), D0=T/(3._pr*pi*eta)
  real(pr), parameter          :: f_fx = dt*D0/T, f_nx = sqrt(2._pr*D0*dt), f_fv = D0/T, f_nv = sqrt(2._pr*D0)
  real(pr), parameter          :: U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  real(pr), parameter          :: P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr

!*******************************************************************************
! Inicializamos las variables
!*******************************************************************************

! Inicializamos los valores de expectación <A>
  Ti = 0._pr
  Ec = 0._pr
  Ep = 0._pr
  Et = 0._pr
  P  = 0._pr
! Inicializamos los valores iniciales <A(0)>
  Ti0 = 0._pr
  Ec0 = 0._pr
  Ep0 = 0._pr
  Et0 = 0._pr
  P0  = 0._pr
! Inicializamos los valores de expectación <A^2>
  Ti2 = 0._pr
  Ec2 = 0._pr
  Ep2 = 0._pr
  Et2 = 0._pr
  P2  = 0._pr
! Inicializamos los valores iniciales <A(0)A(t)>
  Tit = 0._pr
  Ect = 0._pr
  Ept = 0._pr
  Ett = 0._pr
  Pt  = 0._pr
! Inicializamos los historiales de valores
  HTi = 0._pr
  HEc = 0._pr
  HEp = 0._pr
  HEt = 0._pr
  HP  = 0._pr
! Inicializamos las funciones autocorrelación
  CTi = 0._pr
  CEc = 0._pr
  CEp = 0._pr
  CEt = 0._pr
  CP  = 0._pr

!*******************************************************************************
! Termalizamos el sistema
!*******************************************************************************

  call SC(N,L,nl,x)                                       ! Acomodamos las partículas según una red SC
  do i = 1, t_eq                                          ! Nos movemos hasta el paso t_eq
    call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)                                    ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
    call pbc(N,L,x)                                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
  end do

!*******************************************************************************
! Llenamos los vectores historia
!*******************************************************************************

  do nh = 1, nhis
    call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v)
    x(:,:) = xnew(:,:)
    call pbc(N,L,x)
    HTi(nh) = Temperatura(N,v)                                    ! Temperatura cinética instantanea
    HEc(nh) = K(N,v)                                              ! Energía cinética instantanea
    HEp(nh) = U(N,L,r2cut,x) + N*U_tail                           ! Energía potencial instantanea
    HEt(nh) = HEc(nh) + HEp(nh)                                   ! Energía total instantanea
    HP(nh)  = rho*HTi(nh) + pf(N,L,r2cut,x)/(3._pr*Vol) + P_tail  ! Presión instantanea
  end do

!*******************************************************************************
! Calculamos la función Autocorrelación promediando sobre nsim valores
!*******************************************************************************

  do ns = 1, nsim
!  Tomamos los valores al tiempo inicial
    Ti0 = HTi(1)
    Ec0 = HEc(1)
    Ep0 = HEp(1)
    Et0 = HEt(1)
    P0  = HP(1)
!   Sumamos para obtener los valores de expectacion <A> y <A^2>
    Ti  = Ti  + HTi(1)
    Ti2 = Ti2 + HTi(1)**2
    Ec  = Ec  + HEc(1)
    Ec2 = Ec2 + HEc(1)**2
    Ep  = Ep  + HEp(1)
    Ep2 = Ep2 + HEp(1)**2
    Et  = Et  + HEt(1)
    Et2 = Et2 + HEt(1)**2
    P   = P   + HP(1)
    P2  = P2  + HP(1)**2
!   Avanzamos un paso y actualizamos los historiales
    call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v)
    x(:,:) = xnew(:,:)
    call pbc(N,L,x)
    do nh = 1, nhis-1
      HTi(nh) = HTi(nh+1)
      HEc(nh) = HEc(nh+1)
      HEp(nh) = HEp(nh+1)
      HEt(nh) = HEt(nh+1)
      HP(nh)  = HP(nh+1)
    end do
    HTi(nhis) = Temperatura(N,v)
    HEc(nhis) = K(N,v)
    HEp(nhis) = U(N,L,r2cut,x) + N*U_tail
    HEt(nhis) = HEc(nhis) + HEp(nhis)
    HP(nhis)  = rho*HTi(nhis) + pf(N,L,r2cut,x)/(3._pr*Vol) + P_tail
!   Actualizada la historia, obtenemos los valores de expectacion <A(0)A(t)>
    do nh = 1, nhis
      Tit(nh) = Tit(nh) + Ti0*HTi(nh)
      Ect(nh) = Ect(nh) + Ec0*HEc(nh)
      Ept(nh) = Ept(nh) + Ep0*HEp(nh)
      Ett(nh) = Ett(nh) + Et0*HEt(nh)
      Pt(nh)  = Pt(nh ) + P0*HP(nh)
    end do
  end do

! Normalizamos
  Ti  = Ti/real(nsim,pr)
  Ec  = Ec/real(nsim,pr)
  Ep  = Ep/real(nsim,pr)
  Et  = Et/real(nsim,pr)
  P   = P /real(nsim,pr)
  Ti2 = Ti2/real(nsim,pr)
  Ec2 = Ec2/real(nsim,pr)
  Ep2 = Ep2/real(nsim,pr)
  Et2 = Et2/real(nsim,pr)
  P2  = P2 /real(nsim,pr)
  Tit = Tit/real(nsim,pr)
  Ect = Ect/real(nsim,pr)
  Ept = Ept/real(nsim,pr)
  Ett = Ett/real(nsim,pr)
  Pt  = Pt /real(nsim,pr)

! Anotamos los resultados
  open(250,file='tac.d')
  do nh = 1, nhis
    CTi(nh)=(Tit(nh)-Ti*Ti)/(Ti2-Ti*Ti)
    CEc(nh)=(Ect(nh)-Ec*Ec)/(Ec2-Ec*Ec)
    CEp(nh)=(Ept(nh)-Ep*Ep)/(Ep2-Ep*Ep)
    CEt(nh)=(Ett(nh)-Et*Et)/(Et2-Et*Et)
    CP(nh) =(Pt (nh)-P *P )/(P2 -P *P)
    write(250,*) nh*dt, CTi(nh), CEc(nh), CEp(nh), CEt(nh), CP(nh)
  end do
  close(250)

end program run
