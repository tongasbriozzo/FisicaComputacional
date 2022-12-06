! ______________________________________________________________________________
! Ecuación de Calor
! Este programa emplea los métodos forward Euler, implicito y Crank-Nicolson
! para calcular la temperatura de una barra de aluminio de 1 m de longitud, a
! una temperatura inicial uniforme de 100 ºC con condiciones de contorno de
! Dirichlet nulas.
!
! Ésta es una variante del programa original, diseñana para calcular únicamente
! el error global a t=1800s y el tiempo de cpu requerido
!
! Gaston Briozzo, FaMAFyC, UNC, 22/04/2021
! ______________________________________________________________________________

program err
use precision, only: pr
implicit none

  integer                  :: i, j, k
  integer                  :: nx, nt, np
  real(pr)                 :: x, t, dx, dt
  real(pr)                 :: L, T0, Tb1, Tb2
  real(pr)                 :: Ka, C, rho, D, eta
  real(pr)                 :: ti, tf, tfe, timp, tcn
  real(pr)                 :: EA_fE, EA_imp, EA_CN, Tan
  real(pr)                 :: EG_fE, EG_imp, EG_CN
  real(pr), allocatable    :: T_fE(:,:), T_imp(:,:), T_CN(:,:)
  real(pr), allocatable    :: Tini(:), Tint(:), Tfin(:)

! Definimos los formatos que utilizaremos en los documentos de texto
  5 format(A2, 5x, A2, 5x, A9, 16x, A10, 15x, A9)
  6 format(I5, 2x, I5, 2x, 3(E23.16,2x))

! Inicializamos los parámetros
! La idea es cambiar los valores dx y dt y medir el error global y el tiempo de
! cpu correspondiente a cada uno

open(131,file='eg1800.d')  ! Error global a un dado nx y np
open(132,file='et1800.d')  ! Tiempo de cpu a un dado nx y np
write(131,5) 'nx', 'np', 'GloErr fE', 'GloErr imp', 'GloErr CN'
write(132,5) 'nx', 'np', 'Tiempo fE', 'Tiempo imp', 'Tiempo CN'

do  nx = 10, 300, 1       ! Número de divisiones del intervalo espacial
  nt = 1                  ! Vestigio del programa original
do  np = 180, 90000, 180  ! Numero de divisiones del intervalo temporal

  Ka  = 237._pr
  C   = 900._pr
  rho = 2700._pr

  D = Ka/C/rho

  L   = 1._pr
  T0  = 100._pr
  Tb1 = 0._pr
  Tb2 = 0._pr

  dx = L/nx
  dt = 1800._pr/real(np,pr)

  eta = dt*D/L/L/dx/dx

! Damos dimensión a los vectores
  allocate(T_fE(0:nx,0:nt), T_imp(0:nx,0:nt), T_CN(0:nx,0:nt), Tini(0:nx), Tint(0:nx), Tfin(0:nx))

! ______________________________________________________________________________
! Primero probamos con el método forward Euler
! ______________________________________________________________________________

! Inicializamos los vectores que usaremos
  T_fE(:,0)  = T0
  Tini(:)    = T0
  Tfin(:)    = T0
  T_fE(0,:)  = Tb1
  Tini(0)    = Tb1
  Tfin(0)    = Tb1
  T_fE(nx,:) = Tb2
  Tini(nx)   = Tb2
  Tfin(nx)   = Tb2

  call cpu_time(ti)

  do j = 1, nt
    do k = 1, np
      call mat_tri(1._pr-2._pr*eta,eta,eta,nx,Tini,Tfin)
      Tini(:) = Tfin(:)
    end do
    T_fE(:,j) = Tfin(:)
  end do

  call cpu_time(tf)
  tfe = tf-ti

! ______________________________________________________________________________
! Ahora probamos el método implicito
! ______________________________________________________________________________

! Inicializamos los vectores que usaremos
  T_imp(:,0)  = T0
  Tini(:)     = T0
  Tfin(:)     = T0
  T_imp(0,:)  = Tb1
  Tini(0)     = Tb1
  Tfin(0)     = Tb1
  T_imp(nx,:) = Tb2
  Tini(nx)    = Tb2
  Tfin(nx)    = Tb2

  call cpu_time(ti)

  do j = 1, nt
    do k = 1, np
      call mat_inv(1._pr+2._pr*eta,-eta,-eta,nx,Tini,Tfin)
      Tini(:) = Tfin(:)
    end do
    T_imp(:,j) = Tfin(:)
  end do

  call cpu_time(tf)
  timp = tf-ti

! ______________________________________________________________________________
! Ahora probamos el método Crank-Nicolson
! ______________________________________________________________________________

! Inicializamos los vectores que usaremos
  T_CN(:,0)  = T0
  Tini(:)    = T0
  Tint(:)    = T0
  Tfin(:)    = T0
  T_CN(0,:)  = Tb1
  Tini(0)    = Tb1
  Tint(0)    = Tb1
  Tfin(0)    = Tb1
  T_CN(nx,:) = Tb2
  Tini(nx)   = Tb2
  Tint(nx)   = Tb2
  Tfin(nx)   = Tb2

  call cpu_time(ti)

  do j = 1, nt
    do k = 1, np
      call mat_tri(2._pr/eta-2._pr, 1._pr, 1._pr,nx,Tini,Tint)
      call mat_inv(2._pr/eta+2._pr,-1._pr,-1._pr,nx,Tint,Tfin)
      Tini(:) = Tfin(:)
    end do
    T_CN(:,j) = Tfin(:)
  end do

  call cpu_time(tf)
  tcn = tf-ti

! ______________________________________________________________________________
! Ahora, calculamos el error global a t=1800s y anotamos el tiempo de cpu
! ______________________________________________________________________________

  call cpu_time(ti)

! Calculamos el error global y el tiempo de cpu a t=1800s
  j = 1
  EG_fE  = 0._pr
  EG_imp = 0._pr
  EG_CN  = 0._pr
    do i = 1, nx-1
      t = real(j,pr)*real(np,pr)*dt
      x = real(i,pr)*L/real(nx,pr)
      Tan = T_an(x,D*t/L/L,129,T0)
      EA_fE  = abs(T_fE(i,j)  - Tan)
      EA_imp = abs(T_imp(i,j) - Tan)
      EA_CN  = abs(T_CN(i,j)  - Tan)
      EG_fE  = EG_fe  + EA_fE/Tan
      EG_imp = EG_imp + EA_imp/Tan
      EG_CN  = EG_CN  + EA_CN/Tan
    end do

  write(131,6) nx, np, EG_fE/real(nx-1,pr), EG_imp/real(nx-1,pr), EG_CN/real(nx-1,pr)
  write(132,6) nx, np, tfe, timp, tcn

  call cpu_time(tf)
  write(*,*) 'Tiempo de Error a t fijo', (tf-ti), nx, np

! Quitamos dimension a los vectores
  deallocate(T_fE, T_imp, T_CN, Tini, Tint, Tfin)

end do
end do

close(131)
close(132)

CONTAINS

! ______________________________________________________________________________
! Función T_an calcula analíticamente la tenmeratura de una barra unidimensional
! homogenea en la posición x, al tiempo t, con una temperatura inicial uniforme
! T0 y temperatura cero en los extremos para todo t, a orden n.
! ______________________________________________________________________________

  function T_an(x,t,n,T0)
  implicit none
    real(pr), intent(in) :: x, t, T0
    integer, intent(in)  :: n
    real(pr)             :: T_an, sum
    real(pr), parameter  :: pi=acos(-1._pr)
    integer              :: i
    sum = 0._pr
    do i = 1, n, 2
      sum = sum + sin(i*pi*x)*exp(-i*i*pi*pi*t)/(i*pi)
    end do
    T_an = 4._pr*T0*sum
  end function T_an

! ______________________________________________________________________________
! Subrutina mat_tri calcula el producto entre una matriz tridiagonal, con
! coeficientes d en la diagonal, a en la diagonal inferior y c en la diagonal
! superior, con el vector Tini, resultando en el vector Tfin.
! Se concideran condiciones de contorno de Dirichlet nulas.
! ______________________________________________________________________________

  subroutine mat_tri(d,a,c,dim,Tini,Tfin)
  implicit none
    real(pr), intent(in)                   :: d, a, c
    integer, intent(in)                    :: dim
    real(pr), intent(in), dimension(0:dim) :: Tini
    real(pr), dimension(0:dim)             :: Tfin
    integer                                :: i

!   Calculamos el vector resultante
    do i = 1, dim-1
      Tfin(i) = a*Tini(i-1) + d*Tini(i) + c*Tini(i+1)
    end do

!   Añadimos las condiciones de contorno
    Tfin(0)   = Tini(0)
    Tfin(dim) = Tini(dim)

  end subroutine mat_tri

! ______________________________________________________________________________
! Subrutina mat_inv calcula el producto entre la inversa de una matriz
! tridiagonal, con coeficientes d en la diagonal, a en la diagonal inferior
! y c en la diagonal superior, con el vector Tini, resultando en el vector Tfin.
! Se concideran condiciones de contorno de Dirichlet nulas.
! ______________________________________________________________________________

  subroutine mat_inv(d,a,c,dim,Tini,Tfin)
  implicit none
    real(pr), intent(in)                   :: d, a, c
    integer, intent(in)                    :: dim
    real(pr), intent(in), dimension(0:dim) :: Tini
    real(pr), dimension(0:dim)             :: Tfin
    real(pr), dimension(1:dim-1)           :: h, p
    integer                                :: i

!   Inicializamos los parámetros h y p
    h(1) = c/d
    p(1) = Tini(1)/d
    do i = 2, dim-1
      h(i) = c/(d-a*h(i-1))
      p(i) = (Tini(i)-a*p(i-1))/(d-a*h(i-1))
    end do

!   Calculamos el vector resultante
    Tfin(dim-1) = p(dim-1)
    do i = dim-2, 1, -1
      Tfin(i) = p(i) - h(i)*Tfin(i+1)
    end do

!   Añadimos las condiciones de contorno
    Tfin(0)   = Tini(0)
    Tfin(dim) = Tini(dim)

  end subroutine mat_inv

end program err
