! ______________________________________________________________________________
! Ecuación de Calor
! Este programa emplea los métodos forward Euler, implicito y Crank-Nicolson
! para calcular la temperatura de una barra de aluminio de 1 m de longitud, a
! una temperatura inicial uniforme de 100 ºC con condiciones de contorno de
! Dirichlet nulas.
!
! Ésta es una variante del programa original, diseñana para calcular únicamente
! el error relativo a t=1800s y el tiempo de cpu requerido
!
! Gaston Briozzo, FaMAFyC, UNC, 22/04/2021
! ______________________________________________________________________________

program test
use precision, only: pr
implicit none

  integer                  :: i, j, k
  integer                  :: nx, nt, np
  real(pr)                 :: x, t, dx, dt
  real(pr)                 :: L, T0, Tb1, Tb2
  real(pr)                 :: Ka, C, rho, D, eta
  real(pr)                 :: ti, tf
  real(pr)                 :: EA_fE, EA_imp, EA_CN, ER_fE, ER_imp, ER_CN, Tan
  real(pr), allocatable    :: T_fE(:,:), T_imp(:,:), T_CN(:,:)
  real(pr), allocatable    :: Tini(:), Tint(:), Tfin(:)

! Definimos los formatos que utilizaremos en los documentos de texto
  5 format(A8, 17x, A9, 16x, A10, 15x, A9)
  6 format(4(E23.16,2x))

! Inicializamos los parámetros
! La idea es cambiar los valores dx y dt hasta obtener un error relativo de
! 10^-4 en el menor tiempo posible para cada método
  nx = 100
  nt = 1
  np = 6000

  Ka  = 237._pr
  C   = 900._pr
  rho = 2700._pr

  D = Ka/C/rho

  L   = 1._pr
  T0  = 100._pr
  Tb1 = 0._pr
  Tb2 = 0._pr

  dx = L/nx
  dt = 0.3_pr

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
    do i = 1, nx-1
      T_fE(i,j) = Tfin(i)
    end do
  end do

  call cpu_time(tf)
  write(*,*) 'Tiempo de forward Euler', (tf-ti)

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
    do i = 1, nx-1
      T_imp(i,j) = Tfin(i)
    end do
  end do

  call cpu_time(tf)
  write(*,*) 'Tiempo de Implicito', (tf-ti)

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
    do i = 1, nx-1
      T_CN(i,j) = Tfin(i)
    end do
  end do

  call cpu_time(tf)
  write(*,*) 'Tiempo Crank-Nicolson', (tf-ti)

! ______________________________________________________________________________
! Ahora, trazamos las curvas de error a t=1800s
! ______________________________________________________________________________

  call cpu_time(ti)

  open(130,file='er1800.d')
  write(130,5) 'Posicion', 'RelErr fE', 'RelErr imp', 'RelErr CN'

! Calculamos el err abs y rel a t=1800s
  j = 1
    do i = 1, nx-1
      t = real(j,pr)*real(np,pr)*dt
      x = real(i,pr)*L/real(nx,pr)
      Tan = T_an(x,D*t/L/L,129,T0)
      EA_fE  = abs(T_fE(i,j)  - Tan)
      EA_imp = abs(T_imp(i,j) - Tan)
      EA_CN  = abs(T_CN(i,j)  - Tan)
      ER_fE  = EA_fE/Tan
      ER_imp = EA_imp/Tan
      ER_CN  = EA_CN/Tan
      write(130,6) x, ER_fE, ER_imp, ER_CN
    end do

  close(130)

  call cpu_time(tf)
  write(*,*) 'Tiempo de Error a t fijo', (tf-ti)

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

end program test
