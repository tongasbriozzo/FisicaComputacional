!*******************************************************************************
!
! Programa p1c
!
! Este programa calcula la magnetizacion de un sistema de ns por ns spines con
! condiciones de contorno periodicas, en funcion de la temperatura, y realiza
! un histograma de los valores obtenidos
!
! Gaston Briozzo, FaMAFyC, UNC, 18/05/2020
!
!*******************************************************************************

program p1c
use precision, only: pr
use mtmod,     only: grnd
implicit none

  integer                       :: i, j, k
  integer                       :: nf, nc
  integer, parameter            :: ns=40, np=100, nsim=1000000, ndes=100
  integer, dimension(1:ns,1:ns) :: s
  real(pr)                      :: V, T, M, det
  real(pr), dimension(1:100)    :: y

! Definimos los formatos que emplearemos en los documentos de texto
  1 format(A13,   12x, A12)
  2 format(E23.16, 2x, E23.16)

! calculamos el volumen del sistema
  V = real(ns*ns,pr)

! Calculamos M para una dada T, y realizamos el histograma sobre nsim simulaciones
! Definimos la temperatura
  T = 2.50_pr
! Inicializamos la Magnetizacion
  M = 0._pr
! Inicializamos el histograma
  y(:) = 0._pr
!  Elegimos la posicion inicial de los spines al azar
  call inicializar(s)
!  Descartamos los primeros ndes pasos de MC, hasta alcanzar el equilibrio
  do j = 1, ndes                    ! Pasos descartados
    do nf = 1, ns                   ! Numero de fila
      do nc = 1, ns                 ! Numero de columna
        call metropolis(nf,nc,s,T)  ! Algritmo de Metrópolis
      end do
    end do
  end do
! Ahora, promediamos sobre nsim datos
  do k = 1, nsim
!   Comprobamos la evolucion del programa
    write(*,*) T, k
!   Dejamos que el sistema evolucione np pasos de MC para disminuir la correlacion entre mediciones
    do j = 1, np
      do nf = 1, ns
        do nc = 1, ns
          call metropolis(nf,nc,s,T)
        end do
      end do
    end do
!   Calculamos la magnetizacion
    M = abs(mu(s))/V  ! <mu>
    write(*,*) M
    det = real(100,pr)*M
    do i = 1, 100
      if (det<real(i,pr)) then
        y(i) = y(i) + 1._pr
        exit
      end if
    end do
  end do
! Normalizamos
  y(:) = y(:)/real(nsim,pr)
! Abrimos el documento de dexto y describimos los datos a anotar
  open(430,file='his.d')
  write(430,1) 'Magnetizacion', 'Probabilidad'
! Anotamos las probabilidades de obtener cada magnetizacion
  do i = 1, 100
    write(430,2) 0.01_pr*real(i-1,pr), 0._pr
    write(430,2) 0.01_pr*real(i-1,pr), y(i)
    write(430,2) 0.01_pr*real(i,pr),   y(i)
  end do

CONTAINS

!*******************************************************************************
! La funcion Ham calcula el hamiltoniano para una dada configuración de spines
!*******************************************************************************
  function Ham(ss)
  implicit none
    integer, dimension(1:ns,1:ns), intent(in) :: ss
    real(pr)                                  :: Ham
    integer                                   :: n_c, n_f
    integer                                   :: sum_h
    sum_h = 0
    do n_c = 1, ns-1
      do n_f = 1, ns-1
        sum_h = sum_h + ss(n_f,n_c)*ss(n_f,n_c+1) + ss(n_f,n_c)*ss(n_f+1,n_c)
      end do
    end do
    do n_f = 1, ns
      sum_h = sum_h + ss(n_f,1)*ss(n_f,ns)
    end do
    do n_c = 1, ns
      sum_h = sum_h + ss(1,n_c)*ss(ns,n_c)
    end do
    do n_f = 1, ns-1
      sum_h = sum_h + ss(n_f,ns)*ss(n_f+1,ns)
    end do
    do n_c = 1, ns-1
      sum_h = sum_h + ss(ns,n_c)*ss(ns,n_c+1)
    end do
    Ham = -real(sum_h,pr)
  end function Ham

!*******************************************************************************
! La funcion mu calcula la magnetizacion particular de una dada configuracion
! de spines
!*******************************************************************************
  function mu(ss)
  implicit none
    integer, dimension(1:ns,1:ns), intent(in) :: ss
    real(pr)                                  :: mu
    integer                                   :: sum_m
    integer                                   :: n_f, n_c
    sum_m = 0
    do n_c = 1, ns
      do n_f = 1, ns
        sum_m = sum_m + ss(n_f,n_c)
      end do
    end do
    mu = real(sum_m,pr)
  end function mu

!*******************************************************************************
! La funcion delta_E calcula la variacion en la energia producida al flipear un
! spin en la fila n_f y columna n_c
!*******************************************************************************
  function delta_E(n_f,n_c,ss)
  implicit none
    integer, intent(in)                       :: n_f, n_c
    integer, dimension(1:ns,1:ns), intent(in) :: ss
    real(pr)                                  :: delta_E
    integer                                   :: sum_d, n_fi, n_ff, n_ci, n_cf
    if (n_f==1) then
      n_fi = ns
    else
      n_fi = n_f-1
    end if
    if (n_c==1) then
      n_ci = ns
    else
      n_ci = n_c-1
    end if
    if (n_f==ns) then
      n_ff = 1
    else
      n_ff = n_f+1
    end if
    if (n_c==ns) then
      n_cf = 1
    else
      n_cf = n_c+1
    end if
    sum_d = ss(n_f,n_c)*(ss(n_fi,n_c)+ss(n_ff,n_c)+ss(n_f,n_ci)+ss(n_f,n_cf))
    delta_E = real(2*sum_d,pr)
  end function delta_E

!*******************************************************************************
! La subrutina flip da vuelta el spin en la fila n_f y columna n_c
!*******************************************************************************
  subroutine flip(n_f,n_c,ss)
  implicit none
    integer, intent(in)                          :: n_f, n_c
    integer, dimension(1:ns,1:ns), intent(inout) :: ss
    integer                                      :: s_i, s_f
    s_i = ss(n_f,n_c)
    if (s_i/=1) then
      s_f = +1
    else if (s_i/=-1) then
      s_f = -1
    end if
    ss(n_f,n_c) = s_f
    return
  end subroutine flip

!*******************************************************************************
! La subrutina metropolis emplea el algoritmo de Metropolis para decidir si se
! flipea o no un spine
!*******************************************************************************
  subroutine metropolis(n_f,n_c,ss,Tem)
  implicit none
    integer, intent(in)                          :: n_f, n_c
    integer, dimension(1:ns,1:ns), intent(inout) :: ss
    real(pr), intent(in)                         :: Tem
    real(pr)                                     :: E_d, u, R
    E_d = delta_E(n_f,n_c,ss)
    if(E_d<=0._pr) then
      call flip(n_f,n_c,ss)
      return
    else
      R = exp(-E_d/Tem)
      u = grnd()
      if (R>=u) then
        call flip(n_f,n_c,ss)
        return
      else
        return
      end if
    end if
    return
  end subroutine metropolis

!*******************************************************************************
! La subrutina inicializar inicializa los spines del sistema de forma aleatoria
!*******************************************************************************
  subroutine inicializar(ss)
  implicit none
    integer, dimension(1:ns,1:ns), intent(inout) :: ss
    integer                                      :: n_f, n_c
    real(pr)                                     :: u
    ss = 0
    do n_f = 1, ns
      do n_c = 1, ns
        u = grnd()
        if (u<0.25_pr .or. 0.75<=u) then
          ss(n_f,n_c) = -1
        else
          ss(n_f,n_c) = +1
        end if
      end do
    end do
    return
  end subroutine inicializar

end program p1c
