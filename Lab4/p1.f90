!*******************************************************************************
!
! Programa p1
!
! Este programa calcula la Energía, el calor especifico, la magnetizacion y la
! susceptibilidad magnetica de un sistema de ns por ns spines con condiciones
! de contorno periodicas, en funcion de los pasos Monte Carlo, para distintas
! temperaturas.
!
! Gaston Briozzo, FaMAFyC, UNC, 18/05/2020
!
!*******************************************************************************

program p1
use precision, only: pr
use mtmod,     only: grnd
implicit none

  integer                            :: i, j
  integer                            :: p_mc, nf, nc
  integer, parameter                 :: ns=40, np_mc=10000, nsim=1000, nt=9
  integer, dimension(1:nt,1:ns,1:ns) :: s
  real(pr)                           :: V
  real(pr), dimension(1:nt)          :: T
  real(pr), dimension(1:nt,1:np_mc)  :: E, M, C, X
  real(pr), dimension(1:nt,1:np_mc)  :: H2, mu2

! Calculamos el volumen del sistema
  V = real(ns*ns,pr)
! Definimos las temperaturas a evaluar
  T(1) = 2.0000_pr
  T(2) = 2.2200_pr
  T(3) = 2.2676_pr
  T(4) = 2.2682_pr
  T(5) = 2.2692_pr
  T(6) = 2.2702_pr
  T(7) = 2.2708_pr
  T(8) = 2.5000_pr
  T(9) = 3.3000_pr

!*******************************************************************************
! Condicion iniciali T0->0, con todos los spines alineados
!*******************************************************************************

! Inicializamos las variables
  E(:,:) = 0._pr
  M(:,:) = 0._pr
  C(:,:) = 0._pr
  X(:,:) = 0._pr
  H2(:,:)  = 0._pr
  mu2(:,:) = 0._pr
! Realizamos nsim simulaciones independientes
  do i = 1, nsim
!   Supervisamos el avance del programa
    write(*,*) 'alineado', i
    s = 1                                    ! Todos los spines up
    do p_mc = 1, np_mc                       ! p_mc indica el paso de Monte Carlo donde nos encontramos
      do nf = 1, ns                          ! nf indica el numero de fila donde nos encontramos
        do nc = 1, ns                        ! nc indica el numero de columna donde nos encontramos
          do j = 1, nt                       ! j indica la temperatura del sistema
            call metropolis(j,nf,nc,s,T)     ! Implementamos el algoritmo de Metropolis
          end do
        end do
      end do
!     Sumamos todos los valores para las nsim simulaciones
      do j = 1, nt ! nos movemos sobre todas las temperaturas
        E(j,p_mc)   = E(j,p_mc)   +  Ham(j,s)      ! <H>
        H2(j,p_mc)  = H2(j,p_mc)  + (Ham(j,s))**2  !<H^2>
        M(j,p_mc)   = M(j,p_mc)   + abs(mu(j,s))   ! <mu>
        mu2(j,p_mc) = mu2(j,p_mc) + (mu(j,s))**2   ! <mu^2>
      end do
    end do
  end do
! Normalizamos
  E(:,:)   = E(:,:)/real(nsim,pr)/V
  H2(:,:)  = H2(:,:)/real(nsim,pr)/(V*V)
  M(:,:)   = M(:,:)/real(nsim,pr)/V
  mu2(:,:) = mu2(:,:)/real(nsim,pr)/(V*V)
! Calculamos calor especifico y susceptibilidad magnetica
  do i = 1, nt ! nos movemos sobre todas las temperaturas
    C(i,:) = (H2(i,:)  - E(i,:)**2)*V/(T(i)*T(i))
    X(i,:) = (mu2(i,:) - M(i,:)**2)*V/T(i)
  end do
! Abrimos los documentos de texto
  open(411,file='sa1.d')
  open(412,file='sa2.d')
  open(413,file='sa3.d')
  open(414,file='sa4.d')
  open(415,file='sa5.d')
  open(416,file='sa6.d')
  open(417,file='sa7.d')
  open(418,file='sa8.d')
  open(419,file='sa9.d')
! Anotamos los valores obtenidos
  do j = 1, nt
    do p_mc = 1, np_mc
      write(410+j,*) p_mc, E(j,p_mc), C(j,p_mc), M(j,p_mc), X(j,p_mc)
    end do
    close(410+j)
  end do


!*******************************************************************************
! Condicion iniciali T0->inf, con todos los spines desalineados
!*******************************************************************************

! Identico al desarrollo anterior
  E(:,:) = 0._pr
  M(:,:) = 0._pr
  C(:,:) = 0._pr
  X(:,:) = 0._pr
  H2(:,:)  = 0._pr
  mu2(:,:) = 0._pr

  do i = 1, nsim
    write(*,*) 'desalineado', i
!   La subrutina inicializar emplea el generador de numeros pseudoaleatorios
!   mtmod para generar un estado inicial al azar
    call inicializar(s)
    do p_mc = 1, np_mc
      do nf = 1, ns
        do nc = 1, ns
          do j = 1, nt
            call metropolis(j,nf,nc,s,T)
          end do
        end do
      end do
      do j = 1, nt
        E(j,p_mc)   = E(j,p_mc)   +  Ham(j,s)
        H2(j,p_mc)  = H2(j,p_mc)  + (Ham(j,s))**2
        M(j,p_mc)   = M(j,p_mc)   + abs(mu(j,s))
        mu2(j,p_mc) = mu2(j,p_mc) + (mu(j,s))**2
      end do
    end do
  end do

  E(:,:)   = E(:,:)/real(nsim,pr)/V
  H2(:,:)  = H2(:,:)/real(nsim,pr)/(V*V)
  M(:,:)   = M(:,:)/real(nsim,pr)/V
  mu2(:,:) = mu2(:,:)/real(nsim,pr)/(V*V)

  do i = 1, nt
    C(i,:) = (H2(i,:)  - E(i,:)**2)*V/(T(i)*T(i))
    X(i,:) = (mu2(i,:) - M(i,:)**2)*V/T(i)
  end do

  open(421,file='sd1.d')
  open(422,file='sd2.d')
  open(423,file='sd3.d')
  open(424,file='sd4.d')
  open(425,file='sd5.d')
  open(426,file='sd6.d')
  open(427,file='sd7.d')
  open(428,file='sd8.d')
  open(429,file='sd9.d')

do j = 1, nt
  do p_mc = 1, np_mc
    write(420+j,*) p_mc, E(j,p_mc), C(j,p_mc), M(j,p_mc), X(j,p_mc)
  end do
  close(420+j)
end do

CONTAINS

!*******************************************************************************
! La funcion Ham calcula el hamiltoniano para una dada configuración de spines
!*******************************************************************************
  function Ham(jj,ss)
  implicit none
    integer, intent(in)                            :: jj
    integer, dimension(1:nt,1:ns,1:ns), intent(in) :: ss
    real(pr)                                       :: Ham
    integer                                        :: n_c, n_f
    integer                                        :: sum_h
    sum_h = 0
    do n_c = 1, ns-1
      do n_f = 1, ns-1
        sum_h = sum_h + ss(jj,n_f,n_c)*ss(jj,n_f,n_c+1) + ss(jj,n_f,n_c)*ss(jj,n_f+1,n_c)
      end do
    end do
    do n_f = 1, ns
      sum_h = sum_h + ss(jj,n_f,1)*ss(jj,n_f,ns)
    end do
    do n_c = 1, ns
      sum_h = sum_h + ss(jj,1,n_c)*ss(jj,ns,n_c)
    end do
    do n_f = 1, ns-1
      sum_h = sum_h + ss(jj,n_f,ns)*ss(jj,n_f+1,ns)
    end do
    do n_c = 1, ns-1
      sum_h = sum_h + ss(jj,ns,n_c)*ss(jj,ns,n_c+1)
    end do
    Ham = -real(sum_h,pr)
  end function Ham

!*******************************************************************************
! La funcion mu calcula la magnetizacion particular de una dada configuracion
! de spines
!*******************************************************************************
  function mu(jj,ss)
  implicit none
    integer, intent(in)                            :: jj
    integer, dimension(1:nt,1:ns,1:ns), intent(in) :: ss
    real(pr)                                       :: mu
    integer                                        :: sum_m
    integer                                        :: n_f, n_c
    sum_m = 0
    do n_c = 1, ns
      do n_f = 1, ns
        sum_m = sum_m + ss(jj,n_f,n_c)
      end do
    end do
    mu = real(sum_m,pr)
  end function mu

!*******************************************************************************
! La funcion delta_E calcula la variacion en la energia producida al flipear un
! spin en la fila n_f y columna n_c
!*******************************************************************************
  function delta_E(jj,n_f,n_c,ss)
  implicit none
    integer, intent(in)                            :: jj, n_f, n_c
    integer, dimension(1:nt,1:ns,1:ns), intent(in) :: ss
    real(pr)                                       :: delta_E
    integer                                        :: sum_d, n_fi, n_ff, n_ci, n_cf
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
    sum_d = ss(jj,n_f,n_c)*(ss(jj,n_fi,n_c)+ss(jj,n_ff,n_c)+ss(jj,n_f,n_ci)+ss(jj,n_f,n_cf))
    delta_E = real(2*sum_d,pr)
  end function delta_E

  subroutine flip(jj,n_f,n_c,ss)
  implicit none
    integer, intent(in)                               :: jj, n_f, n_c
    integer, dimension(1:nt,1:ns,1:ns), intent(inout) :: ss
    integer                                           :: s_i, s_f
    s_i = ss(jj,n_f,n_c)
    if (s_i/=1) then
      s_f = +1
    else if (s_i/=-1) then
      s_f = -1
    end if
    ss(jj,n_f,n_c) = s_f
    return
  end subroutine flip

  subroutine metropolis(jj,n_f,n_c,ss,Tem)
  implicit none
    integer, intent(in)                               :: n_f, n_c, jj
    integer, dimension(1:nt,1:ns,1:ns), intent(inout) :: ss
    real(pr), dimension(1:nt), intent(in)             :: Tem
    real(pr)                                          :: E_d, u, R
    E_d = delta_E(jj,n_f,n_c,ss)
    if(E_d<=0._pr) then
      call flip(jj,n_f,n_c,ss)
      return
    else
      R = exp(-E_d/Tem(jj))
      u = grnd()
      if (R>=u) then
        call flip(jj,n_f,n_c,ss)
        return
      else
        return
      end if
    end if
    return
  end subroutine metropolis

  subroutine inicializar(ss)
  implicit none
    integer, dimension(1:nt,1:ns,1:ns), intent(inout) :: ss
    integer                                           :: jj, n_f, n_c
    real(pr)                                          :: u
    ss = 0
    do jj = 1, nt
      do n_f = 1, ns
        do n_c = 1, ns
          u = grnd()
          if (u<0.25_pr .or. 0.75<=u) then
            ss(jj,n_f,n_c) = -1
          else
            ss(jj,n_f,n_c) = +1
          end if
        end do
      end do
    end do
    return
  end subroutine inicializar

end program p1
