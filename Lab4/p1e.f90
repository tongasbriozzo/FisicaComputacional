!*******************************************************************************
!
! Program p1e
!
! Este programa calcula la funcion de autocorrelacion en la energia y la
! magnetizacion en un modelo de Ising 2D de ns por ns spines con condiciones de
! contorno periodicas, en funcion de la temperatura
!
! Gaston Briozzo, FaMAFyC, UNC, 20/05/2021
!
!*******************************************************************************

program p1e
use precision, only: pr
use mtmod,     only: grnd
implicit none

  integer                            :: i, j, k
  integer                            :: p_mc, nf, nc
  integer, parameter                 :: ns=40, np_mc=1000, nsim=1000000, ndes=1000, nt=9
  integer, dimension(1:nt,1:ns,1:ns) :: s
  real(pr)                           :: V
  real(pr), dimension(1:nt)          :: T, E0, M0, E, M, E2, M2
  real(pr), dimension(1:nt,1:np_mc)  :: E0t, M0t, CE, CM, HE, HM

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
! Inicializamos las vartiables
  E(:) = 0._pr
  M(:) = 0._pr
  E0(:) = 0._pr
  M0(:) = 0._pr
  E2(:) = 0._pr
  M2(:) = 0._pr
  E0t(:,:) = 0._pr
  M0t(:,:) = 0._pr
  HE(:,:) = 0._pr
  HM(:,:) = 0._pr

! Tomamos un estado inicial al azar
  call inicializar(s)
! Descartamos los primeros ndes pasos
  do k = 1, ndes
    call pasoMC(s,T)
  end do
! Guardamos los vectores historia
  do p_mc = 1, np_mc
    call pasoMC(s,T)
    do j = 1, nt
      HE(j,p_mc) =    Ham(j,s)
      HM(j,p_mc) = abs(mu(j,s))
    end do
  end do
! promediamos sobre nsim simulaciones
  do i = 1, nsim
    write(*,*) 'simulacion numero', i
!  Tomamos los valores al tiempo inicial
    E0(:) = HE(:,1)
    M0(:) = HM(:,1)
!   Sumamos para obtener los valores de expectacion
    E(:)  = E(:)  +  HE(:,1)
    E2(:) = E2(:) + (HE(:,1))**2
    M(:)  = M(:)  +  HM(:,1)
    M2(:) = M2(:) + (HM(:,1))**2
!   Damos un paso Monte Carlo y actualizamos la historia
    call pasoMC(s,T)
    do p_mc = 1, np_mc
      HE(:,p_mc) = HE(:,p_mc+1)
      HM(:,p_mc) = HM(:,p_mc+1)
    end do
    do j = 1, nt
      HE(j,np_mc) =    Ham(j,s)
      HM(j,np_mc) = abs(mu(j,s))
    end do
!   Actualizada la historia, sumamos ara obtener los valores de expectascion
    do p_mc = 1, np_mc
      E0t(:,p_mc) = E0t(:,p_mc) + E0(:)*HE(:,p_mc)
      M0t(:,p_mc) = M0t(:,p_mc) + M0(:)*HM(:,p_mc)
    end do
  end do

! Normalizamos
  E(:)     = E(:)/real(nsim,pr)/V
  E2(:)    = E2(:)/real(nsim,pr)/(V*V)
  E0t(:,:) = E0t(:,:)/real(nsim,pr)/(V*V)
  M(:)     = M(:)/real(nsim,pr)/V
  M2(:)    = M2(:)/real(nsim,pr)/(V*V)
  M0t(:,:) = M0t(:,:)/real(nsim,pr)/(V*V)
! Abrimos los documentos de texto
  open(451,file='ac1.d')
  open(452,file='ac2.d')
  open(453,file='ac3.d')
  open(454,file='ac4.d')
  open(455,file='ac5.d')
  open(456,file='ac6.d')
  open(457,file='ac7.d')
  open(458,file='ac8.d')
  open(459,file='ac9.d')
! Anotamos los calores obtenidos
  do j = 1, nt                    ! El indice j recorre las temperaturas
    do p_mc = 1, np_mc            ! El indice p_mc recorre los pasos de Monte Carlo
      CE(j,p_mc) = (E0t(j,p_mc)-E(j)**2)/(E2(j)-E(j)**2)  ! Autocorrelacion en la Energia
      CM(j,p_mc) = (M0t(j,p_mc)-M(j)**2)/(M2(j)-M(j)**2)  ! Autocorrelacion en la Magnetizacion
      write(450+j,*) p_mc, CE(j,p_mc), CM(j,p_mc)
    end do
    close(450+j)
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

!*******************************************************************************
! La subrutina pasoMC aplica el algoritmo de Metropolis una vez a cada spin del
! sistema, en orden aleatorio
!*******************************************************************************
  subroutine pasoMC(ss,Tem)
  implicit none
    integer, dimension(1:nt,1:ns,1:ns), intent(inout) :: ss
    real(pr), dimension(1:nt), intent(in)             :: Tem
    integer                                           :: jj, ii, n_f, n_c, contador
    integer, dimension(1:ns,1:ns)                     :: sc
    do jj = 1, nt
      do n_f = 1, ns
        do n_c = 1, ns
          call metropolis(jj,n_f,n_c,ss,Tem)
        end do
      end do
    end do
!***************************************************************************************************
!    NO SE CONSIGUIO HACER FUNCIONAR UN GENERADOR DE PASOS QUE ACOMODE LOS SPINES EN ORDEN ALEATORIO
!***************************************************************************************************
!    sc(:,:) = 0
!    do ii = 1, ns*ns
!      n_f = nint(real(ns,pr)*grnd())
!      n_c = nint(real(ns,pr)*grnd())
!      write(*,*) n_f, n_c
!      do while (sc(n_f,n_c)==1)
!        n_f = mod(n_f,ns) + 1
!        contador = contador + 1
!        if (contador==ns) then
!          n_c = mod(n_c,ns) + 1
!          contador = 0
!        end if
!      end do
!      do jj = 1, nt
!        call metropolis(jj,n_f,n_c,ss,Tem)
!      end do
!      sc(n_f,n_c) = 1
!    end do
    return
  end subroutine pasoMC

end program p1e
