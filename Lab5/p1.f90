!*******************************************************************************
!
!  Dinamica Molecular
!
!  Este programa emplea el metodo velocity-Verlet para describir un sistema de N
!  atomos interactuantes segun el potencial de Lennard-Jones con condiciones
!  periodicas de contorno, con densidad rho, temperatura T y paso temporal dt.
!
!  Gaston Briozzo, FaMAFyC, UNC, 04/06/2021
!
!*******************************************************************************

program p1
use precision, only: pr
use mtmod,     only: grnd
implicit none

  integer                         :: i, j, np, nh, ns
  integer, parameter              :: N=256, t_eq=1000, t_scal=50, t_run=1000, nsim=1000
  real(pr)                        :: T
  real(pr), dimension(1:2000)     ::  Ti,   Ec,  Ep,  Et,  P
  real(pr)                        :: iTi,  iEc, iEp, iEt, iP
  real(pr)                        :: dTi,  dEc, dEp, dEt, dP
  real(pr)                        ::  Tim,  Ecm, Epm, Etm, Pm
  real(pr)                        ::  Ti2,  Ec2, Ep2, Et2, P2
  real(pr), dimension(1:N,1:3)    :: x, v, xnew, vnew
  real(pr), dimension(-50:49,1:3) :: his, hisi
  real(pr), parameter             :: r_cut = 2.5_pr, r2cut=r_cut*r_cut, r_inv=1._pr/r_cut, r_i3=r_inv*r_inv*r_inv
  real(pr), parameter             :: rho = 0.8_pr, Vol=real(N,pr)/rho, L=(Vol)**(1._pr/3._pr)
  real(pr), parameter             :: dt=0.005_pr
  real(pr), parameter             :: pi=acos(-1._pr)
  real(pr), parameter             :: U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  real(pr), parameter             :: P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr

! Inicializamos el sistema
  T = 1.1_pr        ! Temperatura

! Inicializamos los valores instantaneos
  Ti(:) = 0._pr
  Ec(:) = 0._pr
  Ep(:) = 0._pr
  Et(:) = 0._pr
  P(:)  = 0._pr
! Inicializamos las sumatorias
  Tim = 0._pr
  Ti2 = 0._pr
  Ecm = 0._pr
  Ec2 = 0._pr
  Epm = 0._pr
  Ep2 = 0._pr
  Etm = 0._pr
  Et2 = 0._pr
  Pm  = 0._pr
  P2  = 0._pr
! inicializamos el histograma
  his(:,:) = 0._pr

DO ns = 1, nsim          ! Promediamos sobre nsim simulaciones independientes

  write(*,*) 'nsim', ns  ! Supervisamos el funcionamiento del sistema

  call FCC(x)       ! Acomodamos las particulas segùn una FCC
  call ini_v(T,v)   ! Inicializamos las velocidades para que la temperatura resulte T

!*******************************************************************************
! Problema 1a, primeros mil pasos
!*******************************************************************************

  do i = 1, t_eq                      ! Nos movemos hasta el paso t_eq
    call verlet_v(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)                ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
    v(:,:) = vnew(:,:)                ! Por eso hacemos el remplazo manual
    if (mod(i,t_scal)==0) then        ! Cada t_scal pasos, rescaleamos la velocidad para conservar la temperatura (Baño termico)
!      call rescal_v(T,v)             ! La subrutina rescal_v hace exactamente lo mismo que lo de la linea de abajo
      iTi = Temperatura(v)
      v(:,:) = v(:,:)*sqrt(T/iTi)
    end if
!   Calculamos los valores instantaneos de las variables de estado
    iTi = Temperatura(v)
    iEc = K(v)
    iEp = U(x) + N*U_tail
    iEt = iEc + iEp
    iP  = rho*iTi + pf(x)/(3._pr*Vol) + P_tail
!   Promediamos sobre las nsim simulaciones
    Ti(i) = Ti(i) + iTi
    Ec(i) = Ec(i) + iEc
    Ep(i) = Ep(i) + iEp
    Et(i) = Et(i) + iEt
    P(i)  = P(i)  + iP
  end do

!*******************************************************************************
! Problema 1a, segundos mil pasos
! Problema 1b
! Problema 1c
!*******************************************************************************

  do i = t_eq + 1, t_eq+t_run                           ! Corremos sobre t_run pasos
    call verlet_v(x,v,xnew,vnew)
    x(:,:) = xnew(:,:)
    v(:,:) = vnew(:,:)
!   Calculamos los valores instantaneos de las variables de estado
    iTi = Temperatura(v)
    iEc = K(v)
    iEp = U(x) + N*U_tail
    iEt = iEc + iEp
    iP  = rho*iTi + pf(x)/(3._pr*Vol) + P_tail
!   Promediamos sobre las nsim simulaciones
    Ti(i) = Ti(i) + iTi
    Ec(i) = Ec(i) + iEc
    Ep(i) = Ep(i) + iEp
    Et(i) = Et(i) + iEt
    P(i)  = P(i)  + iP
!   Tomamos los valores medios
    Tim = Tim + iTi
    Ti2 = Ti2 + iTi*iTi
    Ecm = Ecm + iEc
    Ec2 = Ec2 + iEc*iEc
    Epm = Epm + iEp
    Ep2 = Ep2 + iEp*iEp
    Etm = Etm + iEt
    Et2 = Et2 + iEt*iEt
    Pm  = Pm  + iP
    P2  = P2  + iP*iP
!   Ahora examinamos las velocidades para construir el histograma del problema 1c
    do np = 1, N
      do j = 1, 3
        do nh = -50, 49
          if (real(nh)<=10._pr*v(np,j).and.10._pr*v(np,j)<real(nh+1,pr)) then
            his(nh,j) = his(nh,j) + 1._pr
            cycle
          end if
        end do
      end do
    end do
  end do

END DO

! Normalizamos los valores instantaneos y los anotamos
  Ti(:) = Ti(:)/real(nsim,pr)
  Ec(:) = Ec(:)/real(nsim,pr)
  Ep(:) = Ep(:)/real(nsim,pr)
  Et(:) = Et(:)/real(nsim,pr)
  P(:)  = P(:)/real(nsim,pr)
  open(511,file='1a.d')
  do i = 1, t_eq+t_run
    write(511,*) i, Ti(i), Ec(i), Ep(i), Et(i), P(i)
  end do
  close(511)

! Normalizamos el histograma y guardamos sus valores concluyendo con el problema 1c
  his(:,:) = his(:,:)/real(N,pr)/real(t_run,pr)/real(nsim,pr)
  open(513,file='his.d')
  do nh = -50, 49
    write(513,*) 0.1_pr*real(nh,pr),   0._pr,     0._pr,     0._pr
    write(513,*) 0.1_pr*real(nh,pr),   his(nh,1), his(nh,2), his(nh,3)
    write(513,*) 0.1_pr*real(nh+1,pr), his(nh,1), his(nh,2), his(nh,3)
    write(513,*) 0.1_pr*real(nh+1,pr), 0._pr,     0._pr,     0._pr
  end do
  close(513)

! Normalizamos los valores medios
  Tim = Tim/real(t_run,pr)/real(nsim,pr)
  Ti2 = Ti2/real(t_run,pr)/real(nsim,pr)
  Ecm = Ecm/real(t_run,pr)/real(nsim,pr)
  Ec2 = Ec2/real(t_run,pr)/real(nsim,pr)
  Epm = Epm/real(t_run,pr)/real(nsim,pr)
  Ep2 = Ep2/real(t_run,pr)/real(nsim,pr)
  Etm = Etm/real(t_run,pr)/real(nsim,pr)
  Et2 = Et2/real(t_run,pr)/real(nsim,pr)
  Pm  = Pm /real(t_run,pr)/real(nsim,pr)
  P2  = P2 /real(t_run,pr)/real(nsim,pr)
! Calculamos las varianzas
  dTi = sqrt(Ti2-Tim*Tim)
  dEc = sqrt(Ec2-Ecm*Ecm)
  dEp = sqrt(Ep2-Epm*Epm)
  dEt = sqrt(Et2-Etm*Etm)
  dP  = sqrt( P2- Pm* Pm)
! Anotamos valores medios y varianzas concluyeno con el problema 1b
  open(512,file='1b.d')
  write(512,*) Tim, Ecm, Epm, Etm, Pm
  write(512,*) dTi, dEc, dEp, dEt, dP
  close(512)

!*******************************************************************************
! Problema 1d
! El procedimiento es identico al realizado para los problemas 1 a, b y c, solo
! que partiendo ahora de una distribucion de velocidades de Maxwell-Boltzmann
!*******************************************************************************

! Inicializamos los valores instantaneos
  Ti(:) = 0._pr
  Ec(:) = 0._pr
  Ep(:) = 0._pr
  Et(:) = 0._pr
  P(:)  = 0._pr
! Inicializamos las sumatorias
  Tim = 0._pr
  Ti2 = 0._pr
  Ecm = 0._pr
  Ec2 = 0._pr
  Epm = 0._pr
  Ep2 = 0._pr
  Etm = 0._pr
  Et2 = 0._pr
  Pm  = 0._pr
  P2  = 0._pr
! inicializamos los histogramas
  his(:,:)  = 0._pr
  hisi(:,:) = 0._pr

DO ns = 1, nsim          ! Promediamos sobre nsim simulaciones independientes

  write(*,*) 'nsim', ns  ! Supervisamos el funcionamiento del sistema

! Inicializamos el sistema
  call FCC(x)         ! Acomodamos las particulas segùn una FCC
! Elegimos las componentes de la velocidad segun una distribucion gaussiana
  do j = 1, 3
    do i = 1, N
      v(i,j) = gasdev()
    end do
  end do
  call vcm(v)
! Calculamos la temperatura cinetica del sistema y renormalizamos las velocidades para obtener la temperatura deseada
  iTi = Temperatura(v)
  v(:,:) = v(:,:)*sqrt(T/iTi)

! Construimos un histograma con la distribucion inicial de velocidades para comprobar su naturaleza gaussiana
  do np = 1, N
    do j = 1, 3
      do nh = -50, 49
        if (real(nh)<=10._pr*v(np,j).and.10._pr*v(np,j)<real(nh+1,pr)) then
          hisi(nh,j) = hisi(nh,j) + 1._pr
          cycle
        end if
      end do
    end do
  end do

! Ahora, realizaremos los primeros t_eq pasos hasta termalizar
  do i = 1, t_eq                      ! Nos movemos hasta el paso t_eq
    call verlet_v(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)                ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
    v(:,:) = vnew(:,:)                ! Por eso hacemos el remplazo manual
    if (mod(i,t_scal)==0) then        ! Cada t_scal pasos, rescaleamos la velocidad para conservar la temperatura (Baño termico)
!      call rescal_v(T,v)             ! La subrutina rescal_v hace exactamente lo mismo que lo de la linea de abajo
      iTi = Temperatura(v)
      v(:,:) = v(:,:)*sqrt(T/iTi)
    end if
!   Calculamos los valores instantaneos de las variables de estado
    iTi = Temperatura(v)
    iEc = K(v)
    iEp = U(x) + N*U_tail
    iEt = iEc + iEp
    iP  = rho*iTi + pf(x)/(3._pr*Vol) + P_tail
!   Promediamos sobre las nsim simulaciones
    Ti(i) = Ti(i) + iTi
    Ec(i) = Ec(i) + iEc
    Ep(i) = Ep(i) + iEp
    Et(i) = Et(i) + iEt
    P(i)  = P(i)  + iP
  end do
! Ahora realizaremos los siguientes t_run pasos, tomando valores medios y varianzas
  do i = t_eq+1, t_eq+t_run                           ! Corremos sobre t_run pasos
    call verlet_v(x,v,xnew,vnew)
    x(:,:) = xnew(:,:)
    v(:,:) = vnew(:,:)
!   Calculamos los valores instantaneos de las variables de estado
    iTi = Temperatura(v)
    iEc = K(v)
    iEp = U(x) + N*U_tail
    iEt = iEc + iEp
    iP  = rho*iTi + pf(x)/(3._pr*Vol) + P_tail
!   Promediamos sobre las nsim simulaciones
    Ti(i) = Ti(i) + iTi
    Ec(i) = Ec(i) + iEc
    Ep(i) = Ep(i) + iEp
    Et(i) = Et(i) + iEt
    P(i)  = P(i)  + iP
!   Tomamos los valores medios
    Tim = Tim + iTi
    Ti2 = Ti2 + iTi*iTi
    Ecm = Ecm + iEc
    Ec2 = Ec2 + iEc*iEc
    Epm = Epm + iEp
    Ep2 = Ep2 + iEp*iEp
    Etm = Etm + iEt
    Et2 = Et2 + iEt*iEt
    Pm  = Pm  + iP
    P2  = P2  + iP*iP
!   Ahora examinamos las velocidades para construir el histograma
    do np = 1, N
      do j = 1, 3
        do nh = -50, 49
          if (real(nh)<=10._pr*v(np,j).and.10._pr*v(np,j)<real(nh+1,pr)) then
            his(nh,j) = his(nh,j) + 1._pr
            cycle
          end if
        end do
      end do
    end do
  end do
  close(515)

END DO

! Normalizamos los valores instantaneos y los anotamos
  Ti(:) = Ti(:)/real(nsim,pr)
  Ec(:) = Ec(:)/real(nsim,pr)
  Ep(:) = Ep(:)/real(nsim,pr)
  Et(:) = Et(:)/real(nsim,pr)
  P(:)  = P(:)/real(nsim,pr)
  open(515,file='igt.d')              !
  do i = 1, t_eq+t_run
    write(515,*) i, Ti(i), Ec(i), Ep(i), Et(i), P(i)
  end do
  close(515)

! Normalizamos el histograma de velocidades iniciales y guardamos los valores obtenidos
  hisi(:,:) = hisi(:,:)/real(N,pr)/real(nsim,pr)
  open(514,file='hg1.d')
  do nh = -50, 49
    write(514,*) 0.1_pr*real(nh,pr),   0._pr,     0._pr,     0._pr
    write(514,*) 0.1_pr*real(nh,pr),   hisi(nh,1), hisi(nh,2), hisi(nh,3)
    write(514,*) 0.1_pr*real(nh+1,pr), hisi(nh,1), hisi(nh,2), hisi(nh,3)
    write(514,*) 0.1_pr*real(nh+1,pr), 0._pr,     0._pr,     0._pr
  end do
  close(514)

! Normalizamos el histograma de velocidades y guardamos sus valores
  his(:,:) = his(:,:)/real(N,pr)/real(t_run,pr)/real(nsim,pr)
  open(516,file='hgt.d')
  do nh = -50, 49
    write(516,*) 0.1_pr*real(nh,pr),   0._pr,     0._pr,     0._pr
    write(516,*) 0.1_pr*real(nh,pr),   his(nh,1), his(nh,2), his(nh,3)
    write(516,*) 0.1_pr*real(nh+1,pr), his(nh,1), his(nh,2), his(nh,3)
    write(516,*) 0.1_pr*real(nh+1,pr), 0._pr,     0._pr,     0._pr
  end do
  close(516)

! Normalizamos los valores medios
  Tim = Tim/real(t_run,pr)/real(nsim,pr)
  Ti2 = Ti2/real(t_run,pr)/real(nsim,pr)
  Ecm = Ecm/real(t_run,pr)/real(nsim,pr)
  Ec2 = Ec2/real(t_run,pr)/real(nsim,pr)
  Epm = Epm/real(t_run,pr)/real(nsim,pr)
  Ep2 = Ep2/real(t_run,pr)/real(nsim,pr)
  Etm = Etm/real(t_run,pr)/real(nsim,pr)
  Et2 = Et2/real(t_run,pr)/real(nsim,pr)
  Pm  = Pm /real(t_run,pr)/real(nsim,pr)
  P2  = P2 /real(t_run,pr)/real(nsim,pr)
! Calculamos las varianzas
  dTi = sqrt(Ti2-Tim*Tim)
  dEc = sqrt(Ec2-Ecm*Ecm)
  dEp = sqrt(Ep2-Epm*Epm)
  dEt = sqrt(Et2-Etm*Etm)
  dP  = sqrt( P2 -Pm* Pm)
! Anotamos valores medios y varianzas
  open(517,file='vgt.d')
  write(517,*) Tim, Ecm, Epm, Etm, Pm
  write(517,*) dTi, dEc, dEp, dEt, dP
  close(517)

!*******************************************************************************
!*******************************************************************************
CONTAINS
!*******************************************************************************
!*******************************************************************************

!*******************************************************************************
! La funcion u(xx) nos da la energia potencial de un sistema de N particulas
! distribuidas segun las posiciones xx bajo un potencial de Lennard-Jones
!*******************************************************************************
  function u(xx)
  implicit none
    integer                                  :: ii, jj
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr)                                 :: u, sum_u
    real(pr)                                 :: dx, dy, dz, r2
    sum_u = 0._pr
    do ii = 2, N
      do jj = 1, ii-1
        dx = xx(jj,1) - xx(ii,1)
        dy = xx(jj,2) - xx(ii,2)
        dz = xx(jj,3) - xx(ii,3)
        dx = dx - L*anint(dx/L)
        dy = dy - L*anint(dy/L)
        dz = dz - L*anint(dz/L)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2==0._pr) then
!          write(*,*) 'Hubo distancia cero entre las particulas', ii, jj
          cycle
        else if (r2>=r2cut) then
          cycle
        end if
        sum_u = sum_u + 1._pr/(r2*r2*r2*r2*r2*r2) - 1._pr/(r2*r2*r2)
      end do
    end do
    u = 4._pr*sum_u
  end function u

!*******************************************************************************
! La funcion K(vv) nos da la energia cinetica de un sistema de N particulas con
! velocidades dadas por el vector vv
!*******************************************************************************
  function K(vv)
  implicit none
    integer                                  :: ii
    real(pr), dimension(1:N,1:3), intent(in) :: vv
    real(pr)                                 :: K, sum_k
    sum_k = 0._pr
    do ii = 1, N
      sum_k = sum_k + vv(ii,1)*vv(ii,1) + vv(ii,2)*vv(ii,2) + vv(ii,3)*vv(ii,3)
    end do
    K = 0.5_pr*sum_k
  end function K

!*******************************************************************************
! La funcion fij nos da la fuerza ejercida sobre la particula ii por la
! particula jj, con posiciones dadas por el vector xx, bajo un potencial de
! Lennard-Jones
!*******************************************************************************
  function fij(ii,jj,xx)
  implicit none
    integer, intent(in)                      :: ii, jj
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr), dimension    (1:3)             :: fij
    real(pr)                                 :: dx, dy, dz, r2
    if (jj==ii) then
      fij = 0._pr
      return
    end if
    dx = xx(ii,1) - xx(jj,1)
    dy = xx(ii,2) - xx(jj,2)
    dz = xx(ii,3) - xx(jj,3)
    dx = dx - L*anint(dx/L)
    dy = dy - L*anint(dy/L)
    dz = dz - L*anint(dz/L)
    r2 = dx*dx + dy*dy + dz*dz
    if (r2==0._pr) then
!      write(*,*) 'Hubo distancia cero entre las particulas', ii, jj
      fij = 0._pr
      return
    else if (r2>=r2cut) then
      fij = 0._pr
      return
    end if
    fij(1) = 48._pr*dx*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fij(2) = 48._pr*dy*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fij(3) = 48._pr*dz*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    return
  end function fij

!*******************************************************************************
! La funcion fi nos da la fuerza total ejercida sobre la particula ii en un
! sistema de N particulas con posiciones dadas por xx bajo un potencial de
! Lennard-Jones
!*******************************************************************************
  function fi(ii,xx)
  implicit none
    integer, intent(in)                      :: ii
    integer                                  :: jj
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr), dimension    (1:3)             :: fi, sum_fi
    sum_fi(:) = 0._pr
    do jj = 1, N
      sum_fi = sum_fi + fij(ii,jj,xx)
    end do
    fi = sum_fi
  end function fi

!*******************************************************************************
! La funcion pf da la suma de los productos escalares entre las fuerzas que se
! ejercen las particulas ii y jj y la distancia entre estas
!*******************************************************************************
  function pf(xx)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr)                                 :: pf, sum_pf
    integer                                  :: ii, jj
    real(pr), dimension(1:3)                 :: ff
    real(pr)                                 :: dx, dy, dz
    sum_pf = 0._pr
    do ii = 2, N
      do jj = 1, ii-1
        dx = xx(ii,1) - xx(jj,1)
        dy = xx(ii,2) - xx(jj,2)
        dz = xx(ii,3) - xx(jj,3)
        dx = dx - L*anint(dx/L)
        dy = dy - L*anint(dy/L)
        dz = dz - L*anint(dz/L)
        ff = fij(ii,jj,xx)
        sum_pf = sum_pf + dx*ff(1) + dy*ff(2) + dz*ff(3)
      end do
    end do
    pf = sum_pf
  end function pf

!*******************************************************************************
! La funcion Temperatura nos da la temperatura cinetica de un sistema de N
! particulas con velocidades dadas por el vector vv
!*******************************************************************************
  function Temperatura(vv)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in) :: vv
    real(pr)                                 :: Temperatura
    integer                                  :: ii
    real(pr)                                 :: sum_t
    sum_t = 0._pr
    do ii = 1, N
      sum_t = sum_t + vv(ii,1)*vv(ii,1) + vv(ii,2)*vv(ii,2) + vv(ii,3)*vv(ii,3)
    end do
    Temperatura = sum_t/real(3*N,pr)
  end function Temperatura

!*******************************************************************************
! La subrutina vcm toma la el conjunto de velocidades de un sistema de N
! particulas y modifica una al azar de modo que la velocidad resultante
! del centro de masa sea nula
!*******************************************************************************
  subroutine vcm(vv)
  implicit none
      real(pr), dimension(1:N,1:3), intent(inout) :: vv
      integer                                     :: ii, nnvv
      real(pr), dimension    (1:3)                :: sum_vv
      sum_vv(:) = 0._pr
      do ii = 1, N
        sum_vv(:) = sum_vv(:) + vv(ii,:)
      end do
      nnvv = 1 + int(real(N-1,pr)*grnd())
      vv(nnvv,:) = vv(nnvv,:) - sum_vv(:)
      return
    end subroutine

!*******************************************************************************
! La subrutina pbc aplica condiciones periodicas de contorno a un sistema de n
! particulas en una caja cubica de lado L
!*******************************************************************************
  subroutine pbc(xx)
  implicit none
    real(pr), dimension(1:N,1:3), intent(inout) :: xx
    integer                                     :: ii, jj
    do ii = 1, N
      do jj = 1, 3
        if (x(ii,jj)<0._pr.or.L<=x(ii,jj)) then
          xx(ii,jj) = xx(ii,jj) - L*floor(xx(ii,jj)/L)
        end if
      end do
    end do
    return
  end subroutine pbc

!*******************************************************************************
! La subrutina verlet_v implementa en metodo Velocity-Verlet para obtener las
! posiciones y velocidades de un sistema dinamico de N particulas en el instante
! siguiente. Esta diseñada de forma tal que tanto la posicion como la velocidad
! del centro de masa del sistema resulten nulas
!*******************************************************************************
  subroutine verlet_v(xx,vv,xxnew,vvnew)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in)  :: xx, vv
    real(pr), dimension(1:N,1:3), intent(out) :: xxnew, vvnew
    real(pr), dimension    (1:3)              :: xn, vn
    integer                                   :: ii, jj, nale
    xn(:) = 0._pr
    vn(:) = 0._pr
    do ii = 1, N
      xxnew(ii,:) = xx(ii,:) + vv(ii,:)*dt + 0.5_pr*fi(ii,xx)*dt*dt
      xn(:) = xn(:) + xxnew(ii,:)
    end do
    nale = 1 + anint(real(N-1,pr)*grnd())  ! Elegimos una particula al azar
    xxnew(nale,:) = xxnew(nale,:) - xn(:)  ! Aseguramos posicion del centro de masa nula
    call pbc(xxnew)                        ! Incluimos las condiciones periodicas de contorno
    do ii = 1, N
      vvnew(ii,:) = vv(ii,:) + 0.5_pr*(fi(ii,xxnew)+fi(ii,xx))*dt
      vn(:) = vn(:) + vvnew(ii,:)
    end do
    nale = 1 + anint(real(N-1,pr)*grnd())  ! Elegimos una particula al azar
    vvnew(nale,:) = vvnew(nale,:) - vn(:)  ! Aseguramos velocidad del centro de masa nula
    return
  end subroutine verlet_v

!*******************************************************************************
! La subrutina rescal_v calcula la temperatura cinetica de un sistema de N
! particulas con velocidades dadas por vv y las rescalea para obtener la
! temperatura deseada
!*******************************************************************************
  subroutine rescal_v(TT,vv)
  implicit none
    real(pr), intent(in)                        :: TT
    real(pr), dimension(1:N,1:3), intent(inout) :: vv
    integer                                     :: ii
    real(pr)                                    :: Tv
    Tv = 0._pr
    do ii = 1, N
      Tv = Tv + vv(ii,1)*vv(ii,1) + vv(ii,2)*vv(ii,2) + vv(ii,3)*vv(ii,3)
    end do
    Tv = Tv/(3._pr*real(N,pr))
    vv(:,:) = vv(:,:)*sqrt(TT/Tv)
    return
  end subroutine rescal_v

!*******************************************************************************
! La subrutina FCC acomoda las posiciones de un sistema de N particulas para
! obtener una red cubica centrada en las caras (FCC)
! Profe mire nomas que beieza de red
!*******************************************************************************
  subroutine FCC(xx)
  implicit none
    real(pr), dimension(1:N,1:3), intent(inout) :: xx
    integer                                     :: ii, jj, kk
    do ii = 1, 8, 2
      do jj = 1, 8, 2
        do kk = 1, 8, 2
          xx(1+2*(kk-1)+8*(jj-1)+32*(ii-1),1) = real(ii,pr)*L/8._pr
          xx(1+2*(kk-1)+8*(jj-1)+32*(ii-1),2) = real(jj,pr)*L/8._pr
          xx(1+2*(kk-1)+8*(jj-1)+32*(ii-1),3) = real(kk,pr)*L/8._pr
          xx(2+2*(kk-1)+8*(jj-1)+32*(ii-1),1) = real(ii,pr)*L/8._pr
          xx(2+2*(kk-1)+8*(jj-1)+32*(ii-1),2) = real(jj+1,pr)*L/8._pr
          xx(2+2*(kk-1)+8*(jj-1)+32*(ii-1),3) = real(kk+1,pr)*L/8._pr
          xx(3+2*(kk-1)+8*(jj-1)+32*(ii-1),1) = real(ii+1,pr)*L/8._pr
          xx(3+2*(kk-1)+8*(jj-1)+32*(ii-1),2) = real(jj,pr)*L/8._pr
          xx(3+2*(kk-1)+8*(jj-1)+32*(ii-1),3) = real(kk+1,pr)*L/8._pr
          xx(4+2*(kk-1)+8*(jj-1)+32*(ii-1),1) = real(ii+1,pr)*L/8._pr
          xx(4+2*(kk-1)+8*(jj-1)+32*(ii-1),2) = real(jj+1,pr)*L/8._pr
          xx(4+2*(kk-1)+8*(jj-1)+32*(ii-1),3) = real(kk,pr)*L/8._pr
        end do
      end do
    end do
    return
  end subroutine FCC

!*******************************************************************************
! La subrutina ini_v inicializa las velocidades de un sistema de N particulas
! empleando numeros al azar para obtener sus modulos y orientaciones. El
! resultado es reescaleato para obtener la temperatura cinetica deseada
!*******************************************************************************
  subroutine ini_v(TT,vv)
  implicit none
    real(pr), intent(in)                        :: TT
    real(pr), dimension(1:N,1:3), intent(inout) :: vv
    integer                                     :: ii, nnvv
    real(pr)                                    :: ang1, ang2, modulo_v, Tv
    real(pr), dimension    (1:3)                :: sum_vv
    sum_vv(:) = 0._pr
    do ii = 1, N
      modulo_v = grnd()
      ang1 = grnd()*2._pr*pi
      ang2 = grnd()*2._pr*pi
      vv(ii,1) = modulo_v*sin(ang1)*cos(ang2)
      vv(ii,2) = modulo_v*sin(ang1)*sin(ang2)
      vv(ii,3) = modulo_v*cos(ang1)
      sum_vv(:) = sum_vv(:) + vv(ii,:)   ! Esta es la velocidad del centro de masa
    end do
    nnvv = 1 + int(real(N-1)*grnd())
    vv(nnvv,:) = vv(nnvv,:) - sum_vv(:)  ! vcm = 0
    Tv = Temperatura(vv)
    vv(:,:) = vv(:,:)*sqrt(TT/Tv)        ! Renormalizamos la distribucion
    return
  end subroutine ini_v

!*******************************************************************************
! La funcion gasdev, extraida del Numerical Recipes, da una distribucion de
! velocidades de Maxwell-Boltzmann
!*******************************************************************************
  function gasdev()
  implicit none
    real(pr) :: gasdev
    integer  :: iset
    real(pr) :: fac, gset, rsq, v1, v2, ran1
    SAVE        iset, gset
    DATA        iset/0/
    if (iset/=1) iset = 0
    if (iset/=0) iset = 1
    if (iset.eq.0) then
1     v1 = 2._pr*grnd()-1._pr
      v2 = 2._pr*grnd()-1._pr
      rsq = v1**2 + v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 1
      fac = sqrt(-2._pr*log(rsq)/rsq)
      gset = v1*fac
      gasdev = v2*fac
      iset = 1
    else
      gasdev = gset
      iset = 0
    endif
    return
  end function gasdev

end program p1
