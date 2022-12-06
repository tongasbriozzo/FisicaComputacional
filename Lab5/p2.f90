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

program p2
use precision, only: pr
use mtmod,     only: grnd
implicit none

  integer                           :: i, j, np, ns
  integer, parameter                :: N=500, t_eq=1000, t_scal=50, t_run=100000, nl=5, nsim=1000, nr=500, nh=100000
  real(pr)                          :: T, Ti
  real(pr), dimension(1:t_eq)       :: S
  real(pr), dimension(1:nr)         :: g, fdr
  real(pr), dimension(1:N,1:3)      :: x, v, xnew, vnew
  real(pr), dimension(1:N,1:3)      :: x1, x2
  real(pr), dimension(1:nh)         :: dcm
  real(pr), dimension(0:nh,1:N,1:3) :: xhis
  real(pr), parameter               :: r_cut = 2.5_pr, r2cut=r_cut*r_cut, r_inv=1._pr/r_cut, r_i3=r_inv*r_inv*r_inv
  real(pr), parameter               :: rho = 0.8_pr, Vol=real(N,pr)/rho, L=(Vol)**(1._pr/3._pr), a=L/((N/4._pr)**(1._pr/3._pr))
  real(pr), parameter               :: dt=0.001_pr
  real(pr), parameter               :: pi=acos(-1._pr)

! Inicializamos el sistema
  T = 1.0_pr        ! Temperatura

!*******************************************************************************
! Problema 2b, parametro de orden cristalino
!*******************************************************************************

! Inicializamos el Parametro de Orden Cristalino
  S(:) = 0._pr

  open(522,file='poc.d')  ! Aqui anotaremos las medidas del sistema en funcion de los pasos evolucionados
  do ns = 1, nsim         ! Promediaremos sobre nsim simulaciones
    write(*,*) 'nsim poc', ns
    call FCC(x)           ! Acomodamos las particulas segùn una FCC
  ! Elegimos las componentes de la velocidad segun una distribucion gaussiana
    do j = 1, 3
      do i = 1, N
        v(i,j) = gasdev()
      end do
    end do
  ! Calculamos la temperatura cinetica del sistema y renormalizamos las velocidades para obtener la temperatura deseada
    Ti = Temperatura(v)
    v(:,:) = v(:,:)*sqrt(T/Ti)
    do i = 1, t_eq                      ! Nos movemos hasta el paso t_eq
      call verlet_v(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)                ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      v(:,:) = vnew(:,:)                ! Por eso hacemos el remplazo manual
      call pbc(x)                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
      call term_b(v)                    ! Termalizamos segun Berendsen
      S(i) = S(i) + poc(x)              ! Sumamos los parametros de orden cristalino
    end do
  end do
! Normalizamos y escribimos los resultados
  S(:) = S(:)/real(nsim,pr)
  do i = 1, t_eq
    write(522,*) i, S(i)                 ! Anotamos los valores obtenidos
  end do
  close(522)

!*******************************************************************************
! Problema 2a, funcion distribucion radial
!*******************************************************************************

  g(:) = 0._pr
  fdr(:) = 0._pr

  do i = 1, t_run
    write(*,*) 't_run', i
    call verlet_v(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)                ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
    v(:,:) = vnew(:,:)                ! Por eso hacemos el remplazo manual
    call pbc(x)                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
    call term_b(v)
    call disrad(x,g)
    fdr(:) = fdr(:) + g(:)
  end do
  fdr(:) = fdr(:)/real(t_run,pr)
  open(521,file='fdr.d')
  do i = 1, nr
    write(521,*) 0.5_pr*L/real(nr,pr)*real(i,pr), fdr(i)
  end do
  close(521)

!*******************************************************************************
! Problema 2c, Coeficiente de difusion
!*******************************************************************************

! Termalizamos el sistema
  call FCC(x)
  do j = 1, 3
    do i = 1, N
      v(i,j) = gasdev()
    end do
  end do
  Ti = Temperatura(v)
  v(:,:) = v(:,:)*sqrt(T/Ti)
  do i = 1, t_eq                      ! Nos movemos hasta el paso t_eq
    call verlet_v(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)                ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
    v(:,:) = vnew(:,:)                ! Por eso hacemos el remplazo manual
    call pbc(x)                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
    call term_b(v)                    ! Termalizamos segun Berendsen
  end do

! Inicializamos
  dcm(:) = 0._pr
  xhis(0,:,:) = x(:,:)

! Llenamos el vector xhis
  do i = 1, nh
    do j = 1, 1
      call verlet_vnp(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)                  ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      v(:,:) = vnew(:,:)                  ! Por eso hacemos el remplazo manual
    end do
    xhis(i,:,:) = xnew(:,:)             ! Recordamos las posiciones del sistema en el instante i
  end do

! Promediamos sobre nsim simulaciones
  do ns = 1, 1000
    write(*,*) 'nsim', ns
!   Calculamos los desplazamientos cuadraticos medios entre el instante 0 e i
    x1(:,:) = xhis(0,:,:)
    do i = 1, nh
      x2 = xhis(i,:,:)
      dcm(i) = dcm(i) + dist2(x1,x2)
    end do
!   Hacemos avanzar cada vector en xhit un paso de verlet
    do i = 0, nh-1
      xhis(i,:,:) = xhis(i+1,:,:)
    end do
    do j = 1, 1
      call verlet_vnp(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)                  ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      v(:,:) = vnew(:,:)                  ! Por eso hacemos el remplazo manual
    end do
    xhis(nh,:,:) = xnew(:,:)
  end do

! Normalizamos
  dcm(:) = dcm(:)/real(1000,pr)

! Anotamos los resultados
  open(523,file='dcm.d')
  do i = 1, nh
    write(523,*) real(i,pr)*dt, dcm(i)
  end do
  close(523)


!*******************************************************************************
!*******************************************************************************
CONTAINS
!*******************************************************************************
!*******************************************************************************

!*******************************************************************************
! La funcion dist2 calcula el desplazamiento cuadratico medio de un sistema de N
! particulas entre los instantes 1 y 2
!*******************************************************************************
  function dist2(xx1,xx2)
  implicit none
    real(pr), dimension(1:N,1:3) :: xx1, xx2
    real(pr)                     :: dist2
    integer                      :: ii
    real(pr)                     :: dx, dy, dz, r2
    r2 = 0._pr
    do ii = 1, N
      dx = xx2(ii,1) - xx1(ii,1)
      dy = xx2(ii,2) - xx1(ii,2)
      dz = xx2(ii,3) - xx1(ii,3)
      r2 = r2 + dx*dx + dy*dy + dz*dz
    end do
    dist2 = r2/real(N,pr)
  end function dist2

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
! La funcion fijnp es una modificacion de la funcion fij que no asume
! condiciones periodicas de contorno
!*******************************************************************************
  function fijnp(ii,jj,xx)
  implicit none
    integer, intent(in)                      :: ii, jj
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr), dimension    (1:3)             :: fijnp
    real(pr)                                 :: dx, dy, dz, r2
    if (jj==ii) then
      fijnp = 0._pr
      return
    end if
    dx = xx(ii,1) - xx(jj,1)
    dy = xx(ii,2) - xx(jj,2)
    dz = xx(ii,3) - xx(jj,3)
    r2 = dx*dx + dy*dy + dz*dz
    if (r2==0._pr) then
      fijnp = 0._pr
      return
    else if (r2>=r2cut) then
      fijnp = 0._pr
      return
    end if
    fijnp(1) = 48._pr*dx*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fijnp(2) = 48._pr*dy*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    fijnp(3) = 48._pr*dz*(1._pr/(r2*r2*r2*r2*r2*r2) - 0.5_pr/(r2*r2*r2))/r2
    return
  end function fijnp

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
! La funcion finp es una modificacion de la funcion fi que no asume condiciones
! de contorno periodicas
!*******************************************************************************
  function finp(ii,xx)
  implicit none
    integer, intent(in)                      :: ii
    integer                                  :: jj
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr), dimension    (1:3)             :: finp, sum_finp
    sum_finp(:) = 0._pr
    do jj = 1, N
      sum_finp = sum_finp + fijnp(ii,jj,xx)
    end do
    finp = sum_finp
  end function finp

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
! La funcion poc calcula el parametro de orden cristalino de un sistema de N
! particulas distribuidas segun el vetor xx para el vector k=(2pi/a)(-1,1,-1)
! en funcion del tiempo
!*******************************************************************************
  function poc(xx)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr)                                 :: poc
    integer                                  :: ii, jj
    real(pr)                                 :: sum_r, sum_i, kr
    sum_r = 0._pr
    sum_i = 0._pr
    do ii = 1, N
      kr = 2._pr*pi*(-xx(ii,1)+xx(ii,2)-xx(ii,3))/a
      sum_r = sum_r + cos(kr)
      sum_i = sum_i + sin(kr)
    end do
    poc = (sum_r*sum_r + sum_i*sum_i)/(N*N)
  end function poc

!*******************************************************************************
! la subrutina disrad da la funcion distribucion radial gg para un sistema de N
! particulas distribuidas segun el vector xx
!*******************************************************************************
  subroutine disrad(xx,gg)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in)  :: xx
    real(pr), dimension(1:nr),    intent(out) :: gg
    integer                                   :: ii, jj, kk
    real                                      :: ll, dx, dy, dz, dd, rmax, rmin
    gg(:) = 0._pr
    ll = 0.5_pr*L/real(nr,pr)
    do ii = 2, N
      do jj = 1, ii-1
        dx = xx(ii,1) - xx(jj,1)
        dy = xx(ii,2) - xx(jj,2)
        dz = xx(ii,3) - xx(jj,3)
        dx = dx - L*anint(dx/L)
        dy = dy - L*anint(dy/L)
        dz = dz - L*anint(dz/L)
        dd = sqrt(dx*dx + dy*dy + dz*dz)
        do kk = 1, nr
          if (real(kk-1,pr)*ll<dd.and.dd<=real(kk,pr)*ll) then
            gg(kk) = gg(kk) + 2._pr
          end if
        end do
      end do
    end do
    do kk = 1, nr
      rmin = real(kk-1,pr)*ll
      rmax = real(kk,pr)*ll
      gg(kk) = gg(kk)/N/(4._pr*pi*rho*(rmax*rmax*rmax-rmin*rmin*rmin)/3._pr)
    end do
    return
  end subroutine disrad

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
!    call pbc(xxnew)                        ! Incluimos las condiciones periodicas de contorno
    do ii = 1, N
      vvnew(ii,:) = vv(ii,:) + 0.5_pr*(fi(ii,xxnew)+fi(ii,xx))*dt
      vn(:) = vn(:) + vvnew(ii,:)
    end do
    nale = 1 + anint(real(N-1,pr)*grnd())  ! Elegimos una particula al azar
    vvnew(nale,:) = vvnew(nale,:) - vn(:)  ! Aseguramos velocidad del centro de masa nula
    return
  end subroutine verlet_v

!*******************************************************************************
! La subrutina verlet_vnp implementa en metodo Velocity-Verlet para obtener las
! posiciones y velocidades de un sistema dinamico de N particulas en el instante
! siguiente. Esta modificada para no asumir condiciones periodicas de borda
!*******************************************************************************
  subroutine verlet_vnp(xx,vv,xxnew,vvnew)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in)  :: xx, vv
    real(pr), dimension(1:N,1:3), intent(out) :: xxnew, vvnew
    integer                                   :: ii, jj
    do ii = 1, N
      xxnew(ii,:) = xx(ii,:) + vv(ii,:)*dt + 0.5_pr*finp(ii,xx)*dt*dt
    end do
    do ii = 1, N
      vvnew(ii,:) = vv(ii,:) + 0.5_pr*(fi(ii,xxnew)+finp(ii,xx))*dt
    end do
    return
  end subroutine verlet_vnp

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
! La subrutina termo_b aplica un baño termico modificando las velocidades de un
! sistema de N particulas, mediante el metodo de Berendsen
!*******************************************************************************
  subroutine term_b(vv)
  implicit none
    real(pr), dimension(1:N,1:3), intent(inout) :: vv
    real(pr)                                    :: lambda, tau, Tv
    Tv = Temperatura(vv)
    tau = 2._pr*dt
    lambda = sqrt( 1._pr + (dt/tau)*((T/Tv) - 1._pr) )
    vv(:,:) = lambda*vv(:,:)
    return
  end subroutine

!*******************************************************************************
! La subrutina FCC acomoda las posiciones de un sistema de N particulas para
! obtener una red cubica centrada en las caras (FCC)
! Profe mire nomas que beieza de red
!*******************************************************************************
  subroutine FCC(xx)
  implicit none
    real(pr), dimension(1:N,1:3), intent(inout) :: xx
    integer                                     :: ii, jj, kk
    real(pr)                                    :: lado
    lado = L/real(2*nl,pr)
    do ii = 1, 2*nl, 2
      do jj = 1, 2*nl, 2
        do kk = 1, 2*nl, 2
          xx(1 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 1) = lado*real(ii  ,pr)
          xx(1 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 2) = lado*real(jj  ,pr)
          xx(1 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 3) = lado*real(kk  ,pr)
          xx(2 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 1) = lado*real(ii  ,pr)
          xx(2 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 2) = lado*real(jj+1,pr)
          xx(2 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 3) = lado*real(kk+1,pr)
          xx(3 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 1) = lado*real(ii+1,pr)
          xx(3 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 2) = lado*real(jj  ,pr)
          xx(3 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 3) = lado*real(kk+1,pr)
          xx(4 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 1) = lado*real(ii+1,pr)
          xx(4 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 2) = lado*real(jj+1,pr)
          xx(4 + 4*(kk-1)/2 + 4*nl*(jj-1)/2 + 4*nl*nl*(ii-1)/2, 3) = lado*real(kk  ,pr)
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
    integer                                     :: ii
    real(pr)                                    :: ang1, ang2, modulo_v, Tv
    do ii = 1, N
      modulo_v = grnd()
      ang1 = grnd()*2._pr*pi
      ang2 = grnd()*2._pr*pi
      vv(ii,1) = modulo_v*sin(ang1)*cos(ang2)
      vv(ii,2) = modulo_v*sin(ang1)*sin(ang2)
      vv(ii,3) = modulo_v*cos(ang1)
    end do
    Tv = Temperatura(vv)
    vv(:,:) = vv(:,:)*sqrt(TT/Tv)
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

end program p2
