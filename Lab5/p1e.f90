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

program p1e
use precision, only: pr
use mtmod,     only: grnd
implicit none

  integer                         :: i, j, np, nc, nh, N, n1, n2, nl, Ncell, m
  integer                         :: teq,          trun,          tscal
  real(pr), parameter             :: t_eq=5.04_pr, t_run=5.04_pr, t_scal=0.252_pr
  real(pr)                        :: T, Ti,   Ec,  Ep,  Et,  P
  real(pr)                        ::   dTi,  dEc, dEp, dEt, dP
  real(pr)                        ::    Tim,  Ecm, Epm, Etm, Pm
  real(pr)                        ::    Ti2,  Ec2, Ep2, Et2, P2
  real(pr), allocatable           :: x(:,:), v(:,:), xnew(:,:), vnew(:,:)
  real(pr), allocatable           :: Head(:), list(:), map(:)
  real(pr), dimension(-50:49,1:3) :: his
  real(pr)                        :: r_cut, r2cut, r_inv, r_i3
  real(pr)                        :: rho, Vol, L
  real(pr)                        :: dt
  real(pr)                        :: U_tail, P_tail
  real(pr), parameter             :: pi=acos(-1._pr)
  real(pr)                        :: tiempo1, tiempo2, tinicializar, ttermalizar, tcorrer, ttotal

! Inicializamos el sistema
  T = 1.1_pr        ! Temperatura

  rho = 0.8_pr
  N   = 256
  nl  = 4
  Vol = real(N,pr)/rho
  L   = (Vol)**(1._pr/3._pr)

  r_cut = 2.5_pr
  r2cut = r_cut*r_cut
  r_inv = 1._pr/r_cut
  r_i3  = r_inv*r_inv*r_inv

  U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr

  allocate(x(1:N,1:3) , v(1:N,1:3) , xnew(1:N,1:3) , vnew(1:N,1:3) )

!*******************************************************************************
! Problema 1e
!*******************************************************************************

  open(521,file='dt.d')
  open(522,file='dt+.d')

  do n1 = 0, 60

    write(*,*) 'dt', n1

    dt    = 0.0001_pr*(1.122018454_pr**n1)
    teq   = int(t_eq/dt)
    trun  = int(t_run/dt)
    tscal = floor(real(teq,pr)/20._pr)

!   Inicializamos el sistema
    call FCC(x)         ! Acomodamos las particulas segùn una FCC
!   Elegimos las componentes de la velocidad segun una distribucion gaussiana
    do nc = 1, 3
      do np = 1, N
        v(np,nc) = gasdev()
      end do
    end do
!   Calculamos la temperatura cinetica del sistema y renormalizamos las velocidades para obtener la temperatura deseada
    Ti = Temperatura(v)
    v(:,:) = v(:,:)*sqrt(T/Ti)

    do i = 1, teq                       ! Nos movemos hasta el paso teq
    call verlet_v(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)                ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      v(:,:) = vnew(:,:)                ! Por eso hacemos el remplazo manual
      if (mod(i,tscal)==0) then         ! Cada tscal pasos, rescaleamos la velocidad para conservar la temperatura (Baño termico)
        Ti = Temperatura(v)
        v(:,:) = v(:,:)*sqrt(T/Ti)
      end if
    end do

!   Inicializamos las sumatorias
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

    do i = 1, trun                           ! Corremos sobre t_run pasos
      call verlet_v(x,v,xnew,vnew)
      x(:,:) = xnew(:,:)
      v(:,:) = vnew(:,:)
!     Calculamos valores instantaneos
      Ti = Temperatura(v)
      Ec = K(v)
      Ep = U(x) + N*U_tail
      Et = Ec + Ep
      P  = rho*Ti + pf(x)/(3._pr*Vol) + P_tail
!     Tomamos los valores medios
      Tim = Tim + Ti
      Ti2 = Ti2 + Ti*Ti
      Ecm = Ecm + Ec
      Ec2 = Ec2 + Ec*Ec
      Epm = Epm + Ep
      Ep2 = Ep2 + Ep*Ep
      Etm = Etm + Et
      Et2 = Et2 + Et*Et
      Pm  = Pm  + P
      P2  = P2  + P*P
    end do

!   Normalizamos los valores medios
    Tim = Tim/real(trun,pr)
    Ti2 = Ti2/real(trun,pr)
    Ecm = Ecm/real(trun,pr)
    Ec2 = Ec2/real(trun,pr)
    Epm = Epm/real(trun,pr)
    Ep2 = Ep2/real(trun,pr)
    Etm = Etm/real(trun,pr)
    Et2 = Et2/real(trun,pr)
    Pm  = Pm /real(trun,pr)
    P2  = P2 /real(trun,pr)
!   Calculamos las varianzas
    dTi = sqrt(Ti2-Tim*Tim)
    dEc = sqrt(Ec2-Ecm*Ecm)
    dEp = sqrt(Ep2-Epm*Epm)
    dEt = sqrt(Et2-Etm*Etm)
    dP  = sqrt( P2 -Pm* Pm)
!   Anotamos valores medios y varianzas concluyeno con el problema 1b

    write(521,*) dt, dEt/dEc, dEt/dEp
    write(522,*) dt, Tim, Ecm, Epm, Etm, Pm
    write(522,*) dt, dTi, dEc, dEp, dEt, dP

  end do

  close(521)
  close(522)

!*******************************************************************************
! Problema 1f
!*******************************************************************************

  dt    = 0.002_pr
  teq   = 2500
  trun  = 2500
  tscal = 125

  open(523,file='rc.d')
  open(524,file='rc+.d')

  do n1 = 0, 60

    write(*,*) 'r_cut', n1

    r_cut = 0.5_pr + 0.1_pr*real(n1,pr)
    r2cut = r_cut*r_cut
    r_inv = 1._pr/r_cut
    r_i3  = r_inv*r_inv*r_inv

    U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
    P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr

!    write(*,*) n1, r_cut, U_tail, P_tail

!   Inicializamos el sistema
    call FCC(x)         ! Acomodamos las particulas segùn una FCC
!   Elegimos las componentes de la velocidad segun una distribucion gaussiana
    do nc = 1, 3
      do np = 1, N
        v(np,nc) = gasdev()
      end do
    end do
!   Calculamos la temperatura cinetica del sistema y renormalizamos las velocidades para obtener la temperatura deseada
    Ti = Temperatura(v)
    v(:,:) = v(:,:)*sqrt(T/Ti)

    do i = 1, teq                       ! Nos movemos hasta el paso teq
      call verlet_v(x,v,xnew,vnew)      ! La subrutina verlet_v nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)                ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      v(:,:) = vnew(:,:)                ! Por eso hacemos el remplazo manual
      if (mod(i,tscal)==0) then         ! Cada tscal pasos, rescaleamos la velocidad para conservar la temperatura (Baño termico)
        Ti = Temperatura(v)
        v(:,:) = v(:,:)*sqrt(T/Ti)
      end if
    end do

!   Inicializamos las sumatorias
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

    do i = 1, trun                           ! Corremos sobre t_run pasos
      call verlet_v(x,v,xnew,vnew)
      x(:,:) = xnew(:,:)
      v(:,:) = vnew(:,:)
!     Calculamos valores instantaneos
      Ti = Temperatura(v)
      Ec = K(v)
      Ep = U(x) + N*U_tail
      Et = Ec + Ep
      P  = rho*Ti + pf(x)/(3._pr*Vol) + P_tail
!     Tomamos los valores medios
      Tim = Tim + Ti
      Ti2 = Ti2 + Ti*Ti
      Ecm = Ecm + Ec
      Ec2 = Ec2 + Ec*Ec
      Epm = Epm + Ep
      Ep2 = Ep2 + Ep*Ep
      Etm = Etm + Et
      Et2 = Et2 + Et*Et
      Pm  = Pm  + P
      P2  = P2  + P*P
    end do

!   Normalizamos los valores medios
    Tim = Tim/real(trun,pr)
    Ti2 = Ti2/real(trun,pr)
    Ecm = Ecm/real(trun,pr)
    Ec2 = Ec2/real(trun,pr)
    Epm = Epm/real(trun,pr)
    Ep2 = Ep2/real(trun,pr)
    Etm = Etm/real(trun,pr)
    Et2 = Et2/real(trun,pr)
    Pm  = Pm /real(trun,pr)
    P2  = P2 /real(trun,pr)
!   Calculamos las varianzas
    dTi = sqrt(Ti2-Tim*Tim)
    dEc = sqrt(Ec2-Ecm*Ecm)
    dEp = sqrt(Ep2-Epm*Epm)
    dEt = sqrt(Et2-Etm*Etm)
    dP  = sqrt( P2 -Pm* Pm)
!   Anotamos valores medios y varianzas concluyeno con el problema 1b

    write(523,*) r_cut, dEt/dEc, dEt/dEp
    write(524,*) r_cut, Tim, Ecm, Epm, Etm, Pm
    write(524,*) r_cut, dTi, dEc, dEp, dEt, dP

  end do

  close(523)
  close(524)

!  call cpu_time(tiempo2)

!  write(*,*) 'Tiempo', tiempo2-tiempo1

!*******************************************************************************
! Problema 1 g
!*******************************************************************************

  r_cut = 2.5_pr
  r2cut = r_cut*r_cut
  r_inv = 1._pr/r_cut
  r_i3  = r_inv*r_inv*r_inv

  U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr

  dt    = 0.005_pr
  teq   = 1000
  trun  = 1000
  tscal = 50

  deallocate(x,         v,         xnew,         vnew)

  open(525,file='tn.d')

  do i = 1, 20

    write(*,*) 'Numero de particulas', 4*i*i*i

    N   = 4*i*i*i
    nl  = i
    Vol = real(N,pr)/rho
    L   = (Vol)**(1._pr/3._pr)
    allocate  (x(1:N,1:3),v(1:N,1:3),xnew(1:N,1:3),vnew(1:N,1:3))

!  Vemos cuanto tiempo toma inicializar cada sistema
    call cpu_time(tiempo1)
      call FCC(x)
      do nc = 1, 3
        do np = 1, N
          v(np,nc) = gasdev()
        end do
      end do
    call cpu_time(tiempo2)
    tinicializar = tiempo2 - tiempo1

!   Vemos cuanto tarda en termalizar el sistema
    call cpu_time(tiempo1)
      call verlet_v(x,v,xnew,vnew)
      x(:,:) = xnew(:,:)
      v(:,:) = vnew(:,:)
      Ti = Temperatura(v)
      Ec = K(v)
      Ep = U(x) + N*U_tail
      Et = Ec + Ep
      P  = rho*Ti + pf(x)/(3._pr*Vol) + P_tail
    call cpu_time(tiempo2)
    ttermalizar = 1000._pr*(tiempo2-tiempo1)
    call cpu_time(tiempo1)
      Ti = Temperatura(v)
      v(:,:) = v(:,:)*sqrt(T/Ti)
    call cpu_time(tiempo2)
    ttermalizar = ttermalizar + 20._pr*(tiempo2-tiempo1)

!   Ahora veamos cuanto tarda en correr
    call cpu_time(tiempo1)
      call verlet_v(x,v,xnew,vnew)
      x(:,:) = xnew(:,:)
      v(:,:) = vnew(:,:)
      Ti = Temperatura(v)
      Ec = K(v)
      Ep = U(x) + N*U_tail
      Et = Ec + Ep
      P  = rho*Ti + pf(x)/(3._pr*Vol) + P_tail
      Tim = Tim + Ti
      Ti2 = Ti2 + Ti*Ti
      Ecm = Ecm + Ec
      Ec2 = Ec2 + Ec*Ec
      Epm = Epm + Ep
      Ep2 = Ep2 + Ep*Ep
      Etm = Etm + Et
      Et2 = Et2 + Et*Et
      Pm  = Pm  + P
      P2  = P2  + P*P
    call cpu_time(tiempo2)
    tcorrer = 2000._pr*(tiempo2-tiempo1)

    ttotal = tinicializar + ttermalizar + tcorrer

    write(525,*) i, N, tinicializar, ttermalizar, tcorrer, ttotal

    deallocate(x, v, xnew, vnew)

  end do

  close(525)

!*******************************************************************************
! Problema 1h, linked list
!*******************************************************************************

  r_cut = 2.5_pr
  r2cut = r_cut*r_cut
  r_inv = 1._pr/r_cut
  r_i3  = r_inv*r_inv*r_inv
  U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  dt    = 0.005_pr
  teq   = 1000
  trun  = 1000
  tscal = 50

  open(526,file='ll.d')

  do i = 5, 20

    write(*,*) 'Numero de particulas', 4*i*i*i

    N   = 4*i*i*i
    nl  = i
    Vol = real(N,pr)/rho
    L   = (Vol)**(1._pr/3._pr)
    m = floor(L/r_cut)
    ncell = m*m*m
    allocate  (x(1:N,1:3),v(1:N,1:3),xnew(1:N,1:3),vnew(1:N,1:3))
    allocate(head(1:ncell),list(1:N),map(1:13*ncell))
    call maps

!  Vemos cuanto tiempo toma inicializar cada sistema
    call cpu_time(tiempo1)
      call FCC(x)
      call pbc(x)
      do nc = 1, 3
        do np = 1, N
          v(np,nc) = gasdev()
        end do
      end do
    call cpu_time(tiempo2)
    tinicializar = tiempo2 - tiempo1

!   Vemos cuanto tarda en termalizar el sistema
    call cpu_time(tiempo1)
      call links
      call verlet_vll(x,v,xnew,vnew)
      x(:,:) = xnew(:,:)
      v(:,:) = vnew(:,:)
      call pbc(x)
      Ti = Temperatura(v)
      Ec = K(v)
      Ep = Ull(x) + N*U_tail
      Et = Ec + Ep
      P  = rho*Ti + pfll(x)/(3._pr*Vol) + P_tail
    call cpu_time(tiempo2)
    ttermalizar = 1000._pr*(tiempo2-tiempo1)
    call cpu_time(tiempo1)
      Ti = Temperatura(v)
      v(:,:) = v(:,:)*sqrt(T/Ti)
    call cpu_time(tiempo2)
    ttermalizar = ttermalizar + 20._pr*(tiempo2-tiempo1)

!   Ahora veamos cuanto tarda en correr
    call cpu_time(tiempo1)
      call links
      call verlet_vll(x,v,xnew,vnew)
      x(:,:) = xnew(:,:)
      v(:,:) = vnew(:,:)
      call pbc(x)
      Ti = Temperatura(v)
      Ec = K(v)
      Ep = Ull(x) + N*U_tail
      Et = Ec + Ep
      P  = rho*Ti + pfll(x)/(3._pr*Vol) + P_tail
      Tim = Tim + Ti
      Ti2 = Ti2 + Ti*Ti
      Ecm = Ecm + Ec
      Ec2 = Ec2 + Ec*Ec
      Epm = Epm + Ep
      Ep2 = Ep2 + Ep*Ep
      Etm = Etm + Et
      Et2 = Et2 + Et*Et
      Pm  = Pm  + P
      P2  = P2  + P*P
    call cpu_time(tiempo2)
    tcorrer = 2000._pr*(tiempo2-tiempo1)

    ttotal = tinicializar + ttermalizar + tcorrer

    write(526,*) i, N, tinicializar, ttermalizar, tcorrer, ttotal

    deallocate(x,         v,         xnew,         vnew)
    deallocate(head,list,map)

  end do

  close(526)

!*******************************************************************************
!*******************************************************************************
CONTAINS
!*******************************************************************************
!*******************************************************************************

!*******************************************************************************
! Primero haremos un bloque con las funciones para linked-list
!*******************************************************************************

!*******************************************************************************
! La funcion i_cell se emplea en el mapeo inicial
!*******************************************************************************
  function i_cell(ix,iy,iz)
  implicit none
    integer, intent(in) :: ix, iy, iz
    integer             :: i_cell
    i_cell = 1 + mod(ix-1+m,m) &
               + mod(iy-1+m,m)*m &
               + mod(iz-1+m,m)*m*m
  end function i_cell

!*******************************************************************************
! La subrutina maps mapea los cuadrantes vecinos
!*******************************************************************************
  subroutine maps
  implicit none
    integer :: ix, iy, iz, imap
    do ix = 1, M
      do iy = 1, M
        do iz = 1, M
          imap = (i_cell(ix,iy,iz)-1)*13
          map(imap + 1 ) = i_cell(ix+1,iy  ,iz  )
          map(imap + 2 ) = i_cell(ix+1,iy+1,iz  )
          map(imap + 3 ) = i_cell(ix  ,iy+1,iz  )
          map(imap + 4 ) = i_cell(ix-1,iy+1,iz  )
          map(imap + 5 ) = i_cell(ix+1,iy  ,iz-1)
          map(imap + 6 ) = i_cell(ix+1,iy+1,iz-1)
          map(imap + 7 ) = i_cell(ix  ,iy+1,iz-1)
          map(imap + 8 ) = i_cell(ix-1,iy+1,iz-1)
          map(imap + 9 ) = i_cell(ix+1,iy  ,iz+1)
          map(imap + 10) = i_cell(ix+1,iy+1,iz+1)
          map(imap + 11) = i_cell(ix  ,iy+1,iz+1)
          map(imap + 12) = i_cell(ix-1,iy+1,iz+1)
          map(imap + 13) = i_cell(ix  ,iy  ,iz+1)
        end do
      end do
    end do
  end subroutine maps

!*******************************************************************************
! La subrutina links genera la lista de particulas a visitar
!*******************************************************************************
  subroutine links
  implicit none
    real(pr) :: Lc_inv
    integer  :: ii, icell
    head(:) = 0
    Lc_inv = real(M,pr)/L
    do ii = 1, N
      icell = 1 + int((x(ii,1)-0.5_pr*L/M)*Lc_inv) &
                + int((x(ii,2)-0.5_pr*L/M)*Lc_inv)*M &
                + int((x(ii,3)-0.5_pr*L/M)*Lc_inv)*M*M
      list(ii) = head(icell)
      head(icell) = ii
    end do
  end subroutine links

!*******************************************************************************
! La funcion Uij da la energia potencial entre las patriculas i y j
!*******************************************************************************
  function Uij(ii,jj,xx)
  implicit none
    integer, intent(in)                      :: ii, jj
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr)                                 :: Uij
    real(pr)                                 :: dx, dy, dz, r2
    dx = xx(jj,1) - xx(ii,1)
    dy = xx(jj,2) - xx(ii,2)
    dz = xx(jj,3) - xx(ii,3)
    dx = dx - L*anint(dx/L)
    dy = dy - L*anint(dy/L)
    dz = dz - L*anint(dz/L)
    r2 = dx*dx + dy*dy + dz*dz
    if (r2==0._pr) then
!     write(*,*) 'Hubo distancia cero entre las particulas', ii, jj
      Uij = 0._pr
    else if (r2>=r2cut) then
      Uij = 0._pr
    end if
    Uij= 4._pr*(1._pr/(r2*r2*r2*r2*r2*r2) - 1._pr/(r2*r2*r2))
  end function Uij

!*******************************************************************************
! La funcion Ull calcula la energia potencial del sistema empleando el metodo
! linked list
!*******************************************************************************
  function Ull(xx)
  implicit none
    integer                                  :: ii, jj, icell, jcell, jcell0, nabor
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr)                                 :: Ull, sum_ull
    real(pr)                                 :: dx, dy, dz, r2
    sum_ull = 0._pr
    do icell = 1, Ncell
      ii = head(icell)
      do while (ii/=0)
        jj = list(ii)
        do while(jj/=0)
          sum_ull = sum_ull + Uij(ii,jj,xx)
          jj = list(jj)
        end do
        jcell0 = 13*(icell-1)
        do nabor = 1, 13
          jcell = map(jcell0+nabor)
          jj = head(jcell)
          do while (jj/=0)
            sum_ull = sum_ull + Uij(ii,jj,xx)
            jj = list(jj)
          end do
        end do
        ii = list(ii)
      end do
    end do
    Ull = sum_ull
  end function Ull

!*******************************************************************************
! La funcion fll calcula las fuerzas ejercidas sobre cada particula por las
! demas empleando el metodo linked list
!*******************************************************************************
  function fll(xx)
  implicit none
    integer                                  :: ii, jj, icell, jcell, jcell0, nabor
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr), dimension(1:N,1:3)             :: fll, sum_fll
    sum_fll(:,:) = 0._pr
    do icell = 1, ncell
      ii = head(icell)
      do while (ii/=0)
        jj = list(ii)
        do while(jj/=0)
          sum_fll(ii,:) = sum_fll(ii,:) + fij(ii,jj,xx)
          sum_fll(jj,:) = sum_fll(jj,:) + fij(ii,jj,xx)
          jj = list(jj)
        end do
        jcell0 = 13*(icell-1)
        do nabor = 1, 13
          jcell = map(jcell0+nabor)
          jj = head(jcell)
          do while (jj/=0)
            sum_fll(ii,:) = sum_fll(ii,:) + fij(ii,jj,xx)
            sum_fll(jj,:) = sum_fll(jj,:) + fij(ii,jj,xx)
            jj = list(jj)
          end do
        end do
        ii = list(ii)
      end do
    end do
    fll(:,:) = sum_fll(:,:)
  end function fll

!*******************************************************************************
! La funcion pfll da la suma de los productos escalares entre las fuerzas que se
! ejercen las particulas ii y jj y la distancia entre estas, empleando linked list
!*******************************************************************************
    function pfll(xx)
    implicit none
      integer                                  :: ii, jj, icell, jcell, jcell0, nabor
      real(pr), dimension(1:N,1:3), intent(in) :: xx
      real(pr)                                 :: pfll, sum_pfll
      sum_pfll = 0._pr
      do icell = 1, ncell
        ii = head(icell)
        do while (ii/=0)
          jj = list(ii)
          do while(jj/=0)
            sum_pfll = pfij(ii,jj,xx)
            jj = list(jj)
          end do
          jcell0 = 13*(icell-1)
          do nabor = 1, 13
            jcell = map(jcell0+nabor)
            jj = head(jcell)
            do while (jj/=0)
              sum_pfll = pfij(ii,jj,xx)
              jj = list(jj)
            end do
          end do
          ii = list(ii)
        end do
      end do
      pfll = sum_pfll
    end function pfll

!*******************************************************************************
! La subrutina verlet_vll implementa en metodo Velocity-Verlet para obtener las
! posiciones y velocidades de un sistema dinamico de N particulas en el instante
! siguiente. Esta diseñada de forma tal que tanto la posicion como la velocidad
! del centro de masa del sistema resulten nulas. Esta emplea linked list
!*******************************************************************************
  subroutine verlet_vll(xx,vv,xxnew,vvnew)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in)  :: xx, vv
    real(pr), dimension(1:N,1:3), intent(out) :: xxnew, vvnew
    real(pr), dimension(1:N,1:3)              :: fv, fn
    real(pr), dimension    (1:3)              :: xn, vn
    integer                                   :: ii, jj, nale
    fv = fll(xx)
    xn(:) = 0._pr
    vn(:) = 0._pr
    do ii = 1, N
      xxnew(ii,:) = xx(ii,:) + vv(ii,:)*dt + 0.5_pr*fv(ii,:)*dt*dt
      xn(:) = xn(:) + xxnew(ii,:)
    end do
    nale = 1 + anint(real(N-1,pr)*grnd())  ! Elegimos una particula al azar
    xxnew(nale,:) = xxnew(nale,:) - xn(:)  ! Aseguramos posicion del centro de masa nula
    fn = fll(xxnew)
    call pbc(xxnew)                        ! Incluimos las condiciones periodicas de contorno
    do ii = 1, N
      vvnew(ii,:) = vv(ii,:) + 0.5_pr*(fn(ii,:)+fv(ii,:))*dt
      vn(:) = vn(:) + vvnew(ii,:)
    end do
    nale = 1 + anint(real(N-1,pr)*grnd())  ! Elegimos una particula al azar
    vvnew(nale,:) = vvnew(nale,:) - vn(:)  ! Aseguramos velocidad del centro de masa nula
    return
  end subroutine verlet_vll


!*******************************************************************************
! Lo que sigue es el contains empleado para los problemas anteriores
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
! La funcion pfij da la suma de los productos escalares entre las fuerzas que se
! ejercen las particulas ii y jj y la distancia entre estas
!*******************************************************************************
  function pfij(ii,jj,xx)
  implicit none
    real(pr), dimension(1:N,1:3), intent(in) :: xx
    real(pr)                                 :: pfij
    integer, intent(in)                      :: ii, jj
    real(pr), dimension(1:3)                 :: ff
    real(pr)                                 :: dx, dy, dz
    dx = xx(ii,1) - xx(jj,1)
    dy = xx(ii,2) - xx(jj,2)
    dz = xx(ii,3) - xx(jj,3)
    dx = dx - L*anint(dx/L)
    dy = dy - L*anint(dy/L)
    dz = dz - L*anint(dz/L)
    ff = fij(ii,jj,xx)
    pfij = dx*ff(1) + dy*ff(2) + dz*ff(3)
  end function pfij

!*******************************************************************************
! La funcion pf da la suma de los productos escalares entre las fuerzas que se
! ejercen todas las particulas del sistema y la distancia entre estas
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

end program p1e
