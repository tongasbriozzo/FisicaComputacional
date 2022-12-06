!*******************************************************************************
!
!  Programa Errores
!
!  Este programa realiza un barrido de los valores de dt y r_cut, obteniendo las
!  variandas y errores asociados a las variables de estado del sistema
!
!  Gaston Briozzo, Física Computacional, FaMAFyC, UNC, 07/07/2021
!
!*******************************************************************************

program err
use precision, only: pr
use mtmod,     only: grnd
use dbmod
implicit none

  integer                         :: i, n1
  integer                         :: teq,          trun
  integer, parameter              :: nl=8, N=nl*nl*nl
  real(pr), parameter             :: t_eq=5.00_pr, t_run=5.00_pr
  real(pr)                        :: T, Ti,   Ec,  Ep,  Et,  P
  real(pr)                        ::   dTi,  dEc, dEp, dEt, dP
  real(pr)                        ::    Tim,  Ecm, Epm, Etm, Pm
  real(pr)                        ::    Ti2,  Ec2, Ep2, Et2, P2
  real(pr), dimension(1:N,1:3)    :: x, v, xnew
  real(pr)                        :: r_cut, r2cut, r_inv, r_i3
  real(pr)                        :: rho, Vol, L
  real(pr)                        :: dt
  real(pr)                        :: U_tail, P_tail, f_fx, f_nx, f_fv, f_nv
  real(pr), parameter             :: pi=acos(-1._pr)
  real(pr)                        :: D0, eta

! Inicializamos el sistema
  T   = 2.5_pr                 ! Temperatura
  rho = 0.8_pr                 ! Densidad
  eta = 2.87_pr                ! Viscocidad
  D0  = T/(3._pr*pi*eta)       ! Coef de difución a densidad nula
  Vol = real(N,pr)/rho         ! Volumen del sistema
  L   = (Vol)**(1._pr/3._pr)   ! Longitud del sistema
  r_cut = 2.5_pr               ! Radio de corte
  r2cut = r_cut*r_cut
  r_inv = 1._pr/r_cut
  r_i3  = r_inv*r_inv*r_inv

  U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr

!*******************************************************************************
! Dependencia de los errores con dt
!*******************************************************************************

! Abrimos los documentos de texto
  open(201,file='dt.d')  ! Aqui anotaremos los cocientes entre los errores de las energías
  open(202,file='dt+.d') ! Aquí anotaremos valores medios y desviaciones estándar de las variables de estado

! Barremos distintos valores de dt con n1
  do n1 = 0, 60

    write(*,*) 'dt', n1                     ! Supervizamos el avance del programa

    dt    = 0.0001_pr*(1.122018454_pr**n1)  ! Paso temporal
    teq   = int(t_eq/dt)                    ! Tiempo de termalización
    trun  = int(t_run/dt)                   ! Tiempo de corrida

!   Estos factores se repiten el las ecuaciones de Langevin, por lo que los calculamos una sola vez
    f_fx = dt*D0/T
    f_nx = sqrt(2._pr*D0*dt)
    f_fv = D0/T
    f_nv = sqrt(2._pr*D0)

!   Inicializamos el sistema
    call SC(N,L,nl,x)         ! Acomodamos las particulas segùn una red SC

!   Termalizamos
    do i = 1, teq                                           ! Nos movemos hasta el paso teq
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion en funcion de los anteriores
      x(:,:) = xnew(:,:)                                    ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      call pbc(N,L,x)                                       ! Aplicamos las condiciones periódicas de contorno
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

!   Corremos
    do i = 1, trun           ! Corremos sobre t_run pasos
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v)
!     Calculamos valores instantaneos
      Ti = Temperatura(N,v)                              ! Temperatura cinética
      Ec = K(N,v)                                        ! Energía cinética
      Ep = U(N,L,r2cut,x) + N*U_tail                     ! Energía potencial
      Et = Ec + Ep                                       ! Energía total
      P  = rho*Ti + pf(N,L,r2cut,x)/(3._pr*Vol) + P_tail ! Presión
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
!     Actualizamos las posiciones
      x(:,:) = xnew(:,:)
      call pbc(N,L,x)
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
    write(201,*) dt, dEt/dEc, dEt/dEp
    write(202,*) dt, Tim, dTi, Ecm, dEc, Epm, dEp, Etm, dEt, Pm, dP

  end do

  close(201)
  close(202)

!*******************************************************************************
! Dependencia de los errores con r_cut
!*******************************************************************************

! Inicializamos los parámetros
  dt    = 0.001_pr
  teq   = 5000
  trun  = 5000

  f_fx = dt*D0/T
  f_nx = sqrt(2._pr*D0*dt)
  f_fv = D0/T
  f_nv = sqrt(2._pr*D0)

! Abrimos los documentos de texto
  open(203,file='rc.d')  ! Aquí anotaremos los cocientes entre los errores de las energías
  open(204,file='rc+.d') ! Aquí anotaremos los valores medios y las desviaciones estándar de las variables de estado


! Barremos distintos valores de r_cut con n1
  do n1 = 0, 60

    write(*,*) 'r_cut', n1               ! Supervizamos el avance del programa

    r_cut = 0.5_pr + 0.1_pr*real(n1,pr)  ! Radio de corte
    r2cut = r_cut*r_cut
    r_inv = 1._pr/r_cut
    r_i3  = r_inv*r_inv*r_inv

    U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
    P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr


!   Inicializamos el sistema
    call SC(N,L,nl,x)         ! Acomodamos las particulas segùn una red SC

!   Termalizamos
    do i = 1, teq                                           ! Nos movemos hasta el paso teq
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion en funcion de los anteriores
      x(:,:) = xnew(:,:)                                    ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      call pbc(N,L,x)                                       ! Aplicamos las condiciones periódicas de contorno
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

!   Corremos
    do i = 1, trun           ! Corremos sobre t_run pasos
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v)
!     Calculamos valores instantaneos
      Ti = Temperatura(N,v)
      Ec = K(N,v)
      Ep = U(N,L,r2cut,x) + N*U_tail
      Et = Ec + Ep
      P  = rho*Ti + pf(N,L,r2cut,x)/(3._pr*Vol) + P_tail
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
!     Actualizamos las posiciones
      x(:,:) = xnew(:,:)
      call pbc(N,L,x)
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
    write(203,*) r_cut, dEt/dEc, dEt/dEp
    write(204,*) r_cut, Tim, dTi, Ecm, dEc, Epm, dEp, Etm, dEt, Pm, dP

end do

close(203)
close(204)

end program err
