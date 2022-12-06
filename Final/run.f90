!*******************************************************************************
!
!  Programa Run
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

  integer                           :: i, j, np, nh, ns
  integer, parameter                :: nl=8, N=nl*nl*nl, t_eq=500, t_run=500, nsim=1
  real(pr)                          :: iTi,  iEc, iEp, iEt, iP
  real(pr)                          :: dTi,  dEc, dEp, dEt, dP
  real(pr)                          ::  Tim,  Ecm, Epm, Etm, Pm
  real(pr)                          ::  Ti2,  Ec2, Ep2, Et2, P2
  real(pr), dimension(1:t_eq+t_run) ::  Ti,   Ec,  Ep,  Et,  P
  real(pr), dimension(1:N,1:3)      :: x, v, xnew, vnew
  real(pr), dimension(-50:49,1:3)   :: his
  real(pr), parameter               :: pi=acos(-1._pr)
  real(pr), parameter               :: dt=0.002_pr
  real(pr), parameter               :: r_cut = 2.5_pr, r2cut=r_cut*r_cut, r_inv=1._pr/r_cut, r_i3=r_inv*r_inv*r_inv
  real(pr), parameter               :: T=2.5_pr, rho=0.8_pr, eta=2.87_pr
  real(pr), parameter               :: Vol=real(N,pr)/rho, L=(Vol)**(1._pr/3._pr), D0=T/(3._pr*pi*eta)
  real(pr), parameter               :: f_fx = dt*D0/T, f_nx = sqrt(2._pr*D0*dt), f_fv = D0/T, f_nv = sqrt(2._pr*D0)
  real(pr), parameter               :: U_tail = 8._pr*pi*rho*(r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr
  real(pr), parameter               :: P_tail = 16._pr*pi*rho*rho*(2._pr*r_i3*r_i3*r_i3/3._pr-r_i3)/3._pr

!*******************************************************************************
! A continuacion se verifica el correcto funcionamiento del programa
!*******************************************************************************

!! Escribimos las variables paramétricas a modo de control
!  write(*,*) D0, vol, L, f_fx, f_nx, f_fv, f_nv
!  write(*,*) U_tail, P_tail
!
!! Comprobamos el funcionamiento de la red SC
!  call SC(N,L,nl,x)
!  open(301,file='SC.d')
!  do i = 1, N
!    write(301,*) x(i,1), x(i,2), x(i,3)
!  end do
!  close(301)

!*******************************************************************************
! Inicializamos las variables
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
! inicializamos el histograma
  his(:,:) = 0._pr

!*******************************************************************************
! Ahora, termalizaremos el sistema y lo estudiaremos en su estado de equilibrio
! Promediaremos el resultado sobre un total de nsim simulaciones independientes
!*******************************************************************************

  DO ns = 1, nsim              ! Comenzamos la simulación ns

    write(*,*) 'nsim', ns

    call SC(N,L,nl,x)          ! Acomodamos las partículas según una red SC
    call maxwell(N,T,v)

!   Termalizamos
    do i = 1, t_eq             ! Realizamos t_eq pasos
!      write(*,*) 't_eq', i
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,v,xnew,vnew) ! La subrutina langevin nos da los nuevos vectores posicion en funcion de los anteriores
!     Calculamos los valores instantaneos de las variables de estado
      iTi = Temperatura(N,v)                                ! Temperatura cinética
      iEc = K(N,v)                                          ! Energía cinética
      iEp = U(N,L,r2cut,x) + N*U_tail                       ! Energía potencial
      iEt = iEc + iEp                                       ! Energía total
      iP  = rho*iTi + pf(N,L,r2cut,x)/(3._pr*Vol) + P_tail  ! Presión
!     Promediamos sobre las nsim simulaciones
      Ti(i) = Ti(i) + iTi
      Ec(i) = Ec(i) + iEc
      Ep(i) = Ep(i) + iEp
      Et(i) = Et(i) + iEt
      P(i)  = P(i)  + iP
!     Actualizamos las posiciones
      x(:,:) = xnew(:,:)       ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      call pbc(N,L,x)          ! Aplicamos las condiciones periódicas de contorno
      v(:,:) = vnew(:,:)
    end do

!   Corremos
    do i = t_eq + 1, t_eq+t_run ! Realizamos sobre t_run pasos
!      write(*,*) 't_run', i
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,v,xnew,vnew)
!     Calculamos los valores instantaneos de las variables de estado
      iTi = Temperatura(N,v)
      iEc = K(N,v)
      iEp = U(N,L,r2cut,x) + N*U_tail
      iEt = iEc + iEp
      iP  = rho*iTi + pf(N,L,r2cut,x)/(3._pr*Vol) + P_tail
!     Promediamos sobre las nsim simulaciones
      Ti(i) = Ti(i) + iTi
      Ec(i) = Ec(i) + iEc
      Ep(i) = Ep(i) + iEp
      Et(i) = Et(i) + iEt
      P(i)  = P(i)  + iP
!     Tomamos los valores medios
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
!     Actualizamos las posiciones
      x(:,:) = xnew(:,:)       ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      call pbc(N,L,x)          ! Aplicamos las condiciones periódicas de contorno
      v(:,:) = vnew(:,:)
!     Ahora examinamos las velocidades para construir el histograma del problema 1c
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
  open(302,file='ins.d')
  do i = 1, t_eq+t_run
    write(302,*) i*dt, Ti(i), Ec(i), Ep(i), Et(i), P(i)
  end do
  close(302)

! Normalizamos el histograma y guardamos sus valores
  his(:,:) = his(:,:)/real(N,pr)/real(t_run,pr)/real(nsim,pr)
  open(303,file='his.d')
  do nh = -50, 49
    write(303,*) 0.1_pr*real(nh,pr),   0._pr,     0._pr,     0._pr
    write(303,*) 0.1_pr*real(nh,pr),   his(nh,1), his(nh,2), his(nh,3)
    write(303,*) 0.1_pr*real(nh+1,pr), his(nh,1), his(nh,2), his(nh,3)
    write(303,*) 0.1_pr*real(nh+1,pr), 0._pr,     0._pr,     0._pr
  end do
  close(303)

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
  open(304,file='med.d')
  write(304,*) Tim, Ecm, Epm, Etm, Pm
  write(304,*) dTi, dEc, dEp, dEt, dP
  close(304)

end program run
