!*******************************************************************************
!
!  Programa Función Distribución Radial g(r), Parámetro de Órden Cristalino S(k)
!  y Coeficiente de Difusión D
!
!  Este programa esta diseñado para obtener las funciones g(r), S(k,t) y el
!  parámetro D, empleando un modelo de Dinámica Browniana
!
!  Gaston Briozzo, Física Computacional, FaMAFyC, UNC, 08/07/2021
!
!*******************************************************************************

program gsd
use precision, only: pr
use mtmod,     only: grnd
use dbmod
implicit none

  integer                           :: i, j, ns
  integer, parameter                :: nl=8, N=nl*nl*nl, t_eq=500, t_run=100000, nsim=1000, nr=500, nh=100000
  real(pr), dimension(1:t_eq)       :: S
  real(pr), dimension(1:nr)         :: g, fdr
  real(pr), dimension(1:N,1:3)      :: x, v, xnew
  real(pr), dimension(1:N,1:3)      :: x1, x2
  real(pr), dimension(1:nh)         :: dcm
  real(pr), dimension(0:nh,1:N,1:3) :: xhis
  real(pr), dimension(1:3)          :: krr
  real(pr), parameter               :: pi=acos(-1._pr)
  real(pr), parameter               :: dt=0.002_pr
  real(pr), parameter               :: r_cut = 2.5_pr, r2cut=r_cut*r_cut
  real(pr), parameter               :: T=2.5_pr, rho=1.2_pr, eta=2.87_pr
  real(pr), parameter               :: Vol=real(N,pr)/rho, L=(Vol)**(1._pr/3._pr), D0=T/(3._pr*pi*eta)
  real(pr), parameter               :: f_fx = dt*D0/T, f_nx = sqrt(2._pr*D0*dt), f_fv = D0/T, f_nv = sqrt(2._pr*D0)

!*******************************************************************************
! Parametro de Orden Cristalino
!*******************************************************************************

! Inicializamos el Parametro de Orden Cristalino
  S(:) = 0._pr
  krr(:) = 16._pr*pi/L

  open(522,file='poc.d')  ! Aqui anotaremos las medidas del sistema en funcion de los pasos evolucionados

! Promediamos sobre nsim simulaciones independientes
  do ns = 1, nsim
    write(*,*) 'nsim poc', ns                               ! Supervizamos el avance del programa
    call SC(N,L,nl,x)                                       ! Acomodamos las particulas segùn una red SC
    do i = 1, t_eq                                          ! Nos movemos hasta el paso t_eq
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)                                    ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      call pbc(N,L,x)                                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
      S(i) = S(i) + poc(N,krr,x)                            ! Sumamos los parametros de orden cristalino
    end do
  end do

! Normalizamos y escribimos los resultados
  S(:) = S(:)/real(nsim,pr)
  do i = 1, t_eq
    write(522,*) i*dt, S(i)
  end do
  close(522)

!*******************************************************************************
! Funcion Distribucion Radial
!*******************************************************************************

! Inicializamos las distribuciones
  g(:) = 0._pr
  fdr(:) = 0._pr

! Promediamos sobre t_run pasos
  do i = 1, t_run
    write(*,*) 't_run', i                                 ! Supervizamos el avance del programa
    call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)                                    ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
    call pbc(N,L,x)                                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
    call disrad(N,nr,L,pi,rho,x,g)                        ! La función disrad nos da la distribución radial promediada en un dado instante
    fdr(:) = fdr(:) + g(:)                                ! Promediamos la distribución radial sobre cada uno de los t_run instantes
  end do

! Normalizamos y anotamos los resultados
  fdr(:) = fdr(:)/real(t_run,pr)
  open(521,file='fdr.d')
  do i = 1, nr
    write(521,*) 0.5_pr*L/real(nr,pr)*real(i,pr), fdr(i)
  end do
  close(521)

!*******************************************************************************
! Coeficiente de Difusion
!*******************************************************************************

! Termalizamos el sistema
  call SC(N,L,nl,x)                                       ! Acomodamos las partículas según una red SC
  do i = 1, t_eq                                          ! Nos movemos hasta el paso t_eq
    call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)                                    ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
    call pbc(N,L,x)                                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
  end do

! Inicializamos
  dcm(:) = 0._pr       ! Desplazamiento cuadrático medio
  xhis(0,:,:) = x(:,:) ! Historial de posiciones

! Llenamos el vector xhis
  do i = 1, nh
    call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)
    xhis(i,:,:) = xnew(:,:)                               ! Recordamos las posiciones del sistema en el instante i
  end do

! Promediamos sobre nsim simulaciones
  do ns = 1, nsim
    write(*,*) 'nsim', ns ! Supervizamos el avance del programa
!   Calculamos los desplazamientos cuadraticos medios entre el instante 0 e i
    x1(:,:) = xhis(0,:,:)                ! Seleccionamos el instante inicial
    do i = 1, nh
      x2(:,:) = xhis(i,:,:)              ! Seleccionamos también el instante i
      dcm(i) = dcm(i) + dcm_fi(N,x1,x2)  ! La funcion dcm_fi da el desplazamiento cuadrático medio entre las posiciones x1 y x2
    end do
!   Hacemos avanzar cada vector en xhit un paso de verlet
    do i = 0, nh-1
      xhis(i,:,:) = xhis(i+1,:,:)                         ! Cada posición retrocede un instante
    end do
    call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
    x(:,:) = xnew(:,:)
    xhis(nh,:,:) = xnew(:,:)                              ! Actualizamos el último instante
  end do

! Normalizamos y anotamos los resultados
  dcm(:) = dcm(:)/real(nsim,pr)
  open(523,file='dcm.d')
  do i = 1, nh
    write(523,*) real(i,pr)*dt*D0, dcm(i)
  end do
  close(523)

end program gsd
