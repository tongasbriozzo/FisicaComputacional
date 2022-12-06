!*******************************************************************************
!
!  Programa Coeficiente de Difusión en función de Rho
!
!  Este programa esta diseñado para obtener el coeficiente de difusión D en
!  función de la densidad Rho del sistema
!
!  Gaston Briozzo, Física Computacional, FaMAFyC, UNC, 08/07/2021
!
!*******************************************************************************

program cdr
use precision, only: pr
use mtmod,     only: grnd
use dbmod
implicit none

  integer                           :: i, ns, nr
  integer, parameter                :: nl=8, N=nl*nl*nl, t_eq=1000, nsim=10000, nh=10000!, t_run=100000
  real(pr)                          :: rho, Vol, L, cdif, dcdif, cd
  real(pr), dimension(1:N,1:3)      :: x, v, xnew
  real(pr), dimension(1:N,1:3)      :: x1, x2
  real(pr), dimension(0:nh,1:N,1:3) :: xhis
  real(pr), parameter               :: pi=acos(-1._pr)
  real(pr), parameter               :: dt=0.002_pr
  real(pr), parameter               :: r_cut = 2.5_pr, r2cut=r_cut*r_cut
  real(pr), parameter               :: T=2.5_pr, eta=2.87_pr, D0=T/(3._pr*pi*eta)
  real(pr), parameter               :: f_fx = dt*D0/T, f_nx = sqrt(2._pr*D0*dt), f_fv = D0/T, f_nv = sqrt(2._pr*D0)

! Aquí anotaremos el valor del coeficiente de difusión D y su error estandar en función de rho
  open(201,file='cdr.d')

! Calculamos D para distintos valores de rho
  do nr = 1, 48

!   Inicializamos rho y los parametros que dependen de ella
    rho = real(nr,pr)*0.025_pr
    Vol=real(N,pr)/rho
    L=(Vol)**(1._pr/3._pr)

    write(*,*) 'nr', nr, 'rho', rho ! Supervizamos el avance del programa

!   Termalizamos el sistema
    call SC(N,L,nl,x)                                       ! Acomodamos las partículas según una red SC
    do i = 1, t_eq                                          ! Nos movemos hasta el paso t_eq
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)                                    ! No remplazamos automaticamente estos vectores por si llegaramos a necesitar los viejos
      call pbc(N,L,x)                                       ! Aplicamos las condiciones periodicas de contorno al vector posicion x
    end do

!   Inicializamos el historial de posiciones
    xhis(0,:,:) = x(:,:)
!   Inicializamos <D> y <D^2>
    cdif = 0._pr
    dcdif = 0._pr

!   Llenamos el vector xhis
    do i = 1, nh
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)
      xhis(i,:,:) = xnew(:,:)                               ! Recordamos las posiciones del sistema en el instante i
    end do

!   Promediamos sobre nsim simulaciones
    do ns = 1, nsim
!     Calculamos los desplazamientos cuadraticos medios entre el instante 0 e i
      x1(:,:) = xhis(0,:,:)                ! Seleccionamos el instante inicial
      do i = 1001, nh
        x2(:,:) = xhis(i,:,:)              ! Seleccionamos también el instante i
        cd      = dcm_fi(N,x1,x2)/(6._pr*real(i,pr)*dt)/D0
        cdif    = cdif  + cd
        dcdif   = dcdif + cd*cd
      end do
!     Hacemos avanzar cada vector en xhit un paso de verlet
      do i = 0, nh-1
        xhis(i,:,:) = xhis(i+1,:,:)                         ! Cada posición retrocede un instante
      end do
      call langevin(N,L,r2cut,f_fx,f_nx,f_fv,f_nv,x,xnew,v) ! La subrutina langevin nos da los nuevos vectores posicion y velocidad en funcion de los anteriores
      x(:,:) = xnew(:,:)
      xhis(nh,:,:) = xnew(:,:)                              ! Actualizamos el último instante
    end do

! Normalizamos y anotamos los resultados
    cdif  = cdif/real(nsim*(nh-1000),pr)
    dcdif = dcdif/real(nsim*(nh-1000),pr)
    dcdif = sqrt(dcdif - cdif*cdif)
    write(201,*) rho, cdif, dcdif
    write(*,*) rho, cdif, dcdif

  end do

  close(201)

end program cdr
