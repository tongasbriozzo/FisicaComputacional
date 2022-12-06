! ______________________________________________________________________________
! CAMINATA ALEATORIA 2D
!
! Este programa simula una caminata aleatoria bidimensional, con probabilidad
! 1/4 de dar un paso unitario hacia adelante, atras, izda o drcha. Para obtener
! números al azar se emplean las subrutinas ran0, ran2, Marsaglia y Zaman y
! Mersenne Twister. Se simulan d2 caminatas de d1 pasos, cada una partiendo del
! origen, y se calculan trayectorias, desplazamiento cuadrático medio y estadías
! por cuadrante.
!
! Gaston Briozzo, FaMAFyC, UNC, 03/05/2021
! ______________________________________________________________________________

program p2
use precision, only: pr
use randommod, only: ran0, ran2
use mzranmod,  only: rmzran
use mtmod,     only: grnd
implicit none

  integer, parameter                :: d1 = 10000, d2 = 1000000
  integer                           :: i, j, k, sem1, sem2
  integer                           :: xp, yp
  integer, dimension(1:4)           :: x, y
  real(pr), dimension(1:4,1:d1)     :: dcm, dcmf
  real(pr), dimension(1:4,1:4,0:d1) :: cua, cuaf, ec, ecf
  real(pr), dimension(1:4)          :: u, cc

! Definimos los formatos que utilizaremos en los documentos de texto
  1 format (9(A4, 5x))
  2 format (9(I7, 2x))
  3 format (A4, 5x, 4(A6, 15x))
  4 format (I7, 2x, 4(E23.16, 2x))
  5 format (A4, 5x, 4(A7, 18x))
  6 format (I7, 2x, 4(E23.16, 2x))

! Inicializamos los parámetros
! Semillas
  sem1 = 1073741827  ! Es primo
  sem2 = 2073741823  ! Es primo
! Desplazamientos cuadráticos medios
  dcmf(:,:) = 0._pr
! Estadias por cuadrante
  cuaf(:,:,:) = 0._pr
! Error por cuadrante
  ecf(:,:,:) = 0._pr

! Abrimos los documentos de texto y describimos los datos a guardar
! Aqui van los datos de una unica caminata
  open(321,file='cam.d')      ! La caminata
  open(322,file='dcm1.d')     ! DCM
  open(323,file='cua1_r0.d')  ! Cuadrantes en r0
  open(324,file='cua1_r2.d')  ! Cuadrantes en r2
  open(325,file='cua1_mz.d')  ! Cuadrantes en mz
  open(326,file='cua1_mt.d')  ! Cuadrantes en mt
  open(327,file='ec1_r0.d')   ! Error por cuadrante en r0
  open(328,file='ec1_r2.d')   ! Error por cuadrante en r2
  open(329,file='ec1_mz.d')   ! Error por cuadrante en mz
  open(330,file='ec1_mt.d')   ! Error por cuadrante en mt
! Aqui van los datos promediados de d2 caminatas
  open(331,file='dcm.d')      ! DCM
  open(332,file='cua_r0.d')   ! Cuadrantes en r0
  open(333,file='cua_r2.d')   ! Cuadrantes en r2
  open(334,file='cua_mz.d')   ! Cuadrantes en mz
  open(335,file='cua_mt.d')   ! Cuadrantes en mt
  open(336,file='ec_r0.d')    ! Error por cuadrante en r0
  open(337,file='ec_r2.d')    ! Error por cuadrante en r2
  open(338,file='ec_mz.d')    ! Error por cuadrante en mz
  open(339,file='ec_mt.d')    ! Error por cuadrante en mt
! Ahora escribimos la primera línea para describir los datos
! Primero los de una sola caminata
  write(321,1) 'Paso', 'x_r0', 'y_r0', 'x_r2', 'y_r2', 'x_mz', 'y_mz', 'x_mt', 'y_mt'
  write(322,3) 'Paso', 'dcm_r0', 'dcm_r2', 'dcm_mz', 'dcm_mt'
  write(323,5) 'Paso', 'cua1_r0', 'cua2_r0', 'cua3_r0', 'cua4_r0'
  write(324,5) 'Paso', 'cua1_r2', 'cua2_r2', 'cua3_r2', 'cua4_r2'
  write(325,5) 'Paso', 'cua1_mz', 'cua2_mz', 'cua3_mz', 'cua4_mz'
  write(326,5) 'Paso', 'cua1_mt', 'cua2_mt', 'cua3_mt', 'cua4_mt'
  write(327,5) 'Paso', 'ec1_r0',  'ec2_r0',  'ec3_r0',  'ec4_r0'
  write(328,5) 'Paso', 'ec1_r2',  'ec2_r2',  'ec3_r2',  'ec4_r2'
  write(329,5) 'Paso', 'ec1_mz',  'ec2_mz',  'ec3_mz',  'ec4_mz'
  write(330,5) 'Paso', 'ec1_mt',  'ec2_mt',  'ec3_mt',  'ec4_mt'
! Despues los de d2 caminatas promediadas
  write(331,3) 'Paso', 'dcm_r0',  'dcm_r2',  'dcm_mz',  'dcm_mt'
  write(332,5) 'Paso', 'cua1_r0', 'cua2_r0', 'cua3_r0', 'cua4_r0'
  write(333,5) 'Paso', 'cua1_r2', 'cua2_r2', 'cua3_r2', 'cua4_r2'
  write(334,5) 'Paso', 'cua1_mz', 'cua2_mz', 'cua3_mz', 'cua4_mz'
  write(335,5) 'Paso', 'cua1_mt', 'cua2_mt', 'cua3_mt', 'cua4_mt'
  write(336,5) 'Paso', 'ec1_r0',  'ec2_r0',  'ec3_r0',  'ec4_r0'
  write(337,5) 'Paso', 'ec1_r2',  'ec2_r2',  'ec3_r2',  'ec4_r2'
  write(338,5) 'Paso', 'ec1_mz',  'ec2_mz',  'ec3_mz',  'ec4_mz'
  write(339,5) 'Paso', 'ec1_mt',  'ec2_mt',  'ec3_mt',  'ec4_mt'

! Repetimos la caminata d2 veces
  do i = 1, d2

!   Posición de partida
    x(:) = 0
    y(:) = 0
!   Desplazamiento, estadias y error por cuadrante en esta caminata
    dcm(:,:)    = 0
    cua(:,:,:)  = 0._pr
    ec(:,:,:)   = 0._pr

!   Calculamos el paso j de cada caminata
    do j = 1, d1

!     obtenemos los números u al azar (pseudo)
      u(1) = ran0(sem1)
      u(2) = ran2(sem2)
      u(3) = rmzran()
      u(4) = grnd()

!    De acuerdo al número u obtenido, decidimos hacia donde irá nuestro caminante
      do k = 1, 4
        call paso(u(k),xp,yp)
        x(k) = x(k) + xp
        y(k) = y(k) + yp
      end do

!     Calculamos los desplazamientos cuadráticos medios (dcm) y promediamos en el paso j, sobre las i caminatas
      dcm(:,j)  = real(x(:)*x(:) + y(:)*y(:),pr)  ! DCM en esta caminata
      dcmf(:,j) = dcmf(:,j) + dcm(:,j)            ! DCM promediado

!     Ahora contabilizamos cuantas veces estuvimos en cada cuadrante
      cua(:,:,j) = cua(:,:,j-1)  ! Acarremamos las estadias del instante anterior y sumamos la nueva
      do k = 1, 4
        call cuadrante(x(k),y(k),cc)
        cua(k,:,j) = cua(k,:,j) + cc(:)
      end do

!     Sumamos las estadias por cuadrante (normalizadas) para la estadia final
      cuaf(:,:,j) = cuaf(:,:,j) + cua(:,:,j)

!     Calculamos el error por cuadrante
      ec(:,:,j)  = abs(cua(:,:,j)-0.25_pr*real(j,pr))*4._pr/real(j,pr)

!     Si estamos en la primer caminata, la guardamos para tener una muestra no promediada
      if (i==1) then
        write(321,2) j, x(1), y(1), x(2), y(2), x(3), y(3), x(4), y(4)
        write(322,4) j, dcm(1,j),   dcm(2,j),   dcm(3,j),   dcm(4,j)
        write(323,6) j, cua(1,1,j), cua(1,2,j), cua(1,3,j), cua(1,4,j)
        write(324,6) j, cua(2,1,j), cua(2,2,j), cua(2,3,j), cua(2,4,j)
        write(325,6) j, cua(3,1,j), cua(3,2,j), cua(3,3,j), cua(3,4,j)
        write(326,6) j, cua(4,1,j), cua(4,2,j), cua(4,3,j), cua(4,4,j)
        write(327,6) j, ec(1,1,j),  ec(1,2,j),  ec(1,3,j),  ec(1,4,j)
        write(328,6) j, ec(2,1,j),  ec(2,2,j),  ec(2,3,j),  ec(2,4,j)
        write(329,6) j, ec(3,1,j),  ec(3,2,j),  ec(3,3,j),  ec(3,4,j)
        write(330,6) j, ec(4,1,j),  ec(4,2,j),  ec(4,3,j),  ec(4,4,j)
      end if

    end do

!   Supervizamos el avance del programa por terminal
    write(*,*) 'i=', i

  end do

! Normalizamos los desplazamientos y estadias
  dcmf(:,:)   = dcmf(:,:)/real(d2,pr)
  cuaf(:,:,:) = cuaf(:,:,:)/real(d1*d2,pr)
! Calculamos el error por cuadrante final
  do j = 1, d1
    ecf(:,:,j)   = abs(cuaf(:,:,j)*real(d1,pr)/real(j,pr)-0.25_pr)*4._pr
  end do

! Anotamos los datos obtenidos
  do j = 1, d1
    write(331,4) j,   dcmf(1,j),   dcmf(2,j),   dcmf(3,j),   dcmf(4,j)
    write(332,6) j, cuaf(1,1,j), cuaf(1,2,j), cuaf(1,3,j), cuaf(1,4,j)
    write(333,6) j, cuaf(2,1,j), cuaf(2,2,j), cuaf(2,3,j), cuaf(2,4,j)
    write(334,6) j, cuaf(3,1,j), cuaf(3,2,j), cuaf(3,3,j), cuaf(3,4,j)
    write(335,6) j, cuaf(4,1,j), cuaf(4,2,j), cuaf(4,3,j), cuaf(4,4,j)
    write(336,6) j, ecf(1,1,j),  ecf(1,2,j),  ecf(1,3,j),  ecf(1,4,j)
    write(337,6) j, ecf(2,1,j),  ecf(2,2,j),  ecf(2,3,j),  ecf(2,4,j)
    write(338,6) j, ecf(3,1,j),  ecf(3,2,j),  ecf(3,3,j),  ecf(3,4,j)
    write(339,6) j, ecf(4,1,j),  ecf(4,2,j),  ecf(4,3,j),  ecf(4,4,j)
  end do

CONTAINS

! ______________________________________________________________________________
! Subrutina paso
! Esta subrrutina seleccina, segun la entrada uu, hacia donde se dara el paso
! ______________________________________________________________________________

  subroutine paso(uu,xx,yy)
  implicit none
    real(pr), intent(in) :: uu
    integer              :: xx, yy
    xx = 0
    yy = 0
    if (uu<=0.25_pr) then
      xx = -1
    else if (uu<=0.5_pr) then
      xx = +1
    else if (uu<=0.75_pr) then
      yy = -1
    else
      yy = +1
    end if
  end subroutine paso

! ______________________________________________________________________________
! Subrutina cuadrante
! Esta subrutina nos dice en que cuadrante se encuentra el caminante
! ______________________________________________________________________________

  subroutine cuadrante(xx,yy,kk)
  implicit none
    integer,                  intent(in)  :: xx, yy
    real(pr), dimension(1:4), intent(out) :: kk
    kk(:) = 0._pr
    if      (yy>0  .and. xx>0) then
      kk(1) = 1._pr
    else if (yy>0  .and. xx<0) then
      kk(2) = 1._pr
    else if (yy<0  .and. xx<0) then
      kk(3) = 1._pr
    else if (yy<0  .and. xx>0) then
      kk(4) = 1._pr
    else if (xx==0 .and. yy>0) then
      kk(1) = 0.5_pr
      kk(2) = 0.5_pr
    else if (xx==0 .and. yy<0) then
      kk(3) = 0.5_pr
      kk(4) = 0.5_pr
    else if (yy==0 .and. xx>0) then
      kk(1) = 0.5_pr
      kk(4) = 0.5_pr
    else if (yy==0 .and. xx<0) then
      kk(2) = 0.5_pr
      kk(3) = 0.5_pr
    else if (xx==0 .and. yy==0) then
      kk(1) = 0.25_pr
      kk(2) = 0.25_pr
      kk(3) = 0.25_pr
      kk(4) = 0.25_pr
    end if
  end subroutine cuadrante

end program p2
