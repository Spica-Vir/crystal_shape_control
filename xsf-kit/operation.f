      module operation
        implicit none
        private
!       A2BR : Conversion rate from Angstrom to Bohr
!       NGDM : The maximum number of grid points along X/Y/Z
!        real,parameter    :: A2BR = 1.889726128 disabled
        integer,parameter :: NGDM = 1000

        public  :: planar_avg,shift_origin
        private :: plane_area,line_integ

        contains
        subroutine planar_avg(LATT,ATCOORD,ORG,BOX,GRID,AVGVEC,SHIFT,
     &                        AREA,DIST,AVG1D,INT1D)
!         Calculate the planar averaged data along the given direction
!         AVGVEC  : 1,2,3, lattice vectors along which the planar averaged data is computed
!         SHIFT   : float Length of shift of origin, in Å
!         AREA    : float Cross-sectional area of averaged plane, in Å^2
!         DIST    : NGDAVG * 1 Plane distances, in Å
!         AVG1D   : NGDAVG * 1 Planar averaged data, surface area normalized to 1
!         INT1D   : NGDAVG * 1 Integrated 1D profile, surface area normalized to 1
          real,dimension(:,:),intent(in)   :: ATCOORD
          real,dimension(3),intent(in)     :: ORG
          real,dimension(3,3),intent(in)   :: LATT,BOX
          real,dimension(:,:,:),intent(in) :: GRID
          integer,intent(in)               :: AVGVEC
          real,intent(in)                  :: SHIFT
          integer                          :: NGDAVG,NGDX,NGDY,NGDZ
          integer                          :: I,J,K
          real                             :: TDIST,DDIST
          real,intent(out)                 :: AREA
          real,dimension(:),allocatable,intent(out) :: DIST,AVG1D,INT1D

          NGDX = size(GRID,dim=1)
          NGDY = size(GRID,dim=2)
          NGDZ = size(GRID,dim=3)
          NGDAVG = size(GRID,dim=AVGVEC)
          AREA = plane_area(BOX,AVGVEC)

          TDIST = (BOX(1,AVGVEC)**2
     &           + BOX(2,AVGVEC)**2
     &           + BOX(3,AVGVEC)**2)**0.5
          DDIST = TDIST / NGDAVG

          allocate(DIST(NGDAVG),AVG1D(NGDAVG),INT1D(NGDAVG))
          do I = 1,NGDAVG
! Data point at the middle of a slice
            DIST(I) = DDIST * (I - 0.5) + ORG(NGDAVG)
! Initialize AVG1D
            AVG1D(I) = 0.
          enddo

          if (AVGVEC == 1) then
            do I = 1,NGDX
              do K = 1,NGDZ
                do J = 1,NGDY
                  AVG1D(I) = AVG1D(I) + GRID(I,J,K) / NGDY / NGDZ
                end do
              end do
            end do
          else if (AVGVEC == 2) then
            do J = 1,NGDY
              do K = 1,NGDZ
                do I = 1,NGDX
                  AVG1D(J) = AVG1D(J) + GRID(I,J,K) / NGDX / NGDZ
                end do
              end do
            end do
          else if (AVGVEC == 3) then
            do K = 1,NGDZ
              do J = 1,NGDY
                do I = 1,NGDX
                  AVG1D(K) = AVG1D(K) + GRID(I,J,K) / NGDX / NGDY
                end do
              end do
            end do
          end if

          call shift_origin (LATT,ATCOORD,AVGVEC,SHIFT,DIST,AVG1D)
          print*,'Planar averaged data calculated along ', AVGVEC

! Get integration
          INT1D = line_integ(AVG1D)
          print*,'Integration calculated along          ', AVGVEC
        end subroutine planar_avg

        subroutine shift_origin (LATT,ATCOORD,AVGVEC,SHIFT,DIST,AVG1D)
!         Shift origin of slab along the averaged direction
!         SHIFT : Shifting length. In Angstrom
          real,dimension(3,3),intent(in)  :: LATT
          real,dimension(:,:),intent(in)  :: ATCOORD
          integer,intent(in)              :: AVGVEC
          real,intent(in)                 :: SHIFT
          integer                         :: NATOM,NPT,I,J
          real                            :: FRAC,FRACMI=2.,FRACMX=-2.,
     &                                       FRACMID,LENVEC=0.,DTPDT,
     &                                       DISP,TMP
          real,dimension(3)               :: VEC
          real,dimension(:),intent(inout) :: DIST,AVG1D
!         Shift geometry
          NATOM = size(ATCOORD,dim=2)
          do I = 1,3
            VEC(I) = LATT(I,AVGVEC)
            LENVEC = LENVEC + VEC(I)**2
          end do
          LENVEC = LENVEC**0.5
          do I = 1,NATOM
            DTPDT = 0.
            do J = 1,3
              DTPDT = DTPDT + VEC(J) / LENVEC * ATCOORD(J,I)
            enddo
            FRAC = DTPDT / LENVEC
            if (FRAC < FRACMI) then
              FRACMI = FRAC
            endif
            if (FRAC > FRACMX) then
              FRACMX = FRAC
            endif
          end do
!         Shift data grid
          FRACMID = (FRACMX + FRACMI) / 2
          DISP = -LENVEC * FRACMID
          NPT = size(DIST)
          do I = 1,NPT
            FRAC = DIST(I) / LENVEC - FRACMID
            if (FRAC > 0.5) then
              DIST(I) = DIST(I) + DISP - LENVEC
            else if (FRAC < -0.5) then
              DIST(I) = DIST(I) + DISP + LENVEC
            else
              DIST(I) = DIST(I) + DISP
            endif
            DIST(I) = DIST(I) + SHIFT
          end do
!         Sort DIST,AVG1D
          do I = 1,NPT-1
            do J = I+1,NPT
              if (DIST(I) > DIST(J)) then
                TMP = DIST(I)
                DIST(I) = DIST(J)
                DIST(J) = TMP
                TMP = AVG1D(I)
                AVG1D(I) = AVG1D(J)
                AVG1D(J) = TMP
              endif
            end do
          end do
        end subroutine

        function plane_area(BOX,AVGVEC) result(AREA)
!         Calculate the area of the lattice plane defined by 2 base vectors other than AVGVEC
!         AREA : Area of the plane, in Å^2
          real,dimension(3,3),intent(in) :: BOX
          integer,intent(in)             :: AVGVEC
          real,dimension(3)              :: CROSP
          integer                        :: PVEC1,PVEC2
          real                           :: AREA

          if (AVGVEC == 1) then
            PVEC1 = 2
            PVEC2 = 3
          else if (AVGVEC == 2) then
            PVEC1 = 1
            PVEC2 = 3
          else if (AVGVEC == 3) then
            PVEC1 = 1
            PVEC2 = 2
          else
            print*,'Averaged direction must be 1, 2, or 3. Exiting.'
            stop
          endif

          CROSP(1)
     &      = BOX(2,PVEC1) * BOX(3,PVEC2) - BOX(3,PVEC1) * BOX(2,PVEC2)
          CROSP(2)
     &      = BOX(3,PVEC1) * BOX(1,PVEC2) - BOX(1,PVEC1) * BOX(3,PVEC2)
          CROSP(3)
     &      = BOX(1,PVEC1) * BOX(2,PVEC2) - BOX(2,PVEC1) * BOX(1,PVEC2)
          AREA
     &      = (CROSP(1)**2 + CROSP(2)**2 + CROSP(3)**2)**0.5

          return
        end function plane_area

        function line_integ(AVG1D) result(INT1D)
!         Integrate the line profile of 1D data
!         AVG1D : NGDAVG * 1 Planar averaged data, surface area normalized to 1
!         INT1D : NGDAVG * 1 Integrated 1D profile, surface area normalized to 1
          real,dimension(:),intent(in)              :: AVG1D
          real                                      :: SUMMED=0.
          integer                                   :: NGAVG,I
          real,dimension(:),allocatable             :: INT1D

          NGAVG = size(AVG1D)
          allocate(INT1D(NGAVG))
          do I = 1,NGAVG
            SUMMED = SUMMED + AVG1D(I)
            INT1D(I) = SUMMED
          end do
          return
        end function line_integ
      end module operation


