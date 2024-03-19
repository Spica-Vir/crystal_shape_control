      module option
        implicit none
        private
        integer,parameter :: NGDM = 1000
        integer,parameter :: NATM = 1000

        public :: option1,option2,option3,option4

        contains
        subroutine option1(INPUT,OUTPUT)
!         Get planar averaged data for a single XDF file
          use io
          use operation

          character(len=80),intent(in)         :: INPUT,OUTPUT
          integer                              :: AVGVEC
          character(len=10)                    :: SHIFTC
          real                                 :: SHIFT=0.0,AREA
          logical                              :: DOSHIFT
          real,dimension(3,3)                  :: LATT,BOX
          character*2,dimension(:),allocatable :: ATLABEL
          real,dimension(:,:),allocatable      :: ATCOORD
          real,dimension(3)                    :: ORG
          real,dimension(:,:,:),allocatable    :: GRID
          real,dimension(:),allocatable        :: DIST,AVG1D,INT1D

          print*,'Please specify the direction of 1D line profile: '
          print*,'  (1-3 only. 1=Lattice vector 1, 2=Lattice vector 2,
     &3=Lattice vector 3)'
          read*,AVGVEC

          if (AVGVEC /= 1 .and. AVGVEC /= 2 .and. AVGVEC /= 3) then
            print*,'Only the line profiles along a,b, or c are allowed.'
            stop
          endif

          print*,'Please specify the shift along the averaged direction
     &(Unit: Å): '
          print*,"  Use 'no' for no shift."
          read*,SHIFTC
          if (SHIFTC == 'no') then
            DOSHIFT = .false.
            SHIFT= 0.0
          else
            DOSHIFT = .true.
            read(SHIFTC,'(f10.6)') SHIFT
          endif

          call read_3dxsf(INPUT,LATT,ATLABEL,ATCOORD,ORG,BOX,GRID)
          call planar_avg(ORG,BOX,GRID,AVGVEC,AREA,DIST,AVG1D,INT1D)
          if (DOSHIFT) then
            call
     &        shift_origin(LATT,ATCOORD,AVGVEC,SHIFT,DIST,AVG1D,INT1D)
          endif
          call write_1dtxt(OUTPUT,AREA,SHIFT,DIST,AVG1D,INT1D)
        end subroutine option1
!----
        subroutine option2(INPUT0,OUTPUT)
!         Get difference of 3D data
          use io

          character(len=80),intent(in) :: INPUT0,OUTPUT
          character(len=80)   :: INPUT1
          real,dimension(3,3) :: LATT0,LATT1,BOX0,BOX1
          integer             :: NINPUT,I
          integer             :: NGDX0,NGDX1,NGDY0,NGDY1,NGDZ0,NGDZ1
          character*2,dimension(:),allocatable      :: ATLABEL0,ATLABEL1
          real,dimension(:,:),allocatable           :: ATCOORD0,ATCOORD1
          real,dimension(3)                         :: ORG0,ORG1
          real,dimension(:,:,:),allocatable         :: GRID0,GRID1

          print*,'Please specify the number of other 3D XSF files: '
          read*,NINPUT
          call
     &      read_3dxsf(INPUT0,LATT0,ATLABEL0,ATCOORD0,ORG0,BOX0,GRID0)
          NGDX0 = size(GRID0,dim=1)
          NGDY0 = size(GRID0,dim=2)
          NGDZ0 = size(GRID0,dim=3)

          do I = 1,NINPUT
            print*,'File name: '
            read*,INPUT1
            call
     &        read_3dxsf(INPUT1,LATT1,ATLABEL1,ATCOORD1,ORG1,BOX1,GRID1)
            NGDX1 = size(GRID1,dim=1)
            NGDY1 = size(GRID1,dim=2)
            NGDZ1 = size(GRID1,dim=3)
            if (NGDX0 /= NGDX1 .or. NGDY0 /= NGDY1 .or. NGDZ0 /= NGDZ1)
     &      then
              print*,'Inconsistent dimension detected between ',INPUT0,
     &         ' and ',INPUT1,'. Exiting.'
              stop
            endif
            GRID0 = GRID0 - GRID1
          enddo
          call
     &      write_3dxsf(OUTPUT,LATT0,ATLABEL0,ATCOORD0,ORG0,BOX0,GRID0)
        end subroutine option2
!----
        subroutine option3(INPUT0,OUT3,OUT1)
!----     Get difference of multiple 3D grid data and 1D line profile
          character(len=80),intent(in) :: INPUT0,OUT3,OUT1

          call option2(INPUT0,OUT3)
          call option1(OUT3,OUT1)
        end subroutine option3
!----
        subroutine option4(INPUT,OUTPUT)
!----     Normalize the integrated 3D grid data
          use io

          character(len=80),intent(in)         :: INPUT,OUTPUT
          real,dimension(3,3)                  :: LATT,BOX
          character*2,dimension(:),allocatable :: ATLABEL
          real,dimension(:,:),allocatable      :: ATCOORD
          real,dimension(3)                    :: ORG
          real,dimension(:,:,:),allocatable    :: GRID
          real                                 :: GRIDSUM = 0.,NEWSUM
          integer                              ::
     &      NGDX,NGDY,NGDZ,I,J,K,NOPT

          print*,'Please specify your option:'
          print*,'  1. Normalize the integrated value to input value.'
          print*,'  2. Divide the grid data by input value.'
          read*,NOPT
          if (NOPT/=1 .and. NOPT/=2) then
            print*,'Error: Invalid input.'
            stop
          endif
          print*,'Please specify value that 3D data is normalized to: '
          read*,NEWSUM

          call read_3dxsf(INPUT,LATT,ATLABEL,ATCOORD,ORG,BOX,GRID)
          NGDX = size(GRID,dim=1)
          NGDY = size(GRID,dim=2)
          NGDZ = size(GRID,dim=3)

          if (NOPT == 1) then
            do K = 1,NGDZ
              do J = 1,NGDY
                do I = 1,NGDX
                  GRIDSUM = GRIDSUM + GRID(I,J,K)
                end do
              end do
            end do
            do K = 1,NGDZ
              do J = 1,NGDY
                do I = 1,NGDX
                  GRID(I,J,K) = GRID(I,J,K) / GRIDSUM * NEWSUM
                end do
              end do
            end do
          end if

          if (NOPT == 2) then
            do K = 1,NGDZ
              do J = 1,NGDY
                do I = 1,NGDX
                  GRID(I,J,K) = GRID(I,J,K) / NEWSUM
                end do
              end do
            end do
          endif
          call write_3dxsf(OUTPUT,LATT,ATLABEL,ATCOORD,ORG,BOX,GRID)
        end subroutine option4
      end module option
