      module option
        implicit none
        private
        integer,parameter :: NGDM = 1000
        integer,parameter :: NATM = 1000

        public :: optionA1,optionA2,optionB1,optionB2,optionB3,optionB4

        contains
        subroutine optionA1(INPUT,OUTPUT)
!         Get planar averaged data for a single XDF file
          use io
          use operation

          character(len=80),intent(in)         :: INPUT,OUTPUT
          integer                              :: AVGVEC
          character(len=10)                    :: SHIFTC
          real                                 :: SHIFT=0.0,AREA,SHIFTA
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
     &(Unit: fractional unit of data grid): '
          print*,"  Use 'no' for no shift."
          read*,SHIFTC
          if (trim(SHIFTC) == 'no') then
            SHIFT= 0.0
          else
            read(SHIFTC,'(f10.6)') SHIFT
          endif
          if (SHIFT >= 1 .or. SHIFT <= -1) then
            print*,'WARNING: SHIFT should between -1 and 1. Integer part
     & will be removed'
            SHIFT = SHIFT - int(SHIFT)
          endif

          call read_3dxsf(INPUT,LATT,ATLABEL,ATCOORD,ORG,BOX,GRID)
          call planar_avg(LATT,ATCOORD,ORG,BOX,GRID,AVGVEC,SHIFT,
     &                    AREA,DIST,AVG1D,INT1D)
          SHIFTA = (BOX(1,AVGVEC)**2
     &            + BOX(2,AVGVEC)**2 + BOX(3,AVGVEC)**2)**0.5 * SHIFT
          call write_1dtxt(OUTPUT,AREA,SHIFTA,DIST,AVG1D,INT1D)
          deallocate(ATLABEL,ATCOORD,GRID,DIST,AVG1D,INT1D)
        end subroutine optionA1
!----
        subroutine optionA2
!----     Normalize the integrated 3D grid data
          use io

          character(len=80)                    :: INPUT,OUTPUT
          real,dimension(3,3)                  :: LATT,BOX
          character*2,dimension(:),allocatable :: ATLABEL
          real,dimension(:,:),allocatable      :: ATCOORD
          real,dimension(3)                    :: ORG
          real,dimension(:,:,:),allocatable    :: GRID
          real                                 :: GRIDSUM = 0.,NEWSUM
          integer                              ::
     &      NGDX,NGDY,NGDZ,I,J,K,NOPT

          print*,'Please specify the name of 3D XSF input: '
          read*,INPUT
          print*,'Please specify the name of 3D XSF output: '
          read*,OUTPUT
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
          deallocate(ATLABEL,ATCOORD,GRID)
        end subroutine optionA2
!----
        subroutine optionB1(INPUT0,OUTPUT)
!         Get difference of 3D data
          use io
          use operation

          character(len=80),intent(in) :: INPUT0,OUTPUT
          character(len=80)   :: INPUT1
          real,dimension(3,3) :: LATT0,LATT1,BOX0,BOX1
          integer             :: NINPUT,I
          character*2,dimension(:),allocatable      :: ATLABEL0,ATLABEL1
          real,dimension(:,:),allocatable           :: ATCOORD0,ATCOORD1
          real,dimension(3)                         :: ORG0,ORG1
          real,dimension(:,:,:),allocatable         :: GRID0,GRID1

          print*,'Please specify the number of other 3D XSF files: '
          read*,NINPUT
          call
     &      read_3dxsf(INPUT0,LATT0,ATLABEL0,ATCOORD0,ORG0,BOX0,GRID0)

          do I = 1,NINPUT
            print*,'File name: '
            read*,INPUT1
            call
     &        read_3dxsf(INPUT1,LATT1,ATLABEL1,ATCOORD1,ORG1,BOX1,GRID1)
            call compare_grid(ORG0,BOX0,GRID0,ORG1,BOX1,GRID1)
            GRID0 = GRID0 - GRID1
          end do
          call
     &      write_3dxsf(OUTPUT,LATT0,ATLABEL0,ATCOORD0,ORG0,BOX0,GRID0)
          deallocate(ATLABEL0,ATLABEL1,ATCOORD0,ATCOORD1,GRID0,GRID1)
        end subroutine optionB1
!----
        subroutine optionB2
!----     Get difference of multiple 3D grid data and 1D line profile
          character(len=80) :: INPUT0,OUT3,OUT1

          print*,'Please specify the name of main 3D XSF file: '
          read*,INPUT0
          print*,'Please specify the name of 3D XSF output: '
          read*,OUT3
          print*,'Please specify the name of 1D line profile file: '
          read*,OUT1

          call optionB1(INPUT0,OUT3)
          call optionA1(OUT3,OUT1)
        end subroutine optionB2
!----
        subroutine optionB3
!----     Get correlation of 2 3D xsf data
          use io
          use operation

          character(len=80)                    :: INPUT1,INPUT2,OUTPUT
          real,dimension(3,3)                  :: LATT,BOX1,BOX2
          character*2,dimension(:),allocatable :: ATLABEL
          real,dimension(:,:),allocatable      :: ATCOORD
          real,dimension(3)                    :: ORG1,ORG2
          real,dimension(:,:,:),allocatable    :: GRID1,GRID2
          integer                              :: I,J,K

          print*,'NOTE: This module analysis correlation of 2 xsf data.'
          print*,'      A txt file for scatter plotting is generated.'
          print*,'Please specify the 3D XSF file as x axis:'
          read*,INPUT1
          print*,'Please specify the 3D XSF file as y axis:'
          read*,INPUT2
          print*,'Please specify the 2D plot file name:'
          read*,OUTPUT

          ! Read and compare dimensions of grid data
          call read_3dxsf(INPUT1,LATT,ATLABEL,ATCOORD,ORG1,BOX1,GRID1)
          call read_3dxsf(INPUT2,LATT,ATLABEL,ATCOORD,ORG2,BOX2,GRID2)
          call compare_grid(ORG1,BOX1,GRID1,ORG2,BOX2,GRID2)
          ! write file
          call write_scatter(OUTPUT,INPUT1,GRID1,INPUT2,GRID2)
          deallocate(ATLABEL,ATCOORD,GRID1,GRID2)
        end subroutine optionB3
!----
        subroutine optionB4
!----     Map isosurfaces of one grid over another.
          use io
          use operation

          character(len=80)    :: INPUT1,INPUT2,OUTPUT
          character(len=20)    :: VTXT1,VTXT2,IVAL
          character            :: C
          real                 :: VMIN,VMAX,VTMP,IMIN,IMAX
          logical              :: ISABS=.false.
          integer              :: I,NPM=0
          integer,dimension(2) :: LOCPM
          real,dimension(3,3)                  :: LATT,BOX1,BOX2
          character*2,dimension(:),allocatable :: ATLABEL
          real,dimension(:,:),allocatable      :: ATCOORD
          real,dimension(3)                    :: ORG1,ORG2
          real,dimension(:,:,:),allocatable    :: GRID1,GRID2

          write(0,"(
     & 'Please specify the reference XSF file for color coding:')")
          read*,INPUT1
          write(0,"(
     & 'Please specify the mapped XSF file for isosurface')")
          read*,INPUT2
          write(0, "('Please specify the range of reference data.',
     & ' Data in mapped grid out of this range will be removed')")
          print*,'  Format: VMIN VMAX. For absolute values, use +-'
          read*,VTXT1,VTXT2
          print*,'Please indicate the isosurface value and direction:'
          print*,'  Values beyond this range will be removed.'
          print*,'  Format: <=IVAL, >=IVAL, <=+-IVAL or >=+-IVAL.'
          read*,IVAL
          print*,'Please specify the 3D XSF output:'
          read*,OUTPUT

!         Get values
          ISABS = .false.
          if (VTXT1(1:2) == '+-') then
            ISABS = .true.
            read(VTXT1(3:),'(E14.6)') VMIN
          else
            read(VTXT1,'(E14.6)') VMIN
          end if
          if (VTXT2(1:2) == '+-') then
            ISABS = .true.
            read(VTXT2(3:),'(E14.6)') VMAX
          else
            read(VTXT2,'(E14.6)') VMAX
          end if
          if (VMIN > VMAX) then
            VTMP = VMIN
            VMIN = VMAX
            VMAX = VTMP
          end if
          write(0, "(/,1X,'Range: ',E14.6,' to ',E14.6)") VMIN,VMAX
          if (ISABS) then
            write(0, "(1X,'Use absolute values...')")
          else
            write(0, "(1X,'Do not use absolute values...')")
          end if
!         Map the data
          call read_3dxsf(INPUT1,LATT,ATLABEL,ATCOORD,ORG1,BOX1,GRID1) ! ref grid
          call read_3dxsf(INPUT2,LATT,ATLABEL,ATCOORD,ORG2,BOX2,GRID2) ! isosurface grid
          call compare_grid(ORG1,BOX1,GRID1,ORG2,BOX2,GRID2)

          if (IVAL(1:2) == '<=') then
            if (IVAL(3:4) == '+-') then
              read(IVAL(5:),'(E14.6)') IMAX
              IMIN = -IMAX
            else
              read(IVAL(3:),'(E14.6)') IMAX
              IMIN = minval(GRID2)
            end if
          else if (IVAL(1:2) == '>=') then
              if (IVAL(3:4) == '+-') then
                read(IVAL(5:),'(E14.6)') IMIN ! 2 separate ranges. Let map_grid to get it
                IMAX = -IMIN
              else
                read(IVAL(3:),'(E14.6)') IMIN
                IMAX = maxval(GRID2)
              end if
          else
            print*,'ERROR: Unknown isosurface value range.'
          end if

          call map_grid(GRID1,VMIN,VMAX,ISABS,IMIN,IMAX,GRID2)
          ! write file
          call write_3dxsf(OUTPUT,LATT,ATLABEL,ATCOORD,ORG2,BOX2,GRID2)
          deallocate(ATLABEL,ATCOORD,GRID1,GRID2)
        end subroutine optionB4
      end module option
