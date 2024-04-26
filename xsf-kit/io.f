      module io
        implicit none
        private
!       NGDM     : Maximum number of grids along x y and z (NGDM*NGDM*NGDM)
!       NATM     : Maximum number of atoms
!       NGDX/Y/Z : Actual rid points along x/y/z
!       NAT      : Actual number of atoms
        integer,parameter :: NGDM = 2000
        integer,parameter :: NATM = 1000
        integer           :: NGDX,NGDY,NGDZ,NAT,NGD,NGDAVG

        public :: read_3dxsf,write_1dtxt,write_3dxsf,write_scatter
!    &    ,NGDX,NGDY,NGDZ,NAT,NGDAVG

        contains
        subroutine read_3dxsf(XSF,LATT,ATLABEL,ATCOORD,ORG,BOX,GRID)
!         Read XCrySDen xsf format 3D grid data.
!         XSF     : Name of XCrySDen XSF file
!         LATT    : 3*3 lattice matrix, in Angstrom
!         ATLABEL : NAT*1 character list of atom species
!         ATCOORD : 3*NAT list of atomic coordinates, in Angstrom
!         ORG     : Origin of data box, in Angstrom
!         BOX     : Matrix of 3D data box, in Angstrom
!         GRID    : NGDX*NGDY*NGDZ grid data. For QE, e/Bohr^3. For VASP, e
          character(len=80),intent(in)                     :: XSF
          character(len=100)                               :: HEADER
          integer                                          :: I,J,K
          real,dimension(3,3),intent(out)                  :: LATT,BOX
          character*2,dimension(:),allocatable,intent(out) :: ATLABEL
          real,dimension(:,:),allocatable,intent(out)      :: ATCOORD
          real,dimension(3),intent(out)                    :: ORG
          real,dimension(:,:,:),allocatable,intent(out)    :: GRID

          open(10,file=XSF,status='OLD',err=1000)
          read(10,'(A)',err=1000,end=1000) HEADER
          do while(index(HEADER,'PRIMVEC') == 0)
            read(10,'(A)',err=1000,end=1000) HEADER
          enddo
!         Read lattice matrix
          read(10,*,err=1000,end=1000) ((LATT(I,J),I=1,3),J=1,3)

!         Read atomic species and coords
          read(10,'(A)',err=1000,end=1000) HEADER
          do while(index(HEADER,'PRIMCOORD') == 0)
            read(10,'(A)',err=1000,end=1000) HEADER
          enddo
          read(10,*,err=1000,end=1000) NAT,I
          if (NAT > NATM) then
            print*,'Too many atoms. Maximum atoms allowed: ',NATM,
     &        '. Exiting'
            stop
          endif
          allocate(ATLABEL(NAT),ATCOORD(3,NAT))
          read(10,*,err=1000,end=1000) (ATLABEL(I),
     &      ATCOORD(1,I),ATCOORD(2,I),ATCOORD(3,I),I=1,NAT)

          print*,'3D geometry data read'

!         Read box where 3D data is defined
          read(10,'(A)',err=1000,end=1000) HEADER
          do while(index(HEADER,'DATAGRID_3D_UNKNOWN') == 0)
            read(10,'(A)',err=1000,end=1000) HEADER
          end do

!         Read 3D grid data
          read(10,*,err=1000,end=1000) NGDX,NGDY,NGDZ
          NGD = NGDX * NGDY * NGDZ
          if (NGDX > NGDM .or. NGDY > NGDM .or. NGDZ > NGDM) then
            print*,'Grid too large. Maximum grid numbers along X/Y/Z: ',
     &        NGDM,'. Exiting.'
            stop
          endif
          read(10,*,err=1000,end=1000) (ORG(I),I=1,3)
          read(10,*,err=1000,end=1000) ((BOX(I,J),I=1,3),J=1,3)
          allocate(GRID(NGDX,NGDY,NGDZ))
          read(10,*,err=1000,end=1001)
     &      (((GRID(I,J,K),I=1,NGDX),J=1,NGDY),K=1,NGDZ)
          print*,'3D grid data read'
          close(10)
          return

1001      print*,'3D grid data read - but probably abnormal',
     &      'termination.';close(10);return
1000      print*,'Error opening or reading the 3D XSF file ',XSF;stop
        end subroutine read_3dxsf
!----
        subroutine write_1dtxt(TXTOUT,AREA,SHIFTA,DIST,AVG1D,INT1D)
!         Write 1D planar-averaged data into a txt file
!         TXTOUT  : Output file name
!         AREA    : Cross-sectional area of averaged plane
!         SHIFTA  : Shifting length in Angstrom
!         DIST    : Displacement (x axis), in Å
!         AVG1D   : Averaged data (y axis). Cross section area normalized to 1
!         INT1D   : Integrated data (y axis). Cross section area normalized to 1
          character(len=80),intent(in) :: TXTOUT
          real,intent(in)              :: AREA,SHIFTA
          real,dimension(:),intent(in) :: DIST,AVG1D,INT1D
          integer                      :: I
          real                         :: MIDDIST

          NGDAVG = size(DIST)

          open(20,file=TXTOUT)
          write(20,'(A24,I12)') '# N Points            = ',NGDAVG
          write(20,200) '# Step in Å           = ',DIST(2)-DIST(1)
          write(20,200) '# Cross-section (Å^2) = ',AREA
          write(20,200) '# Shift of box (Å)    = ',SHIFTA
          write(20,'(A16,4X,A16,4X,A16,4X,A16)')
     &    '#       xAVG(Å)','yAVG(Å^-3)','xINT(Å)','yINT(Å^-2)'

          MIDDIST = (DIST(2) - DIST(1)) / 2
          do I=1,NGDAVG-1
            write(20,202) DIST(I),AVG1D(I),DIST(I)+MIDDIST,INT1D(I)
          end do
          write(20,'(F15.8,4X,E15.8)') DIST(NGDAVG),AVG1D(NGDAVG)
200       format(A25,F12.6)
202       format(F15.8,4X,E15.8,4X,F15.8,4X,E15.8)
          write(20,'(/)')
          close(20)
        end subroutine write_1dtxt
!----
        subroutine write_3dxsf(XSFOUT,LATT,ATLABEL,ATCOORD,ORG,BOX,GRID)
!         Write 3D grid data into a XCrySDen xsf file
!         XSFOUT : Output file name
          character(len=80),intent(in)             :: XSFOUT
          real,dimension(3,3),intent(in)           :: LATT,BOX
          character*2,dimension(:),intent(in)      :: ATLABEL
          real,dimension(:,:),intent(in)           :: ATCOORD
          real,dimension(3),intent(in)             :: ORG
          real,dimension(:,:,:),intent(in)         :: GRID
          integer                                  :: I,J,K

          write(0,"(1X,'Writting into file ',A, '...')") trim(XSFOUT)
          NAT = size(ATLABEL,dim=1)
          NGDX = size(GRID,dim=1)
          NGDY = size(GRID,dim=2)
          NGDZ = size(GRID,dim=3)

          open(21,file=XSFOUT)
          write(21,'(A8)') ' CRYSTAL'
          write(21,'(A8)') ' PRIMVEC'
          write(21,'(3F15.9)') ((LATT(I,J),I=1,3),J=1,3)
          write(21,'(A10)') ' PRIMCOORD'
          write(21,'(2I12)') NAT,1
          do J = 1,NAT
            write(21,'(A6,3F15.9)') ATLABEL(J),ATCOORD(1:3,J)
          enddo
          write(21,'(A)') 'BEGIN_BLOCK_DATAGRID_3D'
          write(21,'(A)') '3D_XSFDATA'
          write(21,'(A)') 'BEGIN_DATAGRID_3D_UNKNOWN'
          write(21,'(3I12)') NGDX,NGDY,NGDZ
          write(21,'(3F10.6)') (ORG(I),I=1,3)
          write(21,'(3F12.6)') ((BOX(I,J),I=1,3),J=1,3)
          NGD = 0
          do K = 1,NGDZ
            do J = 1,NGDY
              do I = 1,NGDX
                NGD = NGD + 1
                if (mod(NGD,5) == 0) then
                  write(21,'(E14.6)') GRID(I,J,K)
                else
                  write(21,'(E14.6,$)') GRID(I,J,K)
                endif
              end do
            end do
          end do
          write(21,'(A)') 'END_DATAGRID_3D'
          write(21,'(A,/,/)') 'END_BLOCK_DATAGRID_3D'

          close(21)
        end subroutine write_3dxsf
!----
        subroutine write_scatter(SCTOUT,INPUT1,GRID1,INPUT2,GRID2)
!         Write scatter data of grid data correlation into a 1D txt file
!         SCTOUT : Output file name
          character(len=80),intent(in)     :: SCTOUT,INPUT1,INPUT2
          real,dimension(:,:,:),intent(in) :: GRID1,GRID2
          integer                          :: I,J,K,NGDX,NGDY,NGDZ

          NGDX = size(GRID1,dim=1)
          NGDY = size(GRID1,dim=2)
          NGDZ = size(GRID1,dim=3)

          open(22,file=SCTOUT)
          write(22,203) '# x axis:      ',trim(INPUT1)
          write(22,203) '# y axis:      ',trim(INPUT2)
          write(22,203) '#       x value','        y value'
          do K = 1,NGDZ
            do J = 1,NGDY
              do I = 1,NGDX
                write(22,204) GRID1(I,J,K),GRID2(I,J,K)
              end do
            end do
          end do
          close(22)
203       format(A15,4X,A)
204       format(E15.8,4X,E15.8)
        end subroutine write_scatter
      end module io
