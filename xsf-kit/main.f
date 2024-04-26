      program xsfkit
!     Process 3D real space grid data in XCrySDen XSF format for various
!     applications. Including:
!     1. Difference among various sets of data, such as differential
!     charge density
!     2. Planar-averaged 1D profile projected along a given lattice
!     vector (a,b,c), such as planar-averaged electrostatic potential
!     3. 3D grid data processing, including normalizing integrated value
!     or value on grid points
!     4. Data analysis, namely plot data set A as x axis and B as y for
!     their correlations
!
!     Note that by default, units reported are consistent with raw data.
!
!      To launch the executable, the user can either:
!      1. Copy the binary executable into work directory
!      2. Specify the full path of input data
!
!     By Spica.Vir
!     ------------------------------------------------------------------
!     Originally edited for VASP5, @NWPU. 5th Apr., 2020
!     ------------------------------------------------------------------
!     Revised for VASP5, @ICL. 25th Mar.; 14th May., 2021
!     ------------------------------------------------------------------
!     Revised for XCrySDen, @ICL. 26th Mar., 2023
!     ------------------------------------------------------------------
!     Update for new functionalities, @ICL, 17th Apr., 2024
!     ------------------------------------------------------------------
        use option

        character(len=10) :: OPTNUM
        character(len=80) :: INPUT,OUTPUT

        print*,'======================================================='
        print*,'                 XSF-KIT BY SPICA.VIR'
        print*,'NOTE: Geometry unit: Å'
        print*,'      Data grid unit: Commmensurate with input'
        print*,'======================================================='
        print*,''
        print*,'A Single XSF data grid'
        print*,'  A1 Planar-averaged line profile and its integration.'
        print*,'  A2 Normalization.'
        print*,'-------------------------------------------------------'
        print*,'B Multiple XSF data grids'
        print*,'  B1 Differences of multiple XSF data grids.'
        print*,'  B2 Data grid difference + line profile / integration.'
        print*,'  B3 Correlation between 2 data grids.'
        print*,'  B4 Map isosurfaces of one grid over another(In Dev).'
        print*,''
        print*,'======================================================='
        print*,'Please enter your option: '
        read*,OPTNUM
        OPTNUM = trim(OPTNUM)

        if (OPTNUM == 'A1') then
          print*,'Please specify the name of 3D XSF file: '
          read*,INPUT
          print*,'Please specify the name of 1D line profile file: '
          read*,OUTPUT
          call optionA1(trim(INPUT),trim(OUTPUT))
        else if (OPTNUM == 'A2') then
          call optionA2
        else if (OPTNUM == 'B1') then
          print*,'Please specify the name of main 3D XSF file: '
          read*,INPUT
          print*,'Please specify the name of 3D XSF output: '
          read*,OUTPUT
          call optionB1(trim(INPUT),trim(OUTPUT))
        else if (OPTNUM == 'B2') then
          call optionB2
        else if (OPTNUM == 'B3') then
          call optionB3
        else if (OPTNUM == 'B4') then
          call optionB4
        else
          print*,'Error: Option not supported. Exiting.'
          stop
        endif
        stop
      end program xsfkit