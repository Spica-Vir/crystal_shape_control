      program opt3dxsf
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
!     Revised for normalization, @ICL, 26th Feb., 2024
!     ------------------------------------------------------------------
!     Revised for line integration, @ICL, 19th Mar., 2024
!     ------------------------------------------------------------------
!     Revised for RDG/NCI analysis, @ICL, 16th Apr., 2024
!     ------------------------------------------------------------------
        use option

        integer           :: OPTNUM
        character(len=80) :: INPUT,INPUT2,OUTPUT,OUTPUT2

        print*,'======================================================='
        print*,'                 XSF-KIT BY SPICA.VIR'
        print*,'NOTE: Geometry unit: Å'
        print*,'      Data grid unit: Commmensurate with input'
        print*,'======================================================='
        print*,
     &    '1. Planar-averaged line profile and integration ',
     &    'of 3D XSF data.'
        print*,'2. 3D XSF data differences of multiple files.'
        print*,'3. 3D data difference + line profile and integration.'
        print*,'4. 3D data normalization.'
        print*,'5. Correlation between 2 3D grid data.'
        print*,'Please enter your option: '
        read*,OPTNUM

        if (OPTNUM == 1) then
          print*,'Please specify the name of 3D XSF file: '
          read*,INPUT
          print*,'Please specify the name of 1D line profile file: '
          read*,OUTPUT
          call option1(INPUT,OUTPUT)
        else if (OPTNUM == 2) then
          print*,'Please specify the name of main 3D XSF file: '
          read*,INPUT
          print*,'Please specify the name of 3D XSF output: '
          read*,OUTPUT
          call option2(INPUT,OUTPUT)
        else if (OPTNUM == 3) then
          print*,'Please specify the name of main 3D XSF file: '
          read*,INPUT
          print*,'Please specify the name of 3D XSF output: '
          read*,OUTPUT
          print*,'Please specify the name of 1D line profile file: '
          read*,OUTPUT2
          call option3(INPUT,OUTPUT,OUTPUT2)
        else if (OPTNUM == 4) then
          print*,'Please specify the name of 3D XSF input: '
          read*,INPUT
          print*,'Please specify the name of 3D XSF output: '
          read*,OUTPUT
          call option4(INPUT,OUTPUT)
        else if (OPTNUM == 5) then
          print*,'NOTE: This module analysis correlation of 2 xsf data.'
          print*,'      A txt file for scatter plotting is generated.'
          print*,'Please specify the 3D XSF input as x axis:'
          read*,INPUT
          print*,'Please specify the 3D XSF input as y axis:'
          read*,INPUT2
          print*,'Please specify the 2D plot file name:'
          read*,OUTPUT
          call option5(INPUT,INPUT2,OUTPUT)
        else
          print*,'Error: Option not supported. Exiting.'
          stop
        endif
        stop
      end program opt3dxsf