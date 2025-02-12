# Install script for directory: /home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build/Src/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX" TYPE FILE FILES
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build/lib/cmake/AMReX/AMReXConfig.cmake"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build/lib/cmake/AMReX/AMReXConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build/Src/libamrex_3d.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ccse-mpi.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Math.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Algorithm.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Any.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Array.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BlockMutex.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Enum.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuComplex.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_SmallMatrix.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ConstexprFor.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Vector.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_TableData.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Tuple.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_TypeList.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Demangle.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Exception.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Extension.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_PODVector.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ParmParse.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Functional.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Stack.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_String.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Utility.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FileSystem.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ValLocPair.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Reduce.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Scan.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Partition.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Morton.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Random.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_RandomEngine.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BLassert.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ArrayLim.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_REAL.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_INT.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_CONSTANTS.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_SPACE.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_DistributionMapping.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ParallelDescriptor.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_OpenMP.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ParallelReduce.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ForkJoin.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ParallelContext.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_VisMFBuffer.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_VisMF.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_AsyncOut.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BackgroundThread.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Arena.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BArena.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_CArena.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_PArena.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_DataAllocator.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BLProfiler.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BLBackTrace.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BLFort.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_NFiles.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_parstream.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ANSIEscCode.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FabConv.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FPC.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_VectorIO.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Print.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_IntConv.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_IOFormat.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Box.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BoxIterator.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Dim3.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_IntVect.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_IndexType.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Loop.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Loop.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Orientation.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Periodicity.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_RealBox.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_RealVect.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BoxList.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BoxArray.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BoxDomain.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FArrayBox.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_IArrayBox.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BaseFab.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Array4.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MakeType.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_TypeTraits.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FabDataType.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FabFactory.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BaseFabUtility.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MultiFab.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MFCopyDescriptor.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_iMultiFab.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FabArrayBase.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MFIter.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FabArray.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FACopyDescriptor.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FabArrayCommI.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FBI.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_PCI.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FabArrayUtility.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_LayoutData.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_CoordSys.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_COORDSYS_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_COORDSYS_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Geometry.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MultiFabUtil.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MultiFabUtilI.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MultiFabUtil_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MultiFabUtil_nd_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MultiFabUtil_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BCRec.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_PhysBCFunct.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BCUtil.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BC_TYPES.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FilCC_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FilCC_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FilFC_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FilFC_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FilND_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_NonLocalBC.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_NonLocalBCImpl.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_PlotFileUtil.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_PlotFileDataImpl.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_FEIntegrator.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_IntegratorBase.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_RKIntegrator.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_TimeIntegrator.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_RungeKutta.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Gpu.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuQualifiers.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuKernelInfo.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuPrint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuAssert.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuTypes.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuControl.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunch.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunch.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunchGlobal.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunchMacrosG.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunchMacrosG.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunchMacrosC.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunchMacrosC.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunchFunctsG.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuLaunchFunctsC.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuError.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuDevice.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuBuffer.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuAtomic.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuUtility.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuAsyncArray.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuElixir.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuMemory.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuRange.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuReduce.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuAllocators.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_GpuContainers.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MFParallelFor.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MFParallelForC.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MFParallelForG.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_TagParallelFor.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_CTOParallelForImpl.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_ParReduce.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_CudaGraph.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Machine.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MemPool.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/AMReX_Parser.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/AMReX_Parser_Exe.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/AMReX_Parser_Y.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/amrex_parser.lex.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/amrex_parser.tab.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/AMReX_IParser.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/AMReX_IParser_Exe.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/AMReX_IParser_Y.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/amrex_iparser.lex.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/Parser/amrex_iparser.tab.nolint.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_LUSolver.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_Slopes_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_BaseFwd.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Base/AMReX_MPMD.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_FabSet.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_BndryRegister.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_Mask.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_MultiMask.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_BndryData.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_BoundCond.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_InterpBndryData.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_LO_BCTYPES.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_InterpBndryData_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_InterpBndryData_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_LOUtil_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_YAFluxRegister.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_YAFluxRegister_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_YAFluxRegister_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_BoundaryFwd.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Boundary/AMReX_EdgeFluxRegister.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_AmrCore.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_Cluster.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_ErrorList.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_FillPatchUtil.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_FillPatchUtil_I.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_FillPatcher.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_FluxRegister.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_InterpBase.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_MFInterpolater.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_Interpolater.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_TagBox.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_AmrMesh.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_FluxReg_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_FluxReg_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_Interp_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_Interp_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_MFInterp_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_MFInterp_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_InterpFaceRegister.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_InterpFaceReg_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_InterpFaceReg_3D_C.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_AmrCoreFwd.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_AmrParGDB.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/AmrCore/AMReX_AmrParticles.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_LevelBld.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_Amr.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_AmrLevel.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_Derive.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_StateData.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_PROB_AMR_F.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_StateDescriptor.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_AuxBoundaryData.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_Extrapolater.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_extrapolater_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_extrapolater_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Amr/AMReX_AmrFwd.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLMG.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLMG_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLMG_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLMGBndry.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLLinOp.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLLinOp_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLCellLinOp.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeLinOp_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLCellABecLap.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLCellABecLap_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLCellABecLap_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLCGSolver.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_PCGSolver.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLABecLaplacian.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLABecLap_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLABecLap_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLALaplacian.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLALap_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLALap_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLPoisson.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLPoisson_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLPoisson_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_GMRES.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_GMRES_MLMG.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_GMRES_MV.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_Smoother_MV.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_Algebra.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_AlgPartition.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_AlgVector.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_SpMatrix.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/AMReX_SpMV.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLMG_2D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLPoisson_2D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLALap_2D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLCurlCurl.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLCurlCurl_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLEBNodeFDLaplacian.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLEBNodeFDLap_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLEBNodeFDLap_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeTensorLaplacian.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeTensorLap_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeTensorLap_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeABecLaplacian.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeABecLap_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeABecLap_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeLaplacian.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeLap_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLNodeLap_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLTensorOp.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLTensor_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/MLMG/AMReX_MLTensor_3D_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/OpenBC/AMReX_OpenBC.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/LinearSolvers/OpenBC/AMReX_OpenBC_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_Particles.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleContainer.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_SparseBins.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParGDB.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_Particle_mod_K.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_TracerParticles.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_NeighborParticles.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_NeighborParticlesI.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_NeighborList.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_Particle.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleInit.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleContainerI.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParIter.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleMPIUtil.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleUtil.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_StructOfArrays.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ArrayOfStructs.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleTile.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_MakeParticle.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_NeighborParticlesCPUImpl.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_NeighborParticlesGPUImpl.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleBufferMap.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleCommunication.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleInterpolators.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleReduce.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleMesh.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleLocator.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleIO.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_DenseBins.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_BinIterator.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleTransformation.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_WriteBinaryParticleData.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleContainerBase.H"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Src/Particle/AMReX_ParticleArray.H"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets.cmake"
         "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build/CMakeFiles/Export/2260e541ece776bcef17e59de6c71ec8/AMReXTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX" TYPE FILE FILES "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build/CMakeFiles/Export/2260e541ece776bcef17e59de6c71ec8/AMReXTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX" TYPE FILE FILES "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-build/CMakeFiles/Export/2260e541ece776bcef17e59de6c71ec8/AMReXTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(CREATE_LINK
           libamrex_3d.a
           "/usr/local/lib/libamrex.a"
           COPY_ON_ERROR SYMBOLIC)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/amrex" TYPE DIRECTORY FILES
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Tools/C_scripts"
    "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Tools/typechecker"
    USE_SOURCE_PERMISSIONS)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/AMReX/AMReXCMakeModules" TYPE DIRECTORY FILES "/home/adamdarx/bnuGRMHD/build/_deps/amrex_code-src/Tools/CMake/" USE_SOURCE_PERMISSIONS)
endif()

