# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( Resultant_enumeration_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0005 NEW )
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  # uncomment to use lrslib
  include_directories (BEFORE ../external/lrslib-042c/include)
  
  include_directories (BEFORE ../external/triangulation/include)
  include_directories (BEFORE ../external/extreme_points_d/include)
  include_directories (BEFORE ../external/kernel_d/include)
  
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../patches/include)
  # uncomment when using LinBox
  #include_directories (BEFORE /home/vissarion/src/LinBox/include)
  #include_directories (BEFORE /home/vissarion/src/linbox-1.2.0)
  #include_directories (BEFORE /home/vissarion/src/fflas-ffpack-1.4.1)
  # uncomment when using Eigen
  #include_directories (BEFORE /home/luis/eigen-eigen-3.0.2)

  #add_definitions(${CMAKE_CXX_FLAGS} "-DLinBoxSrcOnly")
  add_definitions(${CMAKE_CXX_FLAGS} "-DPRINT_INFO")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_LINBOX_DET")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_CGAL_DET")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_EIGEN_DET")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DHASH_STATISTICS")
  add_definitions(${CMAKE_CXX_FLAGS} "-DLOG_DET_TIME")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HASHED_DETERMINANTS")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_KERNEL_ORIENTATION")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_KERNEL_ORIENTATION_FOR_CONVEX_HULL")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_GAUSSIAN_ELIMINATION")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_CLEAR_DET_HASH")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_ONLY_CAYLEY_DET_HASH")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_EXTREME_SPECIALIZED_POINTS_ONLY")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_ORIENTATION_DET")
  add_definitions(${CMAKE_CXX_FLAGS} "-DCOMPUTE_VOL")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DRESTRICTED_RES")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_LRSLIB")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DRANDOM_RES")
  
  #link_directories(/home/vissarion/workspace/Triangulation_hdoh/Triangulation/respol/external/lrslib-042c/include)
  #target_link_libraries(../external/lrslib-042c/include/lrsmp.c)
  
  # MAPLE
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_MAPLE_CONVEX_HULL")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DMAPLE_EXECUTABLE=/workspace/maple13/bin/maple")

  # uncomment when using LinBox
  #link_libraries("blas")
  #link_libraries("lapack")
  #link_libraries("givaro")
  target_link_libraries("/usr/lib/liblrs.so")
  create_single_source_cgal_program("res_enum_d.cpp" )
  #ADD_EXECUTABLE(res_enum_d
  #               ../external/lrslib-042c/include/lrsmp.c
  #               ../external/lrslib-042c/include/lrslib.c
  #               res_enum_d.cpp
  #               )
else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
