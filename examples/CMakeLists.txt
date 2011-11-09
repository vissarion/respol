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

  option ( RESPOL_USE_LRS "use LRS to compute convex hulls" FALSE )
  #option ( RESPOL_USE_MAPLE "use Maple to check convex hulls" FALSE )
  #set ( RESPOL_USE_MAPLE_EXEC "/opt/maple13/bin/maple" CACHE FILEPATH
  #      "path to Maple executable" )
  #option ( RESPOL_USE_LINBOX "use LinBox to compute determinants" FALSE )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE ../external/triangulation/include)
  include_directories (BEFORE ../external/spatial_sorting/include)
  include_directories (BEFORE ../external/extreme_points_d/include)
  include_directories (BEFORE ../external/kernel_d/include)

  include_directories (BEFORE ../include)
  include_directories (BEFORE ../patches/include)

  add_definitions(${CMAKE_CXX_FLAGS} "-DPRINT_INFO")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_CGAL_DET")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_EIGEN_DET")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DHASH_STATISTICS")
  add_definitions(${CMAKE_CXX_FLAGS} "-DLOG_DET_TIME")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HASHED_DETERMINANTS")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_SORTED_INDICES")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_KERNEL_ORIENTATION")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_KERNEL_ORIENTATION_FOR_CONVEX_HULL")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_GAUSSIAN_ELIMINATION")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_CLEAR_DET_HASH")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_ONLY_CAYLEY_DET_HASH")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_EXTREME_SPECIALIZED_POINTS_ONLY")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_ORIENTATION_DET")
  add_definitions(${CMAKE_CXX_FLAGS} "-DCOMPUTE_VOL")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DRESTRICTED_RES")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DRANDOM_RES")

  if( RESPOL_USE_LINBOX )
    include_directories (BEFORE /home/vissarion/src/LinBox/include)
    include_directories (BEFORE /home/vissarion/src/linbox-1.2.0)
    include_directories (BEFORE /home/vissarion/src/fflas-ffpack-1.4.1)
    add_definitions(${CMAKE_CXX_FLAGS} "-DLinBoxSrcOnly")
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_LINBOX_DET")
    link_libraries("blas")
    link_libraries("lapack")
    link_libraries("givaro")
  endif( RESPOL_USE_LINBOX )

  if( RESPOL_USE_MAPLE )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_MAPLE_CONVEX_HULL")
    add_definitions(${CMAKE_CXX_FLAGS}
                    "-DMAPLE_EXECUTABLE=${RESPOL_USE_MAPLE_EXEC}")
  endif( RESPOL_USE_MAPLE )

  if( RESPOL_USE_LRS )
    message(WARNING "LRS support is buggy!")
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_LRSLIB")
    add_definitions(${CMAKE_CXX_FLAGS} "-DTIMES")
    add_definitions(${CMAKE_CXX_FLAGS} "-DSIGNALS")
    add_definitions(${CMAKE_CXX_FLAGS} "-DGMP")
    if ( ${CMAKE_SYSTEM_PROCESSOR} MATCHES ".*64.*" )
      add_definitions(${CMAKE_CXX_FLAGS} "-DB64")
    endif ( ${CMAKE_SYSTEM_PROCESSOR} MATCHES ".*64.*" )
    include_directories (BEFORE ../external/lrslib-042c/include)
    set ( RESPOL_AUX_PROGS "../external/lrslib-042c/include/lrslib.c"
                           "../external/lrslib-042c/include/lrsgmp.c"
                           "lrs_functions.c")
    create_single_source_cgal_program("test_lrs.cpp" ${RESPOL_AUX_PROGS})
  endif( RESPOL_USE_LRS )

  create_single_source_cgal_program("res_enum_d.cpp" ${RESPOL_AUX_PROGS})

else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
