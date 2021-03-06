# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Resultant_enumeration_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0005 NEW )
endif()

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  # COMPILATION OPTIONS
  ######################################################################

  include ( CMakeDependentOption )
  option ( RESPOL_USE_FLAT_SET "use boost flat set to store normals" OFF )
  option ( RESPOL_USE_HASHED_DET "hash computed determinant minors" ON )
  CMAKE_DEPENDENT_OPTION ( RESPOL_HASH_STATS
                           "print hashed determinants statistics" OFF
                           "RESPOL_USE_HASHED_DET" OFF )
  CMAKE_DEPENDENT_OPTION ( RESPOL_HASH_SORTED_INDICES
                           "used sorted indices in hash" ON
                           "RESPOL_USE_HASHED_DET" OFF )
  CMAKE_DEPENDENT_OPTION ( RESPOL_HASH_CLEAR
                           "clear the hash when memory consumption is high" ON
                           "RESPOL_USE_HASHED_DET" OFF )
  CMAKE_DEPENDENT_OPTION ( RESPOL_HASH_ONLY_CAYLEY
                           "use hash only for determinants in Cayley space" ON
                           "RESPOL_USE_HASHED_DET" OFF )
  option ( RESPOL_VOL "compute the volume of the output polytope" ON )
  CMAKE_DEPENDENT_OPTION ( RESPOL_HASH_DET_TIME
                           "output the determinant computation time" OFF
                           "RESPOL_USE_HASHED_DET" OFF )
  #option ( RESPOL_USE_LINBOX "use LinBox to compute determinants" OFF )
  #option ( RESPOL_USE_MAPLE "use Maple to check convex hulls" OFF )
  #set ( RESPOL_USE_MAPLE_EXEC
  #      "/opt/maple13/bin/maple"
  #      CACHE
  #      FILEPATH
  #      "path to Maple executable" )
  #option ( RESPOL_USE_LRS "use LRS to compute convex hulls" OFF )
  #CMAKE_DEPENDENT_OPTION ( RESPOL_LRS_DEBUG
  #                         "show info about LRS execution" OFF
  #                         "RESPOL_USE_LRS" OFF )
  #option ( RESPOL_USE_EIGEN "use Eigen to compute determinants" OFF )
  option ( RESPOL_EXTREME_SPECIALIZED_POINTS_ONLY "preprocess the input" OFF )
  option ( RESPOL_BUILD_RANDOMIZED "also build randomized init version" OFF )
  option ( RESPOL_BUILD_CONVERTERS "build input file converters" OFF )

  # INCLUDES
  ######################################################################

  include( ${CGAL_USE_FILE} )
  include( CGAL_CreateSingleSourceCGALProgram )
  include_directories (BEFORE ../external/triangulation/include)
  include_directories (BEFORE ../external/spatial_sorting/include)
  include_directories (BEFORE ../external/extreme_points_d/include)
  include_directories (BEFORE ../external/kernel_d/include)
  include_directories (BEFORE ../external)
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../external/leda/incl)
  include_directories (BEFORE ../patches/include)
  
  # OTHER COMPILATION DEFINITIONS
  ######################################################################

  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_ORIENTATION_DET")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_KERNEL_ORIENTATION")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_KERNEL_ORIENTATION_FOR_CONVEX_HULL")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_GAUSSIAN_ELIMINATION")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DRESTRICTED_RES")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DRANDOM_RES")

  # CMAKE OPTIONS PROCESSING
  ######################################################################

  if( RESPOL_USE_FLAT_SET )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_BOOST_FLAT_SET")
  endif( RESPOL_USE_FLAT_SET )

  if( RESPOL_USE_HASHED_DET )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HASHED_DETERMINANTS")
  else( RESPOL_USE_HASHED_DET )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_CGAL_DET")
  endif( RESPOL_USE_HASHED_DET )

  if( RESPOL_HASH_STATS )
    add_definitions(${CMAKE_CXX_FLAGS} "-DHASH_STATISTICS")
  endif( RESPOL_HASH_STATS )

  if( RESPOL_HASH_SORTED_INDICES )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_SORTED_INDICES")
  endif( RESPOL_HASH_SORTED_INDICES )

  if( RESPOL_EXTREME_SPECIALIZED_POINTS_ONLY )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_EXTREME_SPECIALIZED_POINTS_ONLY")
  endif( RESPOL_EXTREME_SPECIALIZED_POINTS_ONLY )

  # Set the default value of RESPOL_HASH_CLEAR_DETS or use the cached value
  # if it exists.
  if( RESPOL_HASH_CLEAR_DETS )
    set ( RESPOL_HASH_CLEAR_DETS ${RESPOL_HASH_CLEAR_DETS} CACHE INTERNAL
          "cached setting of disabled option" )
  else( RESPOL_HASH_CLEAR_DETS )
    set ( RESPOL_HASH_CLEAR_DETS 1000000 CACHE INTERNAL "default setting" )
  endif( RESPOL_HASH_CLEAR_DETS )

  if( RESPOL_HASH_CLEAR )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_CLEAR_DET_HASH")
    # This line makes this option visible in the GUI.
    set ( RESPOL_HASH_CLEAR_DETS ${RESPOL_HASH_CLEAR_DETS} CACHE INTEGER
          "clear the hash when it stores this amount of determinants" FORCE )
    add_definitions(${CMAKE_CXX_FLAGS}
                    "-DCLEAR_DET_NUMBER=${RESPOL_HASH_CLEAR_DETS}")
  endif( RESPOL_HASH_CLEAR )

  if( RESPOL_HASH_ONLY_CAYLEY )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_ONLY_CAYLEY_DET_HASH")
  endif( RESPOL_HASH_ONLY_CAYLEY )

  if ( RESPOL_VOL )
    add_definitions(${CMAKE_CXX_FLAGS} "-DCOMPUTE_VOL")
  endif ( RESPOL_VOL )

  if( RESPOL_HASH_DET_TIME )
    add_definitions(${CMAKE_CXX_FLAGS} "-DLOG_DET_TIME")
  endif( RESPOL_HASH_DET_TIME )

  if( RESPOL_USE_LINBOX )
    # TODO: let the user choose these directories
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

  #if( RESPOL_USE_LRS )
  #  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_LRSLIB")
  #  add_definitions(${CMAKE_CXX_FLAGS} "-DSIGNALS")
  #  add_definitions(${CMAKE_CXX_FLAGS} "-DGMP")
  #  if ( ${CMAKE_SYSTEM_PROCESSOR} MATCHES ".*64.*" )
  #    add_definitions(${CMAKE_CXX_FLAGS} "-DB64")
  #  endif ( ${CMAKE_SYSTEM_PROCESSOR} MATCHES ".*64.*" )
  #  include_directories (BEFORE ../external/lrslib-042c/include)
  #  set ( RESPOL_AUX_PROGS "../external/lrslib-042c/include/lrslib.c"
  #                         "../external/lrslib-042c/include/lrsgmp.c")
  #  create_single_source_cgal_program("test_lrs.cpp" ${RESPOL_AUX_PROGS})
  #endif( RESPOL_USE_LRS )

  if ( RESPOL_USE_EIGEN )
    add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_EIGEN_DET")
  endif ( RESPOL_USE_EIGEN )

  #if( RESPOL_LRS_DEBUG )
  #  add_definitions(${CMAKE_CXX_FLAGS} "-DTIMES")
  #  add_definitions(${CMAKE_CXX_FLAGS} "-DLRS_DEBUG")
  #endif( RESPOL_LRS_DEBUG )

  if( RESPOL_BUILD_CONVERTERS )
    create_single_source_cgal_program("parse_input.cpp")
  endif( RESPOL_BUILD_CONVERTERS )

  if( RESPOL_BUILD_RANDOMIZED )
    create_single_source_cgal_program("res_enum_d_rand.cpp" ${RESPOL_AUX_PROGS})
  endif( RESPOL_BUILD_RANDOMIZED )

  find_library(LEDA NAMES libleda.a PATHS ../external/leda/)

  if (NOT LEDA)

    message(FATAL_ERROR "This program requires the leda library, and will not be compiled.")

  else ()

    message(STATUS "Library leda found: ${LEDA}")
    link_libraries(${LEDA})

  # COMPILATION
  ######################################################################

  create_single_source_cgal_program("res_enum_d.cpp" ${RESPOL_AUX_PROGS})
  create_single_source_cgal_program("test_suite.cpp")
  #create_single_source_cgal_program("tropli_disc.cpp")
  #create_single_source_cgal_program("compute_ch.cpp")

  enable_testing()
  add_test(NAME test_suite COMMAND "test_suite")
  set_tests_properties(test_suite PROPERTIES FAIL_REGULAR_EXPRESSION "FAIL")
  endif()



else()

  message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
