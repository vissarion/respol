# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( Resultant_enumeration_example )

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
 
if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  
  include( ${CGAL_USE_FILE} )
	
  include( CGAL_CreateSingleSourceCGALProgram )

	# uncomment to use hacked kernel for hashed determinants
  include_directories (BEFORE ../include)
  #include_directories (BEFORE /home/vissarion/src/LinBox/include)
  #include_directories (BEFORE /home/vissarion/src/linbox-1.2.0)
  #include_directories (BEFORE /home/vissarion/src/fflas-ffpack-1.4.1)
  
  #add_definitions(${CMAKE_CXX_FLAGS} "-DLinBoxSrcOnly")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DPRINT_INFO")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_LINBOX_DET")
  #add_definitions(${CMAKE_CXX_FLAGS} "-DHASH_STATISTICS")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HASHED_DETERMINANTS")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_HACKED_KERNEL_ORIENTATION")

  #link_libraries("lapack")
  #link_libraries("givaro")
    
  create_single_source_cgal_program("res_enum_d.cpp" )
	
else()
  
  message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

