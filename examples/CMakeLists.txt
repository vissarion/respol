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

  include_directories (BEFORE ../include)
  include_directories (BEFORE /home/vissarion/src/LinBox/include)
  include_directories (BEFORE /home/vissarion/src/linbox-1.2.0)
  include_directories (BEFORE /home/vissarion/src/fflas-ffpack-1.4.1)
  
  add_definitions(${CMAKE_CXX_FLAGS} "-DLinBoxSrcOnly")
  add_definitions(${CMAKE_CXX_FLAGS} "-DUSE_LINBOX_DET")
  add_definitions(${CMAKE_CXX_FLAGS} "-DHASH_STATISTICS")
  link_libraries("lapack")
  link_libraries("givaro")
    
  #create_single_source_cgal_program( "res_enum_2.cpp" )
  #create_single_source_cgal_program( "res_enum_3_fast.cpp" )
  create_single_source_cgal_program( "res_enum_d.cpp" )
  #create_single_source_cgal_program( "res_enum_d_horn.cpp" )
  #create_single_source_cgal_program( "sec_enum_3_fast.cpp" )
  #create_single_source_cgal_program( "my_quickhull_3.cpp" )
  #create_single_source_cgal_program( "my_quickhull_d.cpp" )
  #create_single_source_cgal_program( "ch_topcom_data.cpp" )

	#target_link_libraries(res_enum_2 "/usr/local/lib/libTOPCOM.a" )
	#target_link_libraries(res_enum_3_fast "/usr/local/lib/libTOPCOM.a" )
	#target_link_libraries(res_enum_d_slow "/usr/local/lib/libTOPCOM.a" )
	#target_link_libraries(res_enum_d_horn "/usr/local/lib/libTOPCOM.a" )
	#target_link_libraries(sec_enum_3_fast "/usr/local/lib/libTOPCOM.a" )
	#target_link_libraries(ch_topcom_data "/usr/local/lib/libTOPCOM.a" )
	
	#target_link_libraries(res_enum_2 "/usr/local/lib/libCHECKREG.a" )
	
else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

