cmake_minimum_required(VERSION 2.8.6) # not sure, simply copied from UWLCM
project(UWLCM_plotters)
cmake_policy(SET CMP0028 NEW)
find_package(HDF5 COMPONENTS CXX HL REQUIRED)
find_package(Boost COMPONENTS units program_options REQUIRED)

############################################################################################
#libcloudphxx
find_package(libcloudph++)
if (NOT libcloudph++_FOUND)
  message(FATAL_ERROR "libclouphxx++ not found")
endif()

############################################################################################
# Blitz++
# with pkgconfig
#find_package(PkgConfig)
#pkg_check_modules(BLITZ blitz>=0.10)
#if (NOT BLITZ_FOUND)
#  message(FATAL_ERROR "Blitz++ library not found or Blitz++ version requirement not met (>=0.10)
#
#* To insall Blitz++, please try:
#*   Debian/Ubuntu: sudo apt-get install libblitz0-dev
#*   Fedora: sudo yum install blitz-devel
#*   Homebrew: brew install blitz
#  ")
#endif()

#message("Blitz include dirs: ${BLITZ_INCLUDE_DIRS}")

# with CMake
find_package(blitz REQUIRED)
############################################################################################

add_library(UWLCM_plotters STATIC)

target_compile_features(${exe} PUBLIC cxx_std_14)
target_link_libraries(${exe} PRIVATE ${HDF5_LIBRARIES})
target_include_directories(${exe} PUBLIC ${HDF5_INCLUDE_DIRS})

# blitz from pkgconfig
#  target_link_libraries(${exe} PRIVATE ${BLITZ_LIBRARIES})
#  target_include_directories(${exe} PUBLIC ${BLITZ_INCLUDE_DIRS})

# blitz from cmake
target_link_libraries(${exe} PRIVATE blitz)

if(TARGET Boost::units AND TARGET Boost::program_options)
  message("linkiing boost targets")
  target_link_libraries(${exe} PRIVATE Boost::units Boost::program_options)
else()
  message("linkiing boost without targets")
  # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
  target_link_libraries(${exe} PRIVATE ${Boost_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${Boost_INCLUDE_DIRS})
endif()

#if(TARGET Boost::iostreams AND TARGET Boost::program_options AND TARGET Boost::filesystem)
#  message("linkiing boost targets")
#  target_link_libraries(${exe} PRIVATE Boost::iostreams Boost::program_options Boost::filesystem)
#else()
#  message("linkiing boost without targets")
#  # we dont link using targets, because they are not set if there is some discrepancy between cmake and boost version
#  target_link_libraries(${exe} PRIVATE ${Boost_LIBRARIES})
#  target_include_directories(${exe} PUBLIC ${Boost_INCLUDE_DIRS})
#endif()


############################################################################################
# installation

install(
  DIRECTORY
    include/UWLCM_plotters
  DESTINATION
    include
)

#install(TARGETS UWLCM_plotters EXPORT UWLCM_plotters-targets DESTINATION bin)
install(EXPORT UWLCM_plotters-targets DESTINATION share/UWLCM_plotters)

# generate and install a config file
include(CMakePackageConfigHelpers)

configure_package_config_file(
  UWLCM_plotters-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/UWLCM_plotters-config.cmake"
  INSTALL_DESTINATION
    share/UWLCM_plotters
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/UWLCM_plotters-config.cmake"
  DESTINATION
    share/UWLCM_plotters
)
############################################################################################



