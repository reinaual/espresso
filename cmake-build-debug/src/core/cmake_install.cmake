# Install script for directory: /home/areinauer/Documents/hiwi/espresso/src/core

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
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so.4" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so.4")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so.4"
         RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd" TYPE SHARED_LIBRARY FILES "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/EspressoConfig.so.4")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so.4" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so.4")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so.4"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so.4")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so"
         RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd" TYPE SHARED_LIBRARY FILES "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/EspressoConfig.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoConfig.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so.4" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so.4")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so.4"
         RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd" TYPE SHARED_LIBRARY FILES "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/EspressoCore.so.4")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so.4" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so.4")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so.4"
         OLD_RPATH "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core:/usr/lib/x86_64-linux-gnu/hdf5/openmpi:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/actor:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/object-in-fluid:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/immersed_boundary:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/shapes:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/constraints:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/accumulators:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/cluster_analysis:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/observables:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/io/writer:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/io/mpiio:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/io/reader:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/bonded_interactions:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/nonbonded_interactions:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/grid_based_algorithms:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/virtual_sites:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/electrostatics_magnetostatics:"
         NEW_RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so.4")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so"
         RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd" TYPE SHARED_LIBRARY FILES "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/EspressoCore.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so"
         OLD_RPATH "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core:/usr/lib/x86_64-linux-gnu/hdf5/openmpi:/usr/lib/x86_64-linux-gnu/openmpi/lib:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/actor:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/object-in-fluid:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/immersed_boundary:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/shapes:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/constraints:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/accumulators:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/cluster_analysis:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/observables:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/io/writer:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/io/mpiio:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/io/reader:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/bonded_interactions:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/nonbonded_interactions:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/grid_based_algorithms:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/virtual_sites:/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/electrostatics_magnetostatics:"
         NEW_RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCore.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so.4" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so.4")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so.4"
         RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd" TYPE SHARED_LIBRARY FILES "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/EspressoCuda.so.4")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so.4" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so.4")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so.4"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so.4")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so"
         RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd" TYPE SHARED_LIBRARY FILES "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/EspressoCuda.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so"
         OLD_RPATH ":::::::::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "/usr/local/lib/python2.7/dist-packages/espressomd")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/python2.7/dist-packages/espressomd/EspressoCuda.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/observables/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/accumulators/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/io/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/actor/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/bonded_interactions/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/grid_based_algorithms/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/immersed_boundary/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/virtual_sites/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/nonbonded_interactions/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/object-in-fluid/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/shapes/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/constraints/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/cluster_analysis/cmake_install.cmake")
  include("/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/src/core/electrostatics_magnetostatics/cmake_install.cmake")

endif()

