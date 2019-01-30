#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh

# test if Python module can be imported
function test_Python() {
  # test espresso installation via `make install DESTDIR=/some/dir`
  assert_return_code "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/pypresso" -c "import sys;sys.path.insert(0, '/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/testsuite/cmake/install//usr/local/lib/python2.7/dist-packages');import espressomd"
  # test espresso installation via `cmake -DCMAKE_INSTALL_PREFIX=/some/dir ..`
  if [ "/usr/local" = "/tmp/espresso-unit-tests" ]
  then
    assert_return_code "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/pypresso" -c "import sys;sys.path.insert(0, '/usr/local/lib/python2.7/dist-packages');import espressomd"
  fi
}

# run tests
run_test_suite

