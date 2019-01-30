#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh
source test_install.sh

# test espresso installation via `make install DESTDIR=/some/dir`
function test_install_DESTDIR() {
  helper_test_install_common "/home/areinauer/Documents/hiwi/espresso/cmake-build-debug/testsuite/cmake/install//usr/local"
}

# run tests
run_test_suite

