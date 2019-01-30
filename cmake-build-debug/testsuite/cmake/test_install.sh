#!/usr/bin/env bash

# test espresso installation
function helper_test_install_common() {
  local root=$1
  local filepaths=("${root}/bin/pypresso" \
                   "${root}/lib/python2.7/dist-packages/espressomd/EspressoCore.so" \
                  )
  if [ "TRUE" = "TRUE" ]
  then
    filepaths+=("${root}/lib/python2.7/dist-packages/espressomd/_init.so" \
                "${root}/lib/python2.7/dist-packages/espressomd/__init__.py"
               )
  fi
  for filepath in ${filepaths[@]}
  do
    assert_file_exists "${filepath}"
  done
}

