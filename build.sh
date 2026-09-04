#!/bin/bash

[ -d mbcclr_utils/bin ] && rm -r mbcclr_utils/bin
mkdir mbcclr_utils/bin    

case $1 in

  osx | macos)
    echo "OSX Build"
    # OSX Build (modify include/lib paths to suit your setup)
    echo "BUILDING READ ASSIGNER"
    clang++ mbcclr_utils/assign_bins.cpp -lomp -fopenmp -lpthread -o mbcclr_utils/bin/assign -I/usr/local/include -L/usr/local/lib -lz -O3
    echo "BUILD FINISHED"
    ;;

  *)
    echo "Linux Build"
    # Linux
    echo "BUILDING READ ASSIGNER"
    g++ mbcclr_utils/assign_bins.cpp -fopenmp -lpthread -o mbcclr_utils/bin/assign -lz -O3
    echo "BUILD FINISHED"
    ;;
esac
