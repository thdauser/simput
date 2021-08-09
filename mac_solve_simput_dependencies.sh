#!/bin/bash -v

#################################################################
# Correcting the dependencies in SIMPUT libraries
# Author: E. Cucchetti (IRAP) & Thomas Dauser (ECAP), 2017-11-13
# If problems remains, contact sixte-support@lists.fau.de
#################################################################

install_prefix=$1
for tool in bin/simputfile bin/simputimg bin/simputlc bin/simputmerge bin/simputmulticell \
    bin/simputmultispec bin/simputpsd bin/simputrotate bin/simputspec bin/simputsrc bin/simputverify \
    lib/libsimput.dylib lib/libhdsp.dylib lib/libhdutils.dylib lib/libape.dylib \
    lib/libhdinit.dylib lib/libhdio.dylib lib/liblabnh.dylib lib/libposstring.dylib \
    lib/libatFunctions.dylib
do
    install_name_tool -change @rpath/libcfitsio.9.dylib $install_prefix/lib/libcfitsio.9.dylib $install_prefix/$tool
    install_name_tool -change libwcs.7.dylib $install_prefix/lib/libwcs.7.dylib $install_prefix/$tool
done
