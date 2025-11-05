#!/bin/bash

#comprehensive check, except unusedFunction
#cppcheck -j 32 --include=$ISSM_DIR/config.h -DHAVE_CONFIG_H -D_HAVE_ADOLC_ -D_HAVE_DAKOTA_ -D_HAVE_JAVASCRIPT_ --enable=warning --enable=style --enable=performance --enable=portability --enable=information --enable=missingInclude $ISSM_DIR/src/c 2> CPPCHECK.err

#unused function only (slow)
$ISSM_DIR/externalpackages/cppcheck/src/cppcheck --include=$ISSM_DIR/config.h -DHAVE_CONFIG_H -D_HAVE_ADOLC_ -D_HAVE_DAKOTA_ -D_HAVE_JAVASCRIPT_ -D_HAVE_SEALEVELCHANGE_ --enable=unusedFunction $ISSM_DIR/src 2> CPPCHECK.err
