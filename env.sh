#!/bin/bash
if [ -z ${PLOT_LIB_DIR+x} ] ; then
	if [ "$0" = "-bash" ]; then
		export PLOT_LIB_DIR=$(dirname "$(pwd)/"$(dirname "$BASH_SOURCE")"/$(basename "$BASH_SOURCE")")
	else
        if [ $? -eq 0 ]; then
            CUR_FILE=$0
    	else
            CUR_FILE=$(pwd)/$0
    	fi
		export PLOT_LIB_DIR=$(dirname "$CUR_FILE")
	fi
	export PATH+=:$PLOT_LIB_DIR/bin:$PLOT_LIB_DIR/scripts
	if [ $(uname) = "Darwin" ]; then
		export DYLD_LIBRARY_PATH+=:$PLOT_LIB_DIR/lib
	else
		export LD_LIBRARY_PATH+=:$PLOT_LIB_DIR/bin
	fi
	alias fcncmake='cd $PLOT_LIB_DIR; make; cd -'
fi