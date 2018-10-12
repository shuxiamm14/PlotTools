#!/bin/bash
if [ -z ${CUR_DIR+x} ] ; then
	if [ "$0" = "-bash" ]; then
		export CUR_DIR=$(dirname "$(pwd)/"$(dirname "$BASH_SOURCE")"/$(basename "$BASH_SOURCE")")
	else
        if [ $? -eq 0 ]; then
            CUR_FILE=$0
    	else
            CUR_FILE=$(pwd)/$0
    	fi
		export CUR_DIR=$(dirname "$CUR_FILE")
	fi
	export PATH+=:$CUR_DIR/bin:$CUR_DIR/scripts
	if [ $(uname) = "Darwin" ]; then
		export DYLD_LIBRARY_PATH+=:$CUR_DIR/lib
	else
		export LD_LIBRARY_PATH+=:$CUR_DIR/bin
	fi
	alias fcncmake='cd $CUR_DIR; make; cd -'
fi