#!/bin/bash
if [ -z ${PLOT_LIB_DIR+x} ] ; then
	SOURCE="${BASH_SOURCE[0]}"
	while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
	  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	  SOURCE="$(readlink "$SOURCE")"
	  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
	done
	export PLOT_LIB_DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
	export PATH+=:$PLOT_LIB_DIR/bin
	if [ $(uname) = "Darwin" ]; then
		export DYLD_LIBRARY_PATH+=:$PLOT_LIB_DIR/lib
	else
		export LD_LIBRARY_PATH+=:$PLOT_LIB_DIR/bin
	fi
	alias plotcd='cd $PLOT_LIB_DIR'
	alias plotmake='cd $PLOT_LIB_DIR/build; make -j4; cd -'
fi