#!/bin/bash

# Stop the easy-dsp daemons in case they are running
/home/lcav/easy-dsp/easy_dsp/board_daemons/stop.sh

# Stop apache
service apache2 stop

# Modify path to local shared libraries
export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH
