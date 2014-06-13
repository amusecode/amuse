#!/bin/sh
touch $1
xterm -geometry 165x15 -title "$1" -e tail -f $1
