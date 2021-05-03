#!/bin/bash

if [ ! -x "$(which aoc 2>&1)" ]; then
  echo aoc not found
  exit 1
fi

mkdir -p /var/scratch/$USER/FPGA
rm -f /tmp/stop

INTERVAL=300

case `hostname` in
  fs5)		NR_PARALLEL=4
		export PARALLEL=11
		;;

  node50[1-4])	NR_PARALLEL=3
		export PARALLEL=7
		;;

  node505)	NR_PARALLEL=10
		export PARALLEL=9
		;;

  dop421)	NR_PARALLEL=5
		export PARALLEL=10
		;;

  tr01)		NR_PARALLEL=8
		export PARALLEL=6
		;;
esac

APPLICATIONS="gridder gridder_gen"

unset DISPLAY
nohup make -f Makefile -k -j$NR_PARALLEL\
  `for i in \`seq $INTERVAL $INTERVAL \\\`expr $NR_PARALLEL \* $INTERVAL - $INTERVAL\\\`\`;do echo sleep.$i;done`\
  `for s in \`seq 400 599\`;do for j in i${APPLICATIONS};do for g in \`seq 1 1\`;do echo /var/scratch/$USER/FPGA/${j}_LU_${g}_${s}/build;done;done;done` 2>&1 &
