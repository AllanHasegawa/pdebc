#!/bin/bash


for cpus in 1 2 3 4 5 6 7 8 9 10 11 12
	do
	echo -n "$cpus"
	for (( i = 0; i < 50; i++ ))
	do
		/usr/bin/time --format "%e" ./run.sh $cpus &> temp
		tempout=$(<./temp)
		echo -n ",$tempout"
		sleep 5
	done
	echo -en "\n"
done
