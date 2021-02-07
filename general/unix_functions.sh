#! /bin/bash
# Author: Jacob Brown
# Email j.brown19@imperial.ac.uk
# Date:   2020-04-08
# Last Modified: 2020-04-24
# Desc: useful defind unix functions


###########################################
# timer function

time_start=$SECONDS

function timer {
	duration=$(($SECONDS - $time_start))
	echo -e "\n..........................\n"
 	echo "Time elapsed: " $duration " sec"
 	echo -e "\n..........................\n"
}


###########################################
# Memory Usage

function memoryUsage {

	echo -e "\n..........................\n"
	free -m | awk 'NR==2{printf "Memory Usage: %s/%sMB (%.2f%%)\n", $3,$2,$3*100/$2 }'
	df -h | awk '$NF=="/"{printf "Disk Usage: %d/%dGB (%s)\n", $3,$2,$5}'
	top -bn1 | grep load | awk '{printf "CPU Load: %.2f\n", $(NF-2)}' 
	echo -e "\n..........................\n"

}