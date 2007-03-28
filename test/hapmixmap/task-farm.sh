#!/bin/bash

# Execute a simple task farm, usable on a local multiprocessor.
# Takes a script file name as an argument, executes each line as
# a separate task, in a manner similar to ICHEC task farm facility.

# Number of concurrent tasks
N="2"

SCRIPT="$1"

if test -z "$SCRIPT"; then
	echo "Please provide script name"
	exit 1
fi


function push_list {
if test -z "$PID_LIST"
then
        PID_LIST="$1"
else
        PID_LIST="$PID_LIST:$1"
fi
}

function pop_list {
RET="$(echo -n $PID_LIST | cut -d: -f1)"
echo -n $RET
}

function list_length {
echo "$PID_LIST" | tr ':' '\n' | wc -l
}

cat "$SCRIPT" \
| while read LINE
do
	echo "Executing $LINE"
	( $LINE ) &
	push_list $!
	if test "$(list_length)" -ge "$N"
	then
		WAIT_PID="$(pop_list)"
		PID_LIST="$(echo -n $PID_LIST | cut --complement -d: -f1)"
		echo "Waiting for $WAIT_PID, list: $PID_LIST"
		wait $WAIT_PID
	fi
done
