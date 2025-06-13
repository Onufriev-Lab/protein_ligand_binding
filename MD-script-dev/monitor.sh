#!/bin/bash

# Define the PID to monitor
target_pid=1356114

# Function to check if the process is still running
is_process_running() {
    ps -p $target_pid > /dev/null 2>&1
}

# Loop to check the process status every hour
while is_process_running; do
    echo "$(date): Process $target_pid is still running. Checking again in 1 hour."
    sleep 1h  # Sleep for 1 hour
done

# Once the process has finished, start run.sh with nohup
echo "$(date): Process $target_pid has finished. Starting run.sh."
nohup ./run.sh > run.log 2>&1 &

