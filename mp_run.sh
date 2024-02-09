#!/bin/bash

config="../config/config.json"
num_processes=3

# Function to handle termination signal
function terminate_script() {
  echo "Terminating script..."
  pkill -P $$  # Send termination signal to all child processes
  exit 1
}

# Trap termination signal
trap terminate_script SIGINT SIGTERM

rm log_*.txt
rm Stage1.root

# Start time for the script
start_time=$(date +%s.%N)

for ((process_num=0; process_num<num_processes; process_num++))
do
  (
    # Start time for the MdmPpacSim process
    mdm_start_time=$(date +%s.%N)
    
    ./MdmPpacSim "$config" "$process_num" > "log_$process_num.txt"
    
    # End time for the MdmPpacSim process
    mdm_end_time=$(date +%s.%N)
    
    mdm_duration=$(awk "BEGIN { printf \"%.2f\", ($mdm_end_time - $mdm_start_time) / 60 }")
    
    echo "Time taken for MdmPpacSim process $process_num: $mdm_duration minutes" >> "log_$process_num.txt"
  ) &
done

# Wait for all background processes to finish
wait

# Start time for the hadd command
hadd_start_time=$(date +%s.%N)

# Merge the root files using hadd
hadd -f Stage1.root Stage1_*.root

# End time for the hadd command
hadd_end_time=$(date +%s.%N)

hadd_duration=$(awk "BEGIN { printf \"%.2f\", ($hadd_end_time - $hadd_start_time) / 60 }")

# Remove individual root files
rm Stage1_*.root

# Total time taken for the script
end_time=$(date +%s.%N)
total_time=$(awk "BEGIN { printf \"%.2f\", ($end_time - $start_time) / 60 }")

echo "Time taken for MdmPpacSim processes: $total_time minutes"
echo "Time taken for hadd command: $hadd_duration minutes"
