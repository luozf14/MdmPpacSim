#!/bin/bash

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
rm *.root

# Start time for the script
start_time=$(date +%s.%N)

for ((process_num=0; process_num<num_processes; process_num++))
do
  (
    # Start time for the MdmPpacSim Stage 1 process
    stage1_start_time=$(date +%s.%N)
    
    ./MdmPpacSim ../config/stage1.json "$process_num" > "log_stage1_$process_num.txt"
    
    # End time for the MdmPpacSim Stage 1 process
    stage1_end_time=$(date +%s.%N)
    
    stage1_duration=$(awk "BEGIN { mins=int(($stage1_end_time - $stage1_start_time) / 60); secs=($stage1_end_time - $stage1_start_time) % 60; printf \"%d min %d sec\", mins, secs }")
    
    echo "Time taken for MdmPpacSim Stage 1 process $process_num: $stage1_duration" >> "log_stage1_$process_num.txt"
    
    # Start time for the MdmPpacSim Stage 2 process
    stage2_start_time=$(date +%s.%N)
    
    ./MdmPpacSim ../config/stage2.json "$process_num" > "log_stage2_$process_num.txt"
    
    # End time for the MdmPpacSim Stage 2 process
    stage2_end_time=$(date +%s.%N)
    
    stage2_duration=$(awk "BEGIN { mins=int(($stage2_end_time - $stage2_start_time) / 60); secs=($stage2_end_time - $stage2_start_time) % 60; printf \"%d min %d sec\", mins, secs }")
    
    echo "Time taken for MdmPpacSim Stage 2 process $process_num: $stage2_duration" >> "log_stage2_$process_num.txt"    
  ) &
done

# Wait for all background processes to finish
wait

# Start time for the hadd command
hadd_start_time=$(date +%s.%N)

# Merge the Stage 1 root files
hadd Stage1.root Stage1_*.root
rm Stage1_*.root

# Merge the Stage 2 root files
hadd Stage2.root Stage2_*.root
rm Stage2_*.root

# End time for the hadd command
hadd_end_time=$(date +%s.%N)

hadd_duration=$(awk "BEGIN { mins=int(($hadd_end_time - $hadd_start_time) / 60); secs=($hadd_end_time - $hadd_start_time) % 60; printf \"%d min %d sec\", mins, secs }")

# Total time taken for the script
end_time=$(date +%s.%N)
total_duration=$(awk "BEGIN { mins=int(($end_time - $start_time) / 60); secs=($end_time - $start_time) % 60; printf \"%d min %d sec\", mins, secs }")

echo "Time taken for MdmPpacSim Stage 1 and Stage 2 processes: $total_duration"
echo "Time taken for hadd command: $hadd_duration"
