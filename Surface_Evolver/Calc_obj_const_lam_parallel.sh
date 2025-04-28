#!/bin/bash

grid_sice=50
lambda=10
interval_number=4
Necks=(1 2 3 4)
interval_numbers=(0 1 2)
lambda_array=(10 50)
Neck_nr=1
#0: h_zero[-0.1, 0.4] -> for the region between bifork
#1: h_zero[-0.2,0.2] -> for energy minimum
#2: h_zero[-0.1,0.1]
#3: h_zero[-0.025,0.012]
#4: h_zero[-0.075,0.012]

# Loop through Neck numbers and run in parallel

for lambda in ${lambda_array[@]}; do
    #for interval_number in ${interval_numbers[@]}; do
        for Neck_nr in ${Necks[@]}; do
    
        pids=()  # Array to track background processes
        for ((j = 0; j < $grid_sice; j++)); do
            {
                obj_file="Min_Necks/Min_Necks_temp/Neck${Neck_nr}_${lambda}_${j}_${interval_number}.obj";
                start_time=$(date +%s%3N)
                echo "[$start_time] Starting creation of $Neck_nr for h_it=$j..."

                # Run a separate instance of Surface Evolver for each (i, j)
                evolver_pipe="Logs/evolver_fifo_${Neck_nr}_${j}"
                mkfifo $evolver_pipe
                nohup evolver -q -x code_Necks_fe/Neck${Neck_nr}.dmp < $evolver_pipe > Logs/evolver_output_${Neck_nr}_${j}.log 2>&1 &

                # Write commands to the pipe
                exec 3> $evolver_pipe
                echo 'read "H.cmd"' >&3
                echo "create_obj_const_lam($j,$lambda, $grid_sice, $interval_number)" >&3
                #echo "create_obj_lambda($j,$lambda, $grid_sice)" >&3
                exec 3>&-
                rm $evolver_pipe



                wait $evolver_pid 

                if [[ -f "$obj_file" ]]; then
                    end_time=$(date +%s%3N)
                    echo "Successfully created obj_file for (j=$j)"
                else
                    echo "[WARNING] Missing $obj_file for (j=$j)" >> Logs/errors.log
                fi
            } &
            pids+=($!)  # Save the process ID of this background job
        done


        # Wait for all background jobs in this row to complete
        for pid in "${pids[@]}"; do
            wait $pid || echo "[WARNING] Process $pid failed." >> Logs/errors.log
        done

        echo "All .obj files for Neck $Neck_nr created."
        done
    
    #done
    echo "All .obj files for Intervall number $interval_number created."
done
echo "All .obj files for lambda  $lambda created."






