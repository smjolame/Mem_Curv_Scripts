#!/bin/bash

# Activate the Conda environment
source /home/bq_jlamers/.local/conda/etc/profile.d/conda.sh  # Adjust the path based on your Conda installation
conda activate Membrane_Curv

grid_sice=50
Neck_nr=5
#Necks=(334 574)
temp_dir="Min_Necks_temp"

# Ensure the temporary directory exists
mkdir -p $temp_dir
mkdir -p Logs



#for Neck_nr in ${Necks[@]}; do
    echo "Processing Neck $Neck_nr..."

    # Process each row sequentially
    for ((i = 0; i < $grid_sice; i++)); do
        echo "Processing row $i for Neck $Neck_nr..."

        # Parallelize the creation of .obj files for all `j` indices in this row
        pids=()  # Array to track background processes
        for ((j = 0; j < $grid_sice; j++)); do
            {
                obj_file="${temp_dir}/Neck${Neck_nr}_${i}${j}.obj"
                start_time=$(date +%s%3N)
                echo "Starting creation of obj_file for (i=$i, j=$j)..."

                # Run a separate instance of Surface Evolver for each (i, j)
                evolver_pipe="Logs/evolver_fifo_${Neck_nr}_${i}_${j}"
                mkfifo $evolver_pipe
                nohup evolver -q -x code_Necks_fe/Neck${Neck_nr}.dmp < $evolver_pipe > Logs/evolver_output_${Neck_nr}_${j}.log 2>&1 &
                evolver_pid=$!
                # Write commands to the pipe
                exec 3> $evolver_pipe
                echo 'read "H.cmd"' >&3
                echo "create_obj_test($i, $j, $grid_sice)" >&3
                exec 3>&-
                rm $evolver_pipe

                wait $evolver_pid 

                if [[ -f "$obj_file" ]]; then
                    end_time=$(date +%s%3N)
                    echo "Successfully created obj_file for (i=$i, j=$j)"
                else
                    echo "[WARNING] Missing $obj_file for (i=$i, j=$j)" >> Logs/errors.log
                fi
            } &
            pids+=($!)  # Save the process ID of this background job
        done

        # Wait for all background jobs in this row to complete
        for pid in "${pids[@]}"; do
            wait $pid || echo "[WARNING] Process $pid failed." >> Logs/errors.log
        done
        echo "All .obj files for Neck $Neck_nr, row $i created."

        # Process all .obj files for this row with Python
        echo "All processes (successful or failed) for row $i finished. Running Python script..."
        python Integrated_distance_not_parallel.py $i $Neck_nr $grid_sice

        # Delete all .obj files for this row after processing is finished
        echo "Deleting .obj files for Neck $Neck_nr, row $i..."
        rm -f ${temp_dir}/Neck${Neck_nr}_${i}*.obj
    done

    echo "Finished processing Neck $Neck_nr."
#done
