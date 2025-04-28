#!/bin/bash

grid_sice=50
lambda=15
interval_number=0
Necks=(332 574)
#0: h_zero[-0.1, 0.4] -> for the region between bifork
#1: h_zero[-0.2,0.2] -> for energy minimum
#2: h_zero[-0.1,0.1]

# Loop through Neck numbers and run in parallel
for Neck_nr in ${Necks[@]}; do

    (
        # Start Surface Evolver in the background
        evolver_pipe="Logs/evolver_fifo_${Neck_nr}"  # Unique pipe for each process
        mkfifo $evolver_pipe
        nohup evolver -x -q code_Necks_fe/code_Neck${Neck_nr}.fe < $evolver_pipe > Logs/evolver_output_${Neck_nr}.log 2>&1 &

        # Write commands to the pipe
        exec 3> $evolver_pipe

        # Initialize Evolver with the H command file
        echo 'read "H.cmd"' >&3
        echo 'edge_fixing()' >&3
        echo "iterate_h_zero_hessian_seek_obj($grid_sice,$lambda,$interval_number)" >&3 

        # Wait for completion
        echo "Waiting for all .obj files for Neck $Neck_nr to be created..."
        echo "All $grid_sice .obj files created for Neck $Neck_nr."

        # Cleanup
        echo "quit" >&3
        exec 3>&-
        rm $evolver_pipe
    ) &  # Run the entire block in the background

done

# Wait for all background processes to finish
wait
echo "All Surface Evolver processes completed."





