#!/bin/bash

Neck_nr=5

python STL_to_fe.py $Neck_nr
python find_boundaries.py $Neck_nr
python create_code_fe.py $Neck_nr