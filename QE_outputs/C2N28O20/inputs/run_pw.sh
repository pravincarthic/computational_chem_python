#!/bin/bash
export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}
exec /root/q-e-qe-7.6/bin/pw.x -nb 8 -in /root/inputs/input.in > /root/output_pt3.out
