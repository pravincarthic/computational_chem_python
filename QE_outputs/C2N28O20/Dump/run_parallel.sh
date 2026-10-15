#!/bin/bash
export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK}
RANK=${OMPI_COMM_WORLD_LOCAL_RANK}
SIZE=${OMPI_COMM_WORLD_LOCAL_SIZE}
T0=$(date +%s)
FILES=("$@")

for (( i=$RANK; i<${#FILES[@]}; i+=$SIZE )); do
    IN_FILE=${FILES[$i]}
    OUT_FILE="${IN_FILE%.in}.out"

    if [ -f "$OUT_FILE" ] && grep -q "JOB DONE" "$OUT_FILE"; then
        echo "Rank $RANK: skipping $IN_FILE, already done"
        continue
    fi

    if (( $(date +%s) - T0 > 36000 )); then
        echo "Rank $RANK: time limit reached, stopping"
        break
    fi

    WORKDIR="work_rank_${RANK}_idx_${i}"
    mkdir -p $WORKDIR && cp $IN_FILE $WORKDIR/in.in

    echo "Rank $RANK: starting $IN_FILE on GPU $RANK"
    /root/q-e-qe-7.6/bin/pw.x -in $WORKDIR/in.in > "$OUT_FILE.tmp"

    if grep -q "JOB DONE" "$OUT_FILE.tmp"; then
        mv "$OUT_FILE.tmp" "$OUT_FILE"
        python3 parse_forces.py "$OUT_FILE" >> forces_rank_${RANK}.jsonl
    else
        echo "Rank $RANK: $IN_FILE did not finish"
        rm -f "$OUT_FILE.tmp"
    fi
    rm -rf $WORKDIR
done
