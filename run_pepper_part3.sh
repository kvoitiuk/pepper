#!/usr/bin/env fish

set SRC_PTH "/data/users/kvoitiuk/pepper_outputs/r10_native_microbial_early_basecaller/"
set NUM_THREADS 30
set MODEL_PTH "model_out/trained_models_12262019_131929/_epoch_23_checkpoint.pkl"
set OUTPUT_PTH "polished_sequences/_epoch_23_output_polished_sequence/"

#Inference
#python3 2_pepper_call_consensus.py \
#-i $SRC_PTH"validation_images/" \
#-m $SRC_PTH$MODEL_PTH  \
#-b 512 \
#-w $NUM_THREADS \
#-o $SRC_PTH$OUTPUT_PTH \
#-g

#Stitch
python3 3_pepper_stitch.py \
-i $SRC_PTH$OUTPUT_PTH"pepper_predictions.hdf" \
-o $SRC_PTH$OUTPUT_PTH \
-t $NUM_THREADS
