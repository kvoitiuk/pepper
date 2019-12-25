#!/usr/bin/env fish

set SRC_PTH "/data/users/kvoitiuk/pepper_outputs/r10_native_microbial_early_basecaller/"
set NUM_THREADS 30

#Inference
python3 2_pepper_call_consensus.py \
-i /data/users/kvoitiuk/pepper_outputs/validation_images \
-m /data/users/kvoitiuk/pepper_outputs/model_out/trained_models_12112019_145832/_epoch_17_checkpoint.pkl \
-b 512 \
-w $NUM_THREADS \
-o /data/users/kvoitiuk/pepper_outputs/polished_sequences/_epoch_17_output_polished_sequence/
-g

#Stitch
python3 3_pepper_stitch.py \
-i /data/users/kvoitiuk/pepper_outputs/polished_sequences/_epoch_17_output_polished_sequence/pepper_predictions.hdf \
-o /data/users/kvoitiuk/pepper_outputs/polished_sequences/_epoch_17_output_polished_sequence/ \
-t $NUM_THREADS
