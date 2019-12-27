#!/usr/bin/env fish

set SRC_PTH "/data/users/kvoitiuk/pepper_outputs/r10_native_microbial_early_basecaller/"
set NUM_THREADS 30


##Train model, validate-------------------------
set TRAIN_FILE $SRC_PTH"train_images/"
set TEST_FILE $SRC_PTH"test_images/"
set OUT_DIR $SRC_PTH"model_out/"

python3 pepper_train.py \
--train_file $TRAIN_FILE \
--test_file $TEST_FILE \
--batch_size 512 \
--epoch_size 100 \
--model_out $OUT_DIR \
--num_workers $NUM_THREADS \
--gpu_mode true


#Pick best model...
