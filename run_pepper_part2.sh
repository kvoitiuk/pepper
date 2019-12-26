#!/usr/bin/env fish

#------fix image generation
#Make images

set DATA_SRC_PTH "/data/users/common/polisher/Microbial_data/r10_native_microbial_early_basecaller/"
set IMAGE_OUTPUT_PTH "/data/users/kvoitiuk/pepper_outputs/r10_native_microbial_early_basecaller/"
set NUM_THREADS 30


#Test images--------------------------------------
echo "Generating Test Images"

set BAM $DATA_SRC_PTH"test_data/Microbial.readsR10NativeG01a_to_shastaR94G305.test_Ecoli.bam"
set DRAFT $DATA_SRC_PTH"test_data/Microbial.shasta_assembly_test_Ecoli.fasta"
set TRUTH_BAM $DATA_SRC_PTH"test_data/Microbial.truth_2_shasta.test_Ecoli.bam"

python3 1_pepper_make_images.py \
      --bam $BAM \
      --draft $DRAFT \
      --truth_bam $TRUTH_BAM \
      --train_mode \
      --output_dir $IMAGE_OUTPUT_PTH"test_images" \
      --threads $NUM_THREADS

echo "Done"

#----------------------------------


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
