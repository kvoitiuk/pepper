#!/usr/bin/env fish

#Make images

DATA_SRC_PTH="/data/users/common/polisher/Microbial_data/r10_native_microbial_early_basecaller/"
IMAGE_OUTPUT_PTH="/data/users/kvoitiuk/pepper_outputs/r10_native_microbial_early_basecaller/"
NUM_THREADS=30


#Validation images--------------------------------------
BAM=$DATA_SRC_PTH+"validation_data/Microbial.readsR10NativeG01a_to_shastaR94G305.validate_StaphAur.bam"
DRAFT=$DATA_SRC_PTH+"validation_data/Microbial.shasta_assembly_validate_staph_aur.fasta"

python3 1_pepper_make_images.py \
      --bam $BAM \
      --draft $DRAFT \
      --output_dir $IMAGE_OUTPUT_PTH+"validation_images/" \
      --threads $NUM_THREADS


#Train images--------------------------------------
BAM=$DATA_SRC_PTH+"train_data/Microbial.readsR10NativeG01a_to_shastaR94G305.train_noEcoli_noStaphAur.bam"
DRAFT=$DATA_SRC_PTH+"train_data/Microbial.shasta_assembly_train_noEcoli_noStaphAur.fasta"
TRUTH_BAM=$DATA_SRC_PTH+"train_data/Microbial.truth_2_shasta.train_noEcoli_noStaphAur.bam"

python3 1_pepper_make_images.py \
      --bam $BAM \
      --draft $DRAFT \
      --truth_bam $TRUTH_BAM \
      --train_mode true \
      --output_dir $IMAGE_OUTPUT_PTH+"train_images/" \
      --threads $NUM_THREADS

#Test images--------------------------------------
BAM=$DATA_SRC_PTH+"test_data/Microbial.readsR10NativeG01a_to_shastaR94G305.test_Ecoli.bam"
DRAFT=$DATA_SRC_PTH+"test_data/Microbial.shasta_assembly_test_Ecoli.fasta"
TRUTH_BAM=$DATA_SRC_PTH+"test_data/Microbial.truth_2_shasta.test_Ecoli.bam"

python3 1_pepper_make_images.py \
      --bam $BAM \
      --draft $DRAFT \
      --truth_bam $TRUTH_BAM \
      --output_dir $IMAGE_OUTPUT_PTH+"test_images" \
      --threads $NUM_THREADS
