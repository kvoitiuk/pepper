#
#!/usr/bin/env fish


#-------------------------

#Train model, validate
python3 pepper_train.py \
--train_file /data/users/kvoitiuk/pepper_outputs/train_images/ \
--test_file /data/users/kvoitiuk/pepper_outputs/test_images/ \
--batch_size 512 \
--epoch_size 100 \
--model_out /data/users/kvoitiuk/pepper_outputs/model_out/ \
--num_workers 30 \
--gpu_mode true


#Pick best model...


#Inference
python3 2_pepper_call_consensus.py \
-i /data/users/kvoitiuk/pepper_outputs/validation_images \
-m /data/users/kvoitiuk/pepper_outputs/model_out/trained_models_12112019_145832/_epoch_17_checkpoint.pkl \
-b 512 \
-w 30 \
-o /data/users/kvoitiuk/pepper_outputs/polished_sequences/_epoch_17_output_polished_sequence/
-g

#Stitch
python3 3_pepper_stitch.py \
-i /data/users/kvoitiuk/pepper_outputs/polished_sequences/_epoch_17_output_polished_sequence/pepper_predictions.hdf \
-o /data/users/kvoitiuk/pepper_outputs/polished_sequences/_epoch_17_output_polished_sequence/ \
-t 30
