#!/bin/bash

######################################################################
# Docker-based OpenSwath+pyP+TRIC analysis of SEC-SWATH-MS data ######
# Moritz Heusel 2019-08-19 ###########################################
######################################################################

######################################################################
# STEP 1: Preparation: ###############################################
######################################################################
# 00) Set up docker
### sudo snap install docker
### sudo usermod -a -G docker $USER
## Test docker installation
docker run hello-world
# get the containerized openswath tool(s)
docker pull openswath/openswath:0.1.2
## Spawn container with working directory mapped to /data/
# remove containers that may exist in these names
# watch out, modifications made inside these containers will be lost!
docker stop openswath
docker rm openswath
# spawn containers on host machine
docker run -u 0 -dit --name openswath -v $PWD/:/data openswath/openswath:0.1.2
# say hello!
docker exec openswath echo hi there, openswath container is happy and alive
# Great, prepare input data:
# 01) git clone analysis framework:
	# https://github.com/heuselm/SECSWATH_PeptideCentricAnalysis.git
# 02) Prepare Spectral/Peptide Query parameter Library:
	# option 1: project-specific library:
	# --> copy spectrast2tsv.tsv prepared as described in Schubert et al. 2015 to /data_library/
	# option2: Use combined human library (Rosenberger 2014):
	# --> donload link:
	# https://db.systemsbiology.net/sbeams/cgi/downloadFile.cgi?name=phl004_canonical_s64_osw.csv;format=tsv;tmp_file=8becf7ae782dd305c0eade59f282bcd1;raw_download=1
	# move to /data_library/ and rename to spectrast2tsv.tsv
# 03) Prepare input data
	# Convert .wiff to mzXML with peak picking / centroiding MS levels 1&2
	# --> copy unfractionated sample SWATH64vw .mzXML to /data_dia/unfractionated_secinput/
	# --> copy SEC fraction sample SWATH64vw .mzXMLs to /data_dia/
######################################################################


######################################################################
# STEP 2: OPENSWATH analysis: ########################################
# Query library peptides in SWATH data by OpenSwathWorkflow ##########
######################################################################
# Note: the following commands have to be run from within ############
# (attached to) the openswath container ##############################
# docker exec works for OpenswathWorkflow but not PyProphet ##########
######################################################################
docker attach openswath

# OpenSwathWorkflow
mkdir /data/results
mkdir /data/results/openswath
mkdir /data/results/openswath/unfractionated_secinput

# Convert líbrary to .pqp
TargetedFileConverter -in /data/data_library/spectrast2tsv.tsv -out /data/data_library/spectrast2tsv.pqp
# Generate decoys
OpenSwathDecoyGenerator -in /data/data_library/spectrast2tsv.pqp -out /data/data_library/spectrast2tsv_td.pqp

# Run OpenSwath on unfractionated sample(s)
for file in /data/data_dia/unfractionated_secinput/*ML; do \
bname=$(echo ${file##*/} | cut -f 1 -d '.'); \
OpenSwathWorkflow \
-in /data/data_dia/unfractionated_secinput/$bname.mzXML \
-tr /data/data_library/spectrast2tsv_td.pqp \
-tr_irt /data/data_library/irtkit.TraML \
-min_upper_edge_dist 1 \
-batchSize 1000 \
-out_osw /data/results/openswath/unfractionated_secinput/unfractionated_secinput.osw \
-Scoring:stop_report_after_feature 5 \
-rt_extraction_window 600 \
-mz_extraction_window 30 \
-ppm \
-threads 12 \
-use_ms1_traces \
-Scoring:Scores:use_ms1_mi \
-Scoring:Scores:use_mi_score ; done

# Run OpenSwath on fractionated samples
for file in /data/data_dia/*ML; do \
bname=$(echo ${file##*/} | cut -f 1 -d '.'); \
OpenSwathWorkflow \
-in /data/data_dia/$bname.mzXML \
-tr /data/data_library/spectrast2tsv_td.pqp \
-tr_irt /data/data_library/irtkit.TraML \
-min_upper_edge_dist 1 \
-batchSize 1000 \
-out_osw /data/results/openswath/$bname.osw \
-Scoring:stop_report_after_feature 5 \
-rt_extraction_window 600 \
-mz_extraction_window 30 \
-ppm \
-threads 12 \
-use_ms1_traces \
-Scoring:Scores:use_ms1_mi \
-Scoring:Scores:use_mi_score ; done

######################################################################
# STEP3: PYPROPHET: ##################################################
# Score and filter OpenSwath results #################################
######################################################################
# Note: This is the "new" pipeline training a global model ###########
# that is then applied to each individual run (stabilized srocing) ###
# Error models are learned and applied on ############################
# precursor, peptide and protein level ###############################
######################################################################
# create result folders
mkdir /data/results/pyprophet
mkdir /data/results/pyprophet/unfractionated_secinput

# Train Model: Subsample over all runs and create jumbo model
################################################
# pyprophet merge --out=/data/results/pyprophet/unfractionated_secinput/subsampled.osw \
# --subsample_ratio=0.5 /data/results/openswath/*.osw
# pyprophet score --threads 4 --in=/data/results/pyprophet/unfractionated_secinput/subsampled.osw \
# --level=ms1ms2
# # for inspection by humans:
# pyprophet export --in=/data/results/pyprophet/unfractionated_secinput/subsampled.osw \
# --out=/data/results/pyprophet/unfractionated_secinput/subsampled.tsv --format=legacy_merged \
# --no-ipf

# Train Model: pyProphet analysis of unfractionated sample 
#####################################################################
pyprophet score --threads 6 --in=/data/results/openswath/unfractionated_secinput/unfractionated_secinput.osw \
--out=/data/results/pyprophet/unfractionated_secinput/model.osw --level=ms1ms2

# Apply: Score all runs using the unfractionated model
#####################################################################
# Merge
pyprophet merge --out=/data/results/pyprophet/allruns.osw \
--subsample_ratio=1 /data/results/openswath/*.osw

# Score using the unfractionated model
pyprophet score --threads 6 --in=/data/results/pyprophet/allruns.osw \
--apply_weights=/data/results/pyprophet/unfractionated_secinput/model.osw --level=ms1ms2

# Export results before pyProphet Qvalue/FDR estimation
pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns_preQvalEstn.tsv --format=legacy_merged \
--no-ipf

# FDR estimation:
#####################################################################
# Create statistical models on peptide and protein level and over
# different contexts (Within each run or over the full dataset)
#####################################################################
# Peptide level
pyprophet \
peptide --in=/data/results/pyprophet/allruns.osw --context=run-specific \
peptide --in=/data/results/pyprophet/allruns.osw --context=experiment-wide \
peptide --in=/data/results/pyprophet/allruns.osw --context=global \
# Protein level
pyprophet \
protein --in=/data/results/pyprophet/allruns.osw --context=run-specific \
protein --in=/data/results/pyprophet/allruns.osw --context=experiment-wide \
protein --in=/data/results/pyprophet/allruns.osw --context=global \

# Export results with experiment-wide and global Qvalue estimation:
pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns_pepQvals_protQvals.tsv --format=legacy_merged \
--no-ipf

pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns_pepQvals_protQvals_matrix.tsv --format=matrix \
--no-ipf

# export scoring pdf reports
pyprophet export \
--in=/data/results/pyprophet/allruns.osw \
--format=score_plots



######################################################################
# STEP 4: TRIC Alignment: ############################################
######################################################################
# Feature alignment based on nonlinear RT alignment ##################
######################################################################
# Roest et al. https://www.nature.com/articles/nmeth.3954 ############
######################################################################
mkdir /data/results/TRIC

feature_alignment.py \
--in /data/results/pyprophet/allruns_pepQvals_protQvals.tsv \
--out /data/results/TRIC/feature_alignment.tsv \
--out_matrix /data/results/TRIC/feature_alignment_matrix.tsv \
--method LocalMST \
--realign_method lowess \
--max_rt_diff 60 \
--mst:useRTCorrection True \
--mst:Stdev_multiplier 3.0 \
--target_fdr 0.05 \
--fdr_cutoff 0.01

exit

# Congrats, you ran the OpenSwathPipeline!
# Check output in results/ folder
# Input for CCprofiler: /data/results/TRIC/feature_alignment.tsv
