#!/bin/bash

######################################################################
# Docker-based OpenSwath+pyP+TRIC analysis of SEC-SWATH-MS data ######
# Moritz Heusel 2019-09-06 ###########################################
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

# 01) Now prepare data analysis folder(s) and input data:
	# cd >myworkdir<; ca. 250 gb should be free
git clone https://github.com/heuselm/SECSWATH_PeptideCentricAnalysis.git
cd SECSWATH_PeptideCentricAnalysis

# 02) Prepare Spectral/Peptide Query parameter Library:
	# option 1: project-specific library:
	# --> copy spectrast2tsv.tsv prepared as described in Schubert et al. 2015 to /data_library/
	# option2: Use combined human library (Rosenberger 2014):
	# --> download link:
	# https://db.systemsbiology.net/sbeams/cgi/downloadFile.cgi?name=phl004_canonical_s64_osw.csv;format=tsv;tmp_file=8becf7ae782dd305c0eade59f282bcd1;raw_download=1
	# Move to data_library and rename
	# mv phl004_canonical_s64_osw.csv data_library/spectrast2tsv.tsv
# 03) Prepare input data
	mkdir data_dia
	mkdir data/data_dia/unfractionated_secinput
	# Convert .wiff to mzXML with peak picking / centroiding MS levels 1&2
	# --> copy unfractionated sample SWATH64vw .mzXML to /data_dia/unfractionated_secinput/
	# --> copy SEC fraction sample SWATH64vw .mzXMLs to /data_dia/

## Spawn OpenSwath container with working directory mapped to /data/

# remove containers that may exist in these names
# watch out, modifications made inside these containers will be lost!
docker stop openswath
docker rm openswath
# spawn container
docker run -u 0 -dit --name openswath -v $PWD/:/data openswath/openswath:0.1.2
# say hello!
docker exec openswath echo hi there, openswath container is happy and alive


######################################################################
# STEP 2: OPENSWATH analysis: ########################################
# Query library peptides in SWATH data by OpenSwathWorkflow ##########
######################################################################
# Note: the following commands have to be run from within ############
# (attached to) the openswath container ##############################
# docker exec works for OpenswathWorkflow but not PyProphet ##########
######################################################################
docker attach openswath
# press enter once to enter the bash commandline in the container

# Set up result folders
mkdir /data/results
mkdir /data/results/openswath
mkdir /data/results/openswath/unfractionated_secinput

# Convert líbrary to .pqp
TargetedFileConverter -in /data/data_library/spectrast2tsv.tsv \
-out /data/data_library/spectrast2tsv.pqp
# Generate decoys
OpenSwathDecoyGenerator -in /data/data_library/spectrast2tsv.pqp \
-out /data/data_library/spectrast2tsv_td.pqp

## OpenSwathWorkflow
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
-threads 6 \
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
-threads 6 \
-use_ms1_traces \
-Scoring:Scores:use_ms1_mi \
-Scoring:Scores:use_mi_score ; done


######################################################################
# STEP3: PYPROPHET: ##################################################
# Score and filter OpenSwath results #################################
######################################################################
# Note: This is the "new" pipeline training a global model ###########
# that is then applied to each individual run (stabilized scoring) ###
######################################################################
# create pyprophet result folders
mkdir /data/results/pyprophet
mkdir /data/results/pyprophet/unfractionated_secinput

# Train Global Model: pyProphet analysis of unfractionated sample 
#####################################################################
pyprophet score --threads 6 --in=/data/results/openswath/unfractionated_secinput/unfractionated_secinput.osw \
--out=/data/results/pyprophet/unfractionated_secinput/model.osw --level=ms1ms2

# Apply global model to score peak groups in all runs evenly
#####################################################################
# Merge
pyprophet merge --out=/data/results/pyprophet/allruns.osw \
--subsample_ratio=1 /data/results/openswath/*.osw

# Score using the unfractionated model
pyprophet score --threads 6 --in=/data/results/pyprophet/allruns.osw \
--apply_weights=/data/results/pyprophet/unfractionated_secinput/model.osw --level=ms1ms2

# Export results before pyProphet Qvalue/FDR estimation
pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns.tsv --format=legacy_merged \
--no-ipf

pyprophet export --in=/data/results/pyprophet/allruns.osw \
--out=/data/results/pyprophet/allruns_matrix.tsv --format=matrix \
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
# create TRIC result folder
mkdir /data/results/TRIC

feature_alignment.py \
--in /data/results/pyprophet/allruns.tsv \
--out /data/results/TRIC/feature_alignment.tsv \
--out_matrix /data/results/TRIC/feature_alignment_matrix.tsv \
--method LocalMST \
--realign_method lowess \
--max_rt_diff 60 \
--mst:useRTCorrection True \
--mst:Stdev_multiplier 3.0 \
--target_fdr 0.05

# exit the docker container
exit

# clean up
docker stop openswath
docker rm openswath

# Congrats, you ran the OpenSwathPipeline!
# Check output in results/ folder
# Input for CCprofiler: /data/results/TRIC/feature_alignment.tsv
