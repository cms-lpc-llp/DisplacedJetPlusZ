#!/bin/bash

#######################
#debugging purposes
#######################
voms-proxy-info --all
ls -l

############################
#define exec and setup cmssw
############################
executable=DisplacedJetPlusZ
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc630

#########################################
#get cmssw
########################################
cmsrel CMSSW_10_2_0

###########################
#get cmssw environment
###########################
cd CMSSW_10_2_0/src/
eval `scram runtime -sh`
echo $1

###########################
#run executable
###########################
cd -
./${executable} --input_list=$1 --output_file=$2

##########################################################
#copy outputfile to /eos space -- define in submitter code
##########################################################

