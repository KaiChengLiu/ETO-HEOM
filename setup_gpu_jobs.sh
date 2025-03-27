#!/bin/bash

# Check for JOBNAME and BATHTYPE
if [ $# -ne 2 ]; then
  echo "Usage: $0 JOBNAME BATHTYPE"
  exit 1
fi

JOBNAME="$1"
BATHTYPE="$2"
FOLDER_NAME="GPU_${JOBNAME}_${BATHTYPE}"

# Set the base home path
HOME_PATH="/home/andrew91411"

# Validate BATHTYPE
if [[ "$BATHTYPE" != "debye_lorentz" && "$BATHTYPE" != "ohmic" && "$BATHTYPE" != "superohmic" ]]; then
  echo "Error: BATHTYPE must be one of: debye_lorentz, ohmic, superohmic"
  exit 1
fi

# Remove existing folder if it already exists
if [ -d "$FOLDER_NAME" ]; then
  rm -rf "$FOLDER_NAME"
fi

# Create directory structure
mkdir -p "$FOLDER_NAME/2d-input"
mkdir -p "$FOLDER_NAME/2d-output"
mkdir -p "$FOLDER_NAME/pbs-script"
mkdir -p "$FOLDER_NAME/pbs-script/pbserr"
mkdir -p "$FOLDER_NAME/pbs-script/pbslog"

# Copy files into the folder
cp "$HOME_PATH/ETO-HEOM/script/IO-scripts/key.key-tmpl" "$FOLDER_NAME/"
cp "$HOME_PATH/ETO-HEOM/gen-2d-spectrum/gen_2d_spectrum.py" "$FOLDER_NAME/"
cp "$HOME_PATH/ETO-HEOM/script/IO-scripts/gen_2d_data.sh" "$FOLDER_NAME/"
cp "$HOME_PATH/ETO-HEOM/script/IO-scripts/clean_2d_data.sh" "$FOLDER_NAME/"
cp "$HOME_PATH/ETO-HEOM/script/gpu-script/submit_jobs.sh" "$FOLDER_NAME/"
cp "$HOME_PATH/ETO-HEOM/etom/${BATHTYPE}_ETOM.py" "$FOLDER_NAME/"

# Modify submit_jobs.sh to hardcode the full folder name as JOBNAME
sed -i "0,/JOBNAME=.*$/s/JOBNAME=.*$/JOBNAME=\"$FOLDER_NAME\"/" "$FOLDER_NAME/submit_jobs.sh"

