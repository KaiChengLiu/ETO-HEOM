#!/bin/bash

# Check for JOBNAME and BATHTYPE
if [ $# -ne 2 ]; then
  echo "Usage: $0 JOBNAME BATHTYPE"
  exit 1
fi

JOBNAME="$1"
BATHTYPE="$2"
FOLDER_NAME="GPU_${JOBNAME}_${BATHTYPE}"

# Set the base home path, it should be the location of ETO-HEOM program 
HOME_PATH="/home/andrew91411/ETO-HEOM"

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
cp "$HOME_PATH/script/IO-scripts/key.key-tmpl" "$FOLDER_NAME/"
cp "$HOME_PATH/gen-2d-spectrum/gen_2d_spectrum.py" "$FOLDER_NAME/"
cp "$HOME_PATH/gen-2d-spectrum/gen_2d_spectrum_example.ipynb" "$FOLDER_NAME/"
cp "$HOME_PATH/script/IO-scripts/gen_2d_data.sh" "$FOLDER_NAME/"
cp "$HOME_PATH/script/IO-scripts/clean_2d_data.sh" "$FOLDER_NAME/"
cp "$HOME_PATH/script/gpu-script/submit_jobs.sh" "$FOLDER_NAME/"
cp "$HOME_PATH/usage/README.md" "$FOLDER_NAME/"
cp "$HOME_PATH/etom/${BATHTYPE}_ETOM.py" "$FOLDER_NAME/"
cp "$HOME_PATH/etom/${BATHTYPE}_ETOM_example.ipynb" "$FOLDER_NAME/"

# Modify submit_jobs.sh to hardcode HOME_PATH and JOBNAME
sed -i "s|^HOME_PATH=.*|HOME_PATH=\"$HOME_PATH\"|" "$FOLDER_NAME/submit_jobs.sh"
sed -i "0,/JOBNAME=.*$/s/JOBNAME=.*$/JOBNAME=\"$FOLDER_NAME\"/" "$FOLDER_NAME/submit_jobs.sh"

