#!/bin/bash
if [ "$(ls -A MITgcm_outputs)" ]; then
  echo "found MITgcm_outputs, skipping download"
else
  echo "downloading MITgcm_outputs.tar.gz"
  wget https://zenodo.org/record/3560264/files/MITgcm_outputs.tar.gz?download=1
  echo "uncompressing MITgcm_outputs.tar.gz"
  tar -zxvf 'MITgcm_outputs.tar.gz?download=1'
  rm 'MITgcm_outputs.tar.gz?download=1'
fi
