#!/bin/bash
_test_dir=$1
#cd $_test_dir/..

echo '---------------------- create fake fields --------------------------------'
if [ "$(ls -A fake_fields)" ]; then
  echo 'found folder fake_fields'
else
  echo 'creating folder fake_fields'
  mkdir fake_fields
fi
if [ -f fake_fields/EXFvwind.0000001224.data ]; then
  echo 'found already existing fake fields'
else
  echo 'creating fake fields'
  #cd fake_fields
  python3 $_test_dir/create_fake_fields.py
  returned_val=$?
  if [[ $returned_val -ne 0 ]] ; then
    echo "running python create_fake_fields.py failed"
    exit 1
  fi
  echo 'done creating fake fields'
fi
echo '---------------------- create fake grids --------------------------------'
if [ "$(ls -A input)" ]; then
  echo 'found folder input'
else
  echo 'creating folder input'
  mkdir input
fi
if [ -f input/GridforLTRANS-Zlevels.nc ]; then
  echo 'found already existing fake grid'
else
  echo 'creating fake grid'
  python3 $_test_dir/create_fake_grids.py
  returned_val=$?
  if [[ $returned_val -ne 0 ]] ; then
    echo "running python create_fake_grids.py failed"
    exit 1
  fi
  echo 'done creating fake grid'
fi
exit 0
