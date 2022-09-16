#!/bin/bash
# takes the following arguments:
# 1. the executable to run
# 2. the input file 
# 3. the reference output file
# 4. the name of the generated output file
# 5. the test directory
# 6. the threshold value of the L2 and Linfty error-norm between ref and generated file

_LTRANS_exe=$1
_LTRANS_input=$2
_ref_output_base_name=$3
_gen_output_base_name=$4
_Adj_file_name=$5
_test_dir=$6
_threshold=$7

#if [ ! -f mesh.${_num_mpi_procs}.bas ]; then
#  echo "ERROR : bas file mesh.${_num_mpi_procs}.bas does not exist"
#  exit 1
#fi

#ln -s ${_test_dir}mesh.${_num_mpi_procs}.bas mesh.${_num_mpi_procs}.bas > /dev/null 2>&1
rm $_LTRANS_input > /dev/null 2>&1
rm output/$_gen_output_base_name > /dev/null 2>&1

echo "cp ${_test_dir}/$_LTRANS_input "
cp ${_test_dir}/$_LTRANS_input $_LTRANS_input  

#if [ ! -L ../input ]; then
#   echo "link ${_test_dir}/../input "
#   ln -s ${_test_dir}/../input ../input  > /dev/null 2>&1
#fi

if [ -f input/$_Adj_file_name ]; then
   sed -i "s/ADJele_file= .False./ADJele_file= .True./" $_LTRANS_input
fi

mkdir output  > /dev/null 2>&1
mkdir output/metadata  > /dev/null 2>&1

_LTRANS_output=${_LTRANS_input}.out
rm output/${_LTRANS_output} > /dev/null 2>&1
rm output/${_gen_output_base_name}_* > /dev/null 2>&1

cd output
echo "running tests in dir $(pwd) : $_LTRANS_exe ../$_LTRANS_input > ${_LTRANS_output}" 
$_LTRANS_exe ../$_LTRANS_input > ${_LTRANS_output} 
returned_val=$?
normal_exit_string=$( grep -q "END LTRANS" ${_LTRANS_output} ; echo $?)

# test if run was successfull 
if [[ $returned_val -ne 0 ]] || [ $normal_exit_string -ne 0 ]
then
    echo "running LTRANS $_LTRANS_input failed, returned_val=$returned_val, normal_exit_string=$normal_exit_string"
    exit 1
fi
  #for n in {1..$_num_mpi_procs}

  # compare reference and new generated output
  _gen_output_full_name=${_gen_output_base_name}
  _ref_output_full_name=${_test_dir}/${_ref_output_base_name}
#echo "paste -d', ' $_ref_output_full_name $_gen_output_full_name"
#paste -d', ' $_ref_output_full_name $_gen_output_full_name
if [ ! -f $_ref_output_full_name ]; then
    echo "file $_ref_output_full_name NOT FOUND"
    exit 1
fi
if [ ! -f $_gen_output_full_name ]; then
    echo "file $_ref_output_full_name NOT FOUND"
    exit 1
fi
  RESULTS=$(paste -d', ' $_ref_output_full_name $_gen_output_full_name | awk -v threshold=$_threshold '
                 function abs(v) {return v < 0 ? -v : v} 
                 function max(v1,v2) {return v1<v2 ? v2 : v1 } 
                 BEGIN{
                     max_err=0; 
                     sum_err2=0;
                     count=1;
                     }
           
                 ($1!=0){
                     err=abs($2-$7)+abs($3-$8)+abs($4-$9); 
                     max_err=max(max_err,err); 
                     sum_err2=sum_err2+err*err ; 
                     #printf("  --- |(%e,%e,%e)-(%e,%e,%e)| -> err=%e ;    max_err=%e, sum_err=%e \n",$2,$3,$4,$7,$8,$9,err,max_err,sum_err2)
                     }
                 ($1==0 && NR!=1){
                     if(sqrt(sum_err2)>=threshold || max_err>=threshold)
                       printf("  >>> iter: %d diff error norms are L2: %e Linfty: %e\n",count,sqrt(sum_err2),max_err)
                     max_err=0; 
                     sum_err2=0;
                     count+=1;
                     }
                 END{
                     if(sqrt(sum_err2)>=threshold || max_err>=threshold)
                        printf("  >>> error! norms are L2: %e Linfty: %e\n",sqrt(sum_err2),max_err)
                     } 
               ' ) 
               #| awk '($8>1E-08 || $10>1E-08){printf("%e %e \n",$8,$10)}' )

  if [[ $(echo -n $RESULTS | wc -c)>0 ]] ; then
    echo " " 
    echo "ERROR some L2 or Linfty error norm between ref and generated files are equals or greater than threshold $_threshold : "
    echo "$RESULTS"
    echo "floating point comparison test FAILED running LTRANS $_LTRANS_input : $_gen_output_full_name differs from $_ref_output_full_name"
    echo " " 
    exit 1
  fi


rm  *.inf > /dev/null 2>&1
if [ $_num_mpi_procs -gt 1 ] ; then
  rm mpi_debug_* > /dev/null 2>&1
fi
exit 0

