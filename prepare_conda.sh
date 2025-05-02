# check if one argument is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: prepare_conda.sh ENV_NAME"

else

    ENV_NAME="$1"

    conda create -n $ENV_NAME python=3.10 ipykernel nb_conda_kernels numpy matplotlib pandas biopython seaborn scikit-learn
    source activate $ENV_NAME

fi
