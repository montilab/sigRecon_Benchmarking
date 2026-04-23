conda create -y -n stack python=3.10
conda activate stack
pip install arc-stack
conda install -y -c conda-forge scanpy python-igraph leidenalg
qsub 01_90th_generation.sh
