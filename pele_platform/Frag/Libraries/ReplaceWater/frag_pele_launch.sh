source ~/.bashrc
conda activate linkers


export PELE="/scratch/PELE-repo/"
export LD_LIBRARY_PATH="/scratch/local_deps/pele_deps/boost-1.52.0/lib/"
export SCHRODINGER="/opt/schrodinger2020-1/"
export PYTHONPATH="/home/ablanco/repos/pele_platform/"
export PATH=$PATH:"/usr/lib64/openmpi/bin/"


python -m pele_platform.main input.yaml
