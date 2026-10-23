#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate haddock3

/trinity/login/rversini/miniconda3/envs/haddock3/bin/haddock3 "$@"
