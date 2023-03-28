#!/bin/sh
/opt/conda/bin/conda init bash
. /root/.bashrc
conda activate
cd /app
python -m unittest discover