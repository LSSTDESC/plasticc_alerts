# Generating alert packets with PLAsTiCC-simulated transients

We want to add on PLAsTiCC simulations on to the existing LSST Alert API. We'll do so by first generating a few fake diaObjects and then modifying the fluxes, before writing them out into an alert stream.

## In order to be able to run this on NERSC, there are a few commands you will need.

* Make sure you have the correct DESC kernels loaded:

source /global/common/software/lsst/common/miniconda/kernels/setup.sh


* Then pull the desc-stack-weekly to get the ap_pipe and ap_association moduels

python /global/common/software/lsst/common/miniconda/start-kernel-cli.py desc-stack-weekly
