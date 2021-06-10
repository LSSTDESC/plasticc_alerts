# Generating alert packets with PLAsTiCC-simulated transients

We want to add on PLAsTiCC simulations on to the existing LSST Alert API. We'll do so by first generating a few fake diaObjects and then modifying the fluxes, before writing them out into an alert stream.

## In order to be able to run this on NERSC, there are a few commands you will need.

* Make sure you have the correct DESC kernels loaded:

source /global/common/software/lsst/common/miniconda/kernels/setup.sh


* Then pull the desc-stack-weekly to get the ap_pipe and ap_association moduels

python /global/common/software/lsst/common/miniconda/start-kernel-cli.py desc-stack-weekly


## Incorporating information from catalog simulations into alerts.

We will use `lsst/alert_packet` to obtain the current schema, use an example of the alert information in the schema, modify the fields we want to modify according to our catalog simulations and serialize the information to avro format.

To do this: we need to install the `alert_packet` code following the instructions in the README (pip install or EUPS serup).  The `validateAvroRoundTrip` will probably not work.

The steps to read in the JSON data, covert it to avro and to deserialize it are shown in `Examples/demo_serialize_alerts.ipynb` 
