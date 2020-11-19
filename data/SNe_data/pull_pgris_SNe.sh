#!/bin/bash
curl https://lsst.ncsa.illinois.edu/sim-data/SNe_data/LC_-2.0_0.2_380.0_800.0_ebvofMW_0.0_vstack.hdf5 --output LC_-2.0_0.2_380.0_800.0_ebvofMW_0.0_vstack.hdf5
curl https://lsst.ncsa.illinois.edu/sim-data/SNe_data/LC_-2.0_0.2_error_model_ebvofMW_0.0_vstack.hdf5 --output LC_-2.0_0.2_error_model_ebvofMW_0.0_vstack.hdf5
curl https://lsst.ncsa.illinois.edu/sim-data/SNe_data/LC_0.0_0.0_380.0_800.0_ebvofMW_0.0_vstack.hdf5 --output LC_0.0_0.0_380.0_800.0_ebvofMW_0.0_vstack.hdf5
curl https://lsst.ncsa.illinois.edu/sim-data/SNe_data/LC_0.0_0.0_error_model_ebvofMW_0.0_vstack.hdf5 --output LC_0.0_0.0_error_model_ebvofMW_0.0_vstack.hdf5
curl https://lsst.ncsa.illinois.edu/sim-data/SNe_data/gamma.hdf5 --output gamma.hdf5

