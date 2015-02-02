# sims_maf_contrib
[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/LSST-nonproject/sims_maf_contrib?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Contributed code for MAF (sims_maf).

To browse some example MAF analyses, check out the **[sims_maf_contrib Wiki!](https://github.com/LSST-nonproject/sims_maf_contrib/wiki)**

If you get stuck, [write us an issue](https://github.com/LSST-nonproject/sims_maf_contrib/issues) and we'll improve this documentation.

## Guidelines for contributors

New metric and stacker classes go in 'mafContrib'.
Driver configuration scripts to use those new metric and stacker classes go in 'examples'.

When contributing new metrics, please do the following:
* Put your name and email in a comment at the top of the file.
* Document your code. This includes a desscription of what your metric does -- if you place this information
  in triple-quotes directly after the class description, it will be available if someone types 
  'help(YourMetricName)' in a python shell. 
* Make sure someone else can tell what is being computed and why.
* Make it clear what the units of the output are.
* Make it clear what values are better or worse.
* Feel free to include a config file that exercises the code (in examples).
    You can use the displayDict to add a caption.

Please submit your code along with a configuration file that runs it so we can easily reproduce your outputs.

## Using the contributed metrics

First, move to where you would like to install the contributed metrics and clone the repo:

    git clone  git@github.com:LSST-nonproject/sims_maf_contrib.git

Make sure you have setup the lsst environment, and declare the package with eups (only need to do this once):

    cd sims_maf_contrib
    eups declare -r . -c 

Setup the package.  

    setup sims_maf_contrib

Now you can run one of the drivers from the examples directory in your workspace - which must contain the required sqlite opsim database file. For example:

    cd work
    ln -s $SIMS_MAF_DATA_DIR/ops1_1140_sqlite.db .
    runDriver.py $SIMS_MAF_CONTRIB_DIR/examples/LensedQuasarTimeDelays.py

In this example, you would set the environment variable `$SIMS_MAF_DATA_DIR` to point at the directory where you keep your opsim database files (your cadence workshop USB stick, for example). You'll want to browse the code that has already been checked in to avoid reinventing the wheel. 
