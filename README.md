# sims_maf_contrib
[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/LSST-nonproject/sims_maf_contrib?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This repository contains user contributed code for MAF (sims_maf), as well as tutorials on the use of MAF.

MAF ipython notebook tutorials can be found in the 'tutorials' directory. Start with the [Index](https://github.com/LSST-nonproject/sims_maf_contrib/blob/master/tutorials/Index.ipynb), which provides information on how to install MAF as well as an index to (some of) the tutorial notebooks.

MAF ipython notebooks demonstrating science applications can be found in the 'science' directory. 
      To browse some example MAF analyses, check out the **[sims_maf_contrib    Wiki!](https://github.com/LSST-nonproject/sims_maf_contrib/wiki)**

If you get stuck, [write us an issue](https://github.com/LSST-nonproject/sims_maf_contrib/issues) and we'll improve this documentation. 

## Guidelines for contributors

New metric and stacker python classes go in 'mafContrib', and an ipython notebook demonstrating and documentating the new code goes into the relevant directory under 'science'.

When contributing new metrics, please be sure to include an ipython notebook documenting and explaining (in words) what your metric was intended to do. Feel free to use the existing notebooks as examples. Put your name and email in your ipython notebook  and a comment near your python class. 


## Using the contributed metrics

First, move to where you would like to install the contributed metrics and clone the repo:

    git clone  git@github.com:LSST-nonproject/sims_maf_contrib.git
OR (to clone via https instead of ssh, if you do not have a github account)

    git clone  https://github.com/LSST-nonproject/sims_maf_contrib.git

Make sure you have setup the lsst environment, and declare the package with eups (you only need to do this once):

    cd sims_maf_contrib
    eups declare sims_maf_contrib -r . -t $USER

Setup the package (you have to do this every time you log into a new shell):

    setup sims_maf_contrib -t $USER -t sims

Now you can run one of the ipython notebooks from the examples in your workspace - which must contain the required sqlite opsim database file. For example:

    cd tutorials
    ln -s [your data directory]/enigma_1189_sqlite.db .
    jupyter notebook
(and then run the relevant notebooks).
