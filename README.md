# sims_maf_contrib
[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/LSST-nonproject/sims_maf_contrib?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Contributed code for MAF (sims_maf).

To browse some example MAF analyses, check out the **[sims_maf_contrib Wiki!](https://github.com/LSST-nonproject/sims_maf_contrib/wiki)**

If you get stuck, [write us an issue](https://github.com/LSST-nonproject/sims_maf_contrib/issues) and we'll improve this documentation. 

## Guidelines for contributors

New metric and stacker classes go in 'mafContrib', and an ipython notebook demonstrating and documentating the new code goes into the relevant directory under 'science'.

When contributing new metrics, please do the following:
* Put your name and email in a comment at the top of the file.
* Document your code. This includes a short description of what your metric does -- if you place this information
  in triple-quotes directly after the class description, it will be available if someone types 
  'help(YourMetricName)' in a python shell. 
* An ipython notebook using your code allows you to include much more documentation, enough so that someone else can  tell what is being computed and why. You can include units, what the returned values are, whether larger or smaller  values are 'better'. 


## Using the contributed metrics

First, move to where you would like to install the contributed metrics and clone the repo:

    git clone  git@github.com:LSST-nonproject/sims_maf_contrib.git

Make sure you have setup the lsst environment, and declare the package with eups (only need to do this once):

    cd sims_maf_contrib
    eups declare -r . -c 

Setup the package.  

    setup sims_maf_contrib

Now you can run one of the ipython notebooks from the examples in your workspace - which must contain the required sqlite opsim database file. For example:

    cd tutorials
    ln -s [your data directory]/enigma_1189_sqlite.db .
    ipython notebook
(and then run the relevant notebooks).
