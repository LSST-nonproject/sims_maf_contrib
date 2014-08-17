sims_maf_contrib
================

Contributed code for MAF (sims_maf)

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
