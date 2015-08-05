Add new stackers and metrics here. 

After placing a new class (for a stacker or metric) in this directory, edit the 
__init__.py file to add a line importing the new metrics/stackers classes from your 
file (YourNewFile.py) as follows: 

from YourNewFile import *

After ensuring the directory mafContrib is in your PYTHONPATH, you will then be
able to access your new metrics/stackers as 
'mafContrib.MyNewMetric' (MyNewMetric being the name of your new metric). 
   (please update your version of MAF if this does not work). 

Please make an ipython notebook in the examples/ directory that demonstrates how your metric works and explains the output.

SOME CODING STANDARDS 
---------------------
It just makes things easier if we have some general guidelines. 

- Your metric file name should be camel case (start with lower case, later 'words' are upper case)
   example: seasonStacker
- You class names should start with upper case
   example: SeasonStacker

- Don't import modules you don't use. 
- Delete code your don't use rather than just commenting it out. 
- Keep line lengths below ~110 characters (this makes viewing on github or in editors much easier).
- Use descriptive variable names, and don't use reserved python names (e.g. not 'filter', 'id')

