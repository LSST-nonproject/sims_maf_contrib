Add new stackers and metrics here. 

After placing a new class (for a stacker or metric) in this directory, edit the 
__init__.py file to add a line importing the new metrics/stackers classes from your 
file (YourNewFile.py) as follows: 

from YourNewFile import *

 ---- (NOT QUITE TRUE YET, BUT WILL BE IN THE NEW RELEASE WHICH SHOULD BE OUT NEXT WEEK)
After ensuring the directory mafContrib is in your PYTHONPATH, you will then be
able to access your new metrics/stackers as 
'mafContrib.MyNewMetric' (MyNewMetric being the name of your new metric). 
