----------
VERA tests
----------

To run the VERA tests, run

::

   python run.py

This will create the mesh cases (which can be changed in run.py), save the input files, and run the problems.

To fit the number of processors, adjust the num_cases setting in run.py and the number_of_threads setting in template_1b/1e.xml. The total number of processors used will be the number of cases times the number of threads. 
