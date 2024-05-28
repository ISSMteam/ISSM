The write_netCDF and read_netCDF modules provide a convenient way to save and restore the state of a model class instance 
in binary format via NetCDF4. This allows users to store the class state on disk and retrieve it later, facilitating seamless 
transitions between Python and MATLAB environments.

To save a model, call either write_netCDF.py or write_netCDF.m depending on whether your class is in matlab or python. 
To read a saved model, call either read_netCDF.py or read_netCDF.m depending on what language you prefer to use the model in.
If you would like to log the names and locations of variables being stored, add the argument verbose = True (verbose = true for matlab).

Usage Instructions:

    Python:
        - Saving a model: 
            from write_netCDF import write_netCDF

            md = bamg(model(), foo.csv, .01)

            write_netCDF(md, 'adress_to_save/../filename.nc')            

        - Reading a model:
            from read_netCDF import read_netCDF

            md = read_netCDF('adress_to_file/../filename.nc')

        Verbose examples:
            write_netCDF(md, adress_to_save/../filename.nc, verbose = True)
            md = read_netCDF(adress_to_file/../filename.nc, verbose = True)

    MATLAB:
        - Saving a model:

            write_netCDF(md, adress_to_save/../filename.nc);

        - Reading a model:

            md = read_netCDF(adress_to_file/../filename.nc);

        Verbose examples:
            write_netCDF(md, adress_to_save/../filename.nc, verbose = true);
	    
          or:

	    write_netCDF(md, adress_to_save/../filename.nc, verbose);
            md = read_netCDF(adress_to_file/../filename.nc, verbose = true);

Dependencies:
    Python: 
        - NumPy 
        - NetCDF4 / NetCDF4.Dataset
        - The model() class
        - results.solution / results.solutionstep / results.resultsdakota
        - inversion.inversion / inversion.m1qn3inversion / inversion.taoinversion

    MATLAB: 
        - The model() class
        - inversion.inversion / inversion.m1qn3inversion / inversion.taoinversion


Additional Information:

There are currently datatypes that both write_netCDF and read_netCDF modules may not be able to handle. These datatypes might 
include lists with multiple datatypes (ie, ['number', 1, 'letter', a, 'color', 'blue']), lists of dicts ect. 

To add functionality for these additional cases, one must simply create a function to handle the case and call it using a 
conditional case within the create_var() function. To read the data from the NetCDF4 file, add the case to the 
copy_variable_data_to_new_model() function in read_netCDF so that the data can be added to a new model() instance.

Known issues:

Unlike Python, MATLAB doesn't utilize subclasses in its model class. This leads to a loss of certain subclass instances. 
For instance, the results.solutionstep() class poses a known issue. In MATLAB, there's no direct equivalent. The fields in 
'md.results' in MATLAB might correspond to instances of resultsdakota(), solution(), or solutionstep() in Python, but 
because those classes don't exist in MATLAB, there is no way for python to know which instance it needs. 

The current workaround, while not theoretically sound, involves searching for the class name string in MATLAB's 'results' 
field names. For instance, 'md.results.TransientSolution' is recorded as a solution() class instance. However, problems arise 
in cases like 'md.results.StressbalanceSolution', where the code notes a solution() instance, while in Python, it should be a 
solutionstep() instance.

So far, there have been no recorded problems swapping a solutionstep() instance for a solution() instance.

Potential solutions are:

    - Restructure both Python and MATLAB solve frameworks. In Python, when creating an md.results.<solutionstep()> instance, 
    embed 'solutionstep' in the class instance name.
        >> This solution is very involved, and would include the tedious modification of >5 files in total
    - Create a hash table linking solutions with their corresponding 'md.results.<class>' for reference when saving models to 
    the netCDF file. 
