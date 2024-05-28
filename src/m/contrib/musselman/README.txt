The write_netCDF and read_netCDF modules provide a convenient way to save and restore the state of a model class instance 
in binary format via NetCDF4. This allows users to store the class state on disk and retrieve it later, facilitating seamless 
transitions between Python and MATLAB environments.

To save a model, call either write_netCDF.py or write_netCDF.m depending on whether your class is in matlab or python. 
To read a saved model, call either read_netCDF.py or read_netCDF.m depending on what language you prefer to use the model in.

Usage Instructions:

    Python:
        - Saving a model: 
            from write_netCDF import write_netCDF

            md = bamg(model(), foo.csv, .01)

            write_netCDF(md, 'md', 'adress_to_save/../filename.nc')

        - Reading a model:
            from read_netCDF import read_netCDF

            md = read_netCDF('adress_to_file/../filename.nc')

    MATLAB:
        - Saving a model: 
            import write_netCDF

            write_netCDF(md, adress_to_save/../filename.nc)

        - Reading a model:
            import read_netCDF

            md = read_netCDF(adress_to_file/../filename.nc)



Dependencies:
    Python: 
        - NumPy 
        - NetCDF4.Dataset
        - The model() class

    MATLAB: 
        - The model() class


Additional Information:

There are currently datatypes that both write_netCDF and read_netCDF modules may not be able to handle. These datatypes might 
include dictionaries with multiple value datatypes, lists with multiple datatypes, lists of dicts ect. 

To add functionality for these additional cases, one must simply create a function to handle the case and call it using a 
conditional case within the create_var() function. To read the data from the NetCDF4 file, add the case to the 
copy_variable_data_to_new_model() function so that the data can be added to a new model() instance.
