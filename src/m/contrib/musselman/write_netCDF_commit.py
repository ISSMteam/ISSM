# imports
import netCDF4
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import time
from os import path, remove
from model import *
from results import *


'''
Given a model, this set of functions will perform the following:
    1. Enter each nested class of the model.
    2. View each attribute of each nested class.
    3. Compare state of attribute in the model to an empty model class.
    4. If states are identical, pass.
    5. Otherwise, create nested groups named after class structure.
    6. Create variable named after class attribute and assign value to it.
'''


def write_netCDF(model_var, model_name: str, filename: str):
    '''
    model_var = class object to be saved
    model_name = name of class instance variable but inside quotation marks: ie if md = model(), then model_name = 'md'
    filename = path and name to save file under
    '''
    # this assigns the name model_name to the class object model_var... very important
    globals()[model_name] = model_var
    
    # Create a NetCDF file to write to
    make_NetCDF(filename)
    
    # Create an instance of an empty model class to compare model_var against
    global empty_model
    empty_model = model()

    # Walk through the model_var class and compare subclass states to empty_model
    walk_through_model(model_var, model_name)

    # in order to handle some subclasses in the results class, we have to utilize this band-aid
    # there will likely be more band-aids added unless a class name library is created with all class names that might be added to a model
    try:
        # if results has meaningful data, save the name of the subclass and class instance
        NetCDF.groups['results']
        results_subclasses_bandaid(model_var)
        # otherwise, ignore
    except KeyError:
        pass
        
    NetCDF.close()
    print('Model successfully saved as NetCDF4')
    


def results_subclasses_bandaid(model_var):
    # since the results class may have nested classes within it, we need to record the name of the 
    # nested class instance variable as it appears in the model that we're trying to save
    quality_control = []
    for class_instance_name in model_var.results.__dict__.keys():
        # for each class instance in results, see which class its from and record that info in the netcdf to recreate structure later
        # check to see if there is a solutionstep class instance
        if isinstance(model_var.results.__dict__[class_instance_name],solutionstep):
            quality_control.append(1)
            write_string_to_netcdf(variable_name=str('solutionstep'), adress_of_child=str(class_instance_name), group=NetCDF.groups['results'])
        # check to see if there is a solution class instance
        if isinstance(model_var.results.__dict__[class_instance_name],solution):
            quality_control.append(1)
            write_string_to_netcdf(variable_name=str('solution'), adress_of_child=str(class_instance_name), group=NetCDF.groups['results'])
        # check to see if there is a resultsdakota class instance
        if isinstance(model_var.results.__dict__[class_instance_name],resultsdakota):
            quality_control.append(1)
            write_string_to_netcdf(variable_name=str('resultsdakota'), adress_of_child=str(class_instance_name), group=NetCDF.groups['results'])
    if len(quality_control) != len(model_var.results.__dict__.keys()):
        print('Error: The class instance within your model.results class is not currently supported by this application')
        print(type(model_var.results.__dict__[class_instance_name]))
    else:
        print('The results class was successfully stored on disk')


    
def make_NetCDF(filename: str):
    # If file already exists delete / rename it
    if path.exists(filename):
        print('File {} allready exist'.format(filename))
    
        # If so, inqure for a new name or to do delete the existing file
        newname = input('Give a new name or "delete" to replace: ')

        if newname == 'delete':
            remove(filename)
        else:
            print(('New file name is {}'.format(newname)))
            filename = newname
    else:
        # Otherwise create the file and define it globally so other functions can call it
        global NetCDF
        NetCDF = Dataset(filename, 'w', format='NETCDF4')
        NetCDF.history = 'Created ' + time.ctime(time.time())
        NetCDF.createDimension('Unlim', None)  # unlimited dimension
        NetCDF.createDimension('float', 1)     # single integer dimension
        NetCDF.createDimension('int', 1)       # single float dimension
    
    print('Successfully created ' + filename)


    
def walk_through_model(model_var, model_name: str):
    # Iterate over first layer of model_var attributes and assume this first layer is only classes
    for group in model_var.__dict__.keys():
        adress = str(model_name + '.' + str(group))
        # Recursively walk through subclasses
        walk_through_subclasses(model_var, adress, model_name)       
        

def walk_through_subclasses(model_var, adress: str, model_name: str):
    # Iterate over each subclass' attributes
    # Use try/except since it's either a class or it's not, no unknown exceptions
    try:
        # enter the subclass, see if it has nested classes and/or attributes
        # then compare attributes between models and write to netCDF if they differ
        # if subclass found, walk through it and repeat
        for child in eval(adress + '.__dict__.keys()'):
            # make a string variable so we can send thru this func again
            adress_of_child = str(adress + '.' + str(child))
            # If the attribute is unchanged, move onto the next layer
            adress_of_child_in_empty_class = 'empty_model' + adress_of_child.removeprefix(str(model_name))
            # using try/except here because sometimes a model can have class instances/attributes that are not 
            # in the framework of an empty model. If this is the case, we move to the except statement
            try:
                if type(child) == type(eval(adress_of_child_in_empty_class)):
                    walk_through_subclasses(model_var, adress_of_child, model_name)
                # If it has been modified, record it in the NetCDF file
                else:
                    create_group(model_var, adress_of_child)
                    walk_through_subclasses(model_var, adress_of_child, model_name)
            except AttributeError:
                create_group(model_var, adress_of_child)
                walk_through_subclasses(model_var, adress_of_child, model_name)
    except Exception as e: print(e)


        
def create_group(model_var, adress_of_child):
    # start by splitting the adress_of_child into its components
    levels_of_class = adress_of_child.split('.')

    # Handle the first layer of the group(s)
    group_name = levels_of_class[1]
    group = NetCDF.createGroup(str(group_name))

    # if the data is nested, create nested groups to match class structure
    if len(levels_of_class) > 3:
        for name in levels_of_class[2:-1]:
            group = group.createGroup(str(name))
    else: pass

    # Lastly, handle the variable(s)
    variable_name = levels_of_class[-1]
    create_var(variable_name, adress_of_child, group)


def create_var(variable_name, adress_of_child, group):
    # There are lots of different variable types that we need to handle from the model class
    
    # This first conditional statement will catch numpy arrays of any dimension and save them
    if isinstance(eval(adress_of_child), np.ndarray):
        write_numpy_array_to_netcdf(variable_name, adress_of_child, group)
    
    # check if it's an int
    elif isinstance(eval(adress_of_child), int):
        variable = group.createVariable(variable_name, int, ('int',))
        variable[:] = eval(adress_of_child)
    
    # or a float
    elif isinstance(eval(adress_of_child), float):
        variable = group.createVariable(variable_name, float, ('float',))
        variable[:] = eval(adress_of_child)

    # or a string
    elif isinstance(eval(adress_of_child), str):
        write_string_to_netcdf(variable_name, adress_of_child, group)
        
    # or an empty list
    elif isinstance(eval(adress_of_child), list) and len(eval(adress_of_child))==0:
        variable = group.createVariable(variable_name, int, ('int',))

    # or a list of strings -- this needs work as it can only handle a list of 1 string
    elif isinstance(eval(adress_of_child),list) and isinstance(eval(adress_of_child)[0],str):
        for string in eval(adress_of_child):
            write_string_to_netcdf(variable_name, string, group)

    # or a regular list
    elif isinstance(eval(adress_of_child), list):
        print(eval(adress_of_child))
        variable = group.createVariable(variable_name, type(eval(adress_of_child)[0]), ('Unlim',))
        variable[:] = eval(adress_of_child)

    # anything else... (will likely need to add more cases; ie dict)
    else:
        try:
            variable = group.createVariable(variable_name, type(eval(adress_of_child)), ('Unlim',))
            variable[:] = eval(adress_of_child)
        except Exception as e: 
            print(e)
            print('Datatype given: ' + str(type(eval(adress_of_child))))

    print('Successfully transferred data from ' + adress_of_child + ' to the NetCDF')
    



def write_string_to_netcdf(variable_name, adress_of_child, group):
    # netcdf and strings dont get along.. we have to do it 'custom':
    # if we hand it an adress we need to do it this way:
    try:
        the_string_to_save = eval(adress_of_child)
        length_of_the_string = len(the_string_to_save)
        numpy_datatype = 'S' + str(length_of_the_string)
        str_out = netCDF4.stringtochar(np.array([the_string_to_save], dtype=numpy_datatype))
    #otherwise we need to treat it like a string:
    except: 
        the_string_to_save = adress_of_child
        length_of_the_string = len(the_string_to_save)
        numpy_datatype = 'S' + str(length_of_the_string)
        str_out = netCDF4.stringtochar(np.array([the_string_to_save], dtype=numpy_datatype))        

    # we'll need to make a new dimension for the string if it doesn't already exist
    name_of_dimension = 'char' + str(length_of_the_string)
    try: 
        group.createDimension(name_of_dimension, length_of_the_string)
    except: pass
    # now we can make a variable in this dimension:
    string = group.createVariable(variable_name, 'S1', (name_of_dimension))
    #finally we can write the variable:
    string[:] = str_out


def write_numpy_array_to_netcdf(variable_name, adress_of_child, group):
    # to make a nested array in netCDF, we have to get the dimensions of the array,
    # create corresponding dimensions in the netCDF file, then we can make a variable
    # in the netCDF with dimensions identical to those in the original array
    
    # start by getting the data type at the lowest level in the array:
    typeis = eval(adress_of_child + '.dtype')
    
    # if the array is 1D, we don't need to do anything fancy
    # sometimes an array has just 1 element in it though, so we need to account for those cases here:
    if len(eval(adress_of_child)) == 1:
        if typeis is np.dtype('float64'):
            variable = group.createVariable(variable_name, typeis, ('float',))
            variable[:] = eval(adress_of_child)            
        elif typeis is np.dtype('int64'):
            variable = group.createVariable(variable_name, typeis, ('int',))
            variable[:] = eval(adress_of_child)            
        else:
            variable = group.createVariable(variable_name, typeis, ('Unlim',))
            variable[:] = eval(adress_of_child)
    
    # this is the 1D case:
    elif len(np.shape(eval(adress_of_child))) == 1: 
        variable = group.createVariable(variable_name, typeis, ('Unlim',))
        variable[:] = eval(adress_of_child)
    
    # But if the array is >1D, we do need to be fancy:
    else:
        # make the dimensions
        dimensions = []
        for dimension in np.shape(eval(adress_of_child)):
            dimensions.append(str('dim' + str(dimension)))
            # if the dimension already exists we can't have a duplicate
            try:
                group.createDimension(str('dim' + str(dimension)), dimension)
            except: pass # this would mean that the dimension already exists

        # create the variable:
        variable = group.createVariable(variable_name, typeis, tuple(dimensions))

        # write the variable:
        variable[:] = eval(adress_of_child)