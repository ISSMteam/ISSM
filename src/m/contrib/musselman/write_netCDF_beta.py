# imports
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import time
from os import path, remove
from model import *


'''
this is like my todo list

Given a model, this set of functions will perform the following:
    1. Enter each nested class of the model.
    2. View each attribute of each nested class.
    3. Compare state of attribute in the model to an empty model class.
    4. If states are identical, pass.
    5. Else, create group named after original class.
    6. Create variable named after nested class attribute and assign value to it.
    7. 
'''

'''
need to add cases for nested arrays, dicts, etc...
To do this I need to: 
    know exactly the data types that are generating problems
'''


def write_netCDF(model_var, model_name: str, filename: str):
    '''
    model_var = class object to be saved
    model_name = name of class instance variable but inside quotation marks: ie if md = model(), then model_name = 'md'
    filename = path and name to save file under
    '''
      
    globals()[model_name] = model_var

    # a sanity check
    print('sanity check to make sure names are defined: ' + str(model_var == eval(model_name)))
    
    # Create a NetCDF file to write to
    make_NetCDF(filename)
    
    # Create an instance of an empty model class to compare model_var against
    global empty_model
    empty_model = model()

    # Walk through the model_var class and compare subclass states to empty_model
    walk_through_model(model_var, model_name)
    
    NetCDF.close()
    

    
def make_NetCDF(filename: str):
    # Check if file already exists
    if path.exists(filename):
        print('File {} allready exist'.format(filename))
    
        # If so, inqure for a new name or to do delete existing file
        newname = input('Give a new name or "delete" to replace: ')

        if newname == 'delete':
            remove(filename)
        else:
            print(('New file name is {}'.format(newname)))
            filename = newname

    # Create file and define it globally (global variables are stored in memory/global namespace)
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
        print(str(group))
        adress = str(model_name + '.' + str(group))
        print(adress)
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
            print('adress_of_child: ' + adress_of_child)
            # If the attribute is unchanged, move onto the next layer
            adress_of_child_in_empty_class = 'empty_model' + adress_of_child.removeprefix(str(model_name))
            print('adress_of_child_in_empty_class: '+ adress_of_child_in_empty_class + '\n')
            # using try/except here because sometimes a model can have class instances/attributes that are not 
            # in the framework of an empty model. If this is the case, we move to the except statement
            try:
                if type(child) == type(eval(adress_of_child_in_empty_class)):
                    print('passed a non-variable\n')
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
    print('entered create for: ' + adress_of_child + '\n')
    print('the type is: ' + str(type(eval(adress_of_child))) + '\n')
    levels_of_class = adress_of_child.split('.')
    print(levels_of_class)

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
        
    elif isinstance(eval(adress_of_child), int):
        print('caught an int!')
        variable = group.createVariable(variable_name, int, ('int',))
        variable[:] = eval(adress_of_child)
        
    elif isinstance(eval(adress_of_child), float):
        print('caught a float!')
        variable = group.createVariable(variable_name, float, ('float',))
        variable[:] = eval(adress_of_child)
        
    else:
        try:
            variable = group.createVariable(variable_name, type(eval(adress_of_child)), ('Unlim',))
            variable = eval(adress_of_child)
        except Exception as e: print(e)

    print('successfully wrote ' + adress_of_child + ' to netcdf file')
    
    
    

def write_numpy_array_to_netcdf(variable_name, adress_of_child, group):
    print('entered write_numpy_array_to_netcdf for: ' + variable_name)
    # to make a nested array in netCDF, we have to get the dimensions of the array,
    # create corresponding dimensions in the netCDF file, then we can make a variable
    # in the netCDF with dimensions identical to those in the original array
    
    # start by getting the data type at the lowest level in the array:
    typeis = eval(adress_of_child + '.dtype')
    print('the type of elements are: ' + str(typeis))
    
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

        print('the dimensions are: ' + str(dimensions))

        # create the variable:
        variable = group.createVariable(variable_name, typeis, tuple(dimensions))
        print('created variable OK')

        # write the variable:
        variable[:] = eval(adress_of_child)

        
        
        
