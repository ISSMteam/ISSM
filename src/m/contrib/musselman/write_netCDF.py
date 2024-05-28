# imports
import netCDF4
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import time
import os
from model import *
from results import *
from m1qn3inversion import m1qn3inversion
from taoinversion import taoinversion
#import OrderedStruct


'''
Given a md, this set of functions will perform the following:
    1. View each attribute of each nested class.
    2. Compare state of attribute in the model to an empty model.
    3. If states are identical, pass. (except for np arrays which will always be saved)
    4. Otherwise, create nested groups named after class structure.
    5. Create variable named after class attribute and assign value to it.
'''


def write_netCDF(md, filename: str, verbose = False):
    if verbose:
        print('Python C2NetCDF4 v1.2.0')
    else: pass
    '''
    md = model class instance to be saved
    filename = path and name to save file under
    verbose = T/F = show or muted log statements. Naturally muted
    '''
    # this is a precaution so that data is not lost
    try:
        # Create a NCData file to write to
        NCData = create_NetCDF(filename, verbose)
        
        # Create an instance of an empty md class to compare md_var against
        empty_model = model()
    
        # Walk through the md class and compare subclass states to empty_model
        walk_through_model(md, empty_model, NCData, verbose)
    
        # in order to handle some subclasses in the results class, we have to utilize this band-aid
        # there will likely be more band-aids added unless a class name library is created with all class names that might be added to a md
        try:
            # if results has meaningful data, save the name of the subclass and class instance
            NCData.groups['results']
            results_subclasses_bandaid(md, NCData, verbose)
            # otherwise, ignore
        except KeyError:
            pass
            
        NCData.close()
        if verbose:
            print('Model successfully saved as NetCDF4')
        else: pass

    # just in case something unexpected happens
    except Exception as e:
        if 'NCData' in locals():
            NCData.close()
        raise e
    

def results_subclasses_bandaid(md, NCData, verbose = False):
    # since the results class may have nested classes within it, we need to record the name of the 
    # nested class instance variable as it appears in the md that we're trying to save
    quality_control = []

    # we save lists of instances to the NCData
    solutions = []
    solutionsteps = []
    resultsdakotas = []
    
    for class_instance_name in md.results.__dict__.keys():
        if verbose:
            print(class_instance_name)
        # for each class instance in results, see which class its from and record that info in the NCData to recreate structure later
        # check to see if there is a solutionstep class instance
        if isinstance(md.results.__dict__[class_instance_name],solutionstep):
            quality_control.append(1)
            solutionsteps.append(class_instance_name)

        # check to see if there is a solution class instance
        if isinstance(md.results.__dict__[class_instance_name],solution):
            quality_control.append(1)
            solutions.append(class_instance_name)

        # check to see if there is a resultsdakota class instance
        if isinstance(md.results.__dict__[class_instance_name],resultsdakota):
            quality_control.append(1)
            resultsdakotas.append(class_instance_name)

    if solutionsteps != []:
        serialize_string(variable_name=str('solutionstep'), address_of_child=solutionsteps, group=NCData.groups['results'], list=True, NCData=NCData, verbose=verbose)

    if solutions != []:
        serialize_string(variable_name=str('solution'), address_of_child=solutions, group=NCData.groups['results'], list=True, NCData=NCData, verbose=verbose)

    if resultsdakotas != []:
        serialize_string(variable_name=str('resultsdakota'), address_of_child=resultsdakotas, group=NCData.groups['results'], list=True, NCData=NCData, verbose=verbose)

    
    if len(quality_control) != len(md.results.__dict__.keys()):
        print('Error: The class instance within your md.results class is not currently supported by this application')
        print(type(md.results.__dict__[class_instance_name]))
    else:
        if verbose:
            print('The results class was successfully stored on disk')
        else: pass


def create_NetCDF(filename: str, verbose = False):
    # If file already exists delete / rename it
    if os.path.exists(filename):
        print('File {} allready exist'.format(filename))
    
        # If so, inqure for a new name or to do delete the existing file
        newname = input('Give a new name or "delete" to replace: ')

        if newname == 'delete':
            os.remove(filename)
        else:
            print(('New file name is {}'.format(newname)))
            filename = newname
    else:
        # Otherwise create the file and define it globally so other functions can call it
        NCData = Dataset(filename, 'w', format='NETCDF4')
        NCData.history = 'Created ' + time.ctime(time.time())
        NCData.createDimension('Unlim', None)  # unlimited dimension
        NCData.createDimension('float', 1)     # single integer dimension
        NCData.createDimension('int', 1)       # single float dimension
    
    if verbose:
        print('Successfully created ' + filename)

    return NCData


def walk_through_model(md, empty_model, NCData, verbose= False):
    # Iterate over first layer of md attributes and assume this first layer is only classes
    for group in md.__dict__.keys():
        address = md.__dict__[group]
        empty_address = empty_model.__dict__[group]
        # we need to record the layers of the md so we can save them to the NCData file
        layers = [group]

        # Recursively walk through subclasses
        walk_through_subclasses(address, empty_address, layers, NCData, empty_model, verbose)       


def walk_through_subclasses(address, empty_address, layers: list, NCData, empty_model, verbose = False):
    # See if we have an object with keys or a not
    try:
        address.__dict__.keys()
        is_object = True
    except: is_object = False # this is not an object with keys

    if is_object:
        # enter the subclass, see if it has nested classes and/or attributes
        # then compare attributes between mds and write to NCData if they differ
        # if subclass found, walk through it and repeat
        for child in address.__dict__.keys():
            # record the current location
            current_layer = layers.copy()
            current_layer.append(child)
            
            # navigate to child in each md
            address_of_child = address.__dict__[child]
            
            # if the current object is a results.<solution> object and has nonzero steps attr it needs special treatment
            if isinstance(address_of_child, solution) and len(address_of_child.steps) != 0:
                create_group(address_of_child, current_layer, is_struct = True, is_special_list = False,  NCData=NCData, verbose = verbose)

            # if the current object is a list of objects (currently only filters for lists/arrays of classes)
            elif isinstance(address_of_child, list) and len(address_of_child) > 0 and hasattr(address_of_child[0], '__dict__'):
                create_group(address_of_child, current_layer, is_struct = False, is_special_list = True, NCData=NCData, verbose = verbose)

            # if the variable is an array, assume it has relevant data (this is because the next line cannot evaluate "==" with an array)
            elif isinstance(address_of_child, np.ndarray):
                create_group(address_of_child, current_layer, is_struct = False, is_special_list = False,  NCData=NCData, verbose = verbose)
            
            # see if the child exists in the empty md. If not, record it in the NCData
            else:
                try: 
                    address_of_child_in_empty_class = empty_address.__dict__[child]
                    # if that line worked, we can see how the mds' attributes at this layer compare:
    
                    # if the attributes are identical we don't need to save anything
                    if address_of_child == address_of_child_in_empty_class:
                        walk_through_subclasses(address_of_child, address_of_child_in_empty_class, current_layer, NCData, empty_model, verbose)
    
                    # If it has been modified, record it in the NCData file
                    else:
                        create_group(address_of_child, current_layer, is_struct = False, is_special_list = False,  NCData=NCData, verbose = verbose)
                        walk_through_subclasses(address_of_child, address_of_child_in_empty_class, current_layer, NCData, empty_model, verbose)
    
                except KeyError: # record in NCData and continue to walk thru md
                    walk_through_subclasses(address_of_child, empty_address, current_layer, NCData, empty_model, verbose)
                    create_group(address_of_child, current_layer, is_struct = False, is_special_list = False,  NCData=NCData, verbose = verbose)
    else: pass


def create_group(address_of_child, layers, is_struct = False, is_special_list = False,  NCData=None, verbose = False):

    # Handle the first layer of the group(s)
    group_name = layers[0]
    
    # try to make a group unless the group is already made
    try:
        group = NCData.createGroup(str(group_name))
    except:
        group = NCData.groups[str(group_name)]

    # need to check if inversion or m1qn3inversion class
    if group_name == 'inversion':
        check_inversion_class(address_of_child, NCData, verbose)
    else: pass

    # if the data is nested in md, create nested groups to match class structure
    if len(layers) > 2:
        for name in layers[1:-1]:
            try:
                group = group.createGroup(str(name))
            except:
                group = NCData.groups[str(name)]
    else: pass

    # Lastly, handle the variable(s)
    if is_struct:
        parent_struct_name = layers[-1]
        serialize_nested_results_struct(parent_struct_name, address_of_child, group, NCData, verbose)

    elif is_special_list:
        list_name = layers[-1]
        serialize_array_of_objects(list_name, address_of_child, group, NCData, verbose)
    
    else:
        variable_name = layers[-1]
        serialize_var(variable_name, address_of_child, group, NCData, verbose)
            

def singleton(func):
    """
    A decorator to ensure a function is only executed once.
    """
    def wrapper(*args, **kwargs):
        if not wrapper.has_run:
            wrapper.result = func(*args, **kwargs)
            wrapper.has_run = True
        return wrapper.result
    wrapper.has_run = False
    wrapper.result = None
    return wrapper
    

@singleton
def check_inversion_class(address_of_child, NCData, verbose = False):
    # need to make sure that we have the right inversion class: inversion, m1qn3inversion, taoinversion
    if isinstance(address_of_child, m1qn3inversion):
        serialize_string(variable_name=str('inversion_class_name'), address_of_child=str('m1qn3inversion'), group=NCData.groups['inversion'], NCData=NCData, verbose = verbose)
        if verbose:
            print('Successfully saved inversion class instance ' + 'm1qn3inversion')
    elif isinstance(address_of_child, taoinversion):
        serialize_string(variable_name=str('inversion_class_name'), address_of_child=str('taoinversion'), group=NCData.groups['inversion'], NCData=NCData, verbose = verbose)
        if verbose:
            print('Successfully saved inversion class instance ' + 'taoinversion')
    else:
        serialize_string(variable_name=str('inversion_class_name'), address_of_child=str('inversion'), group=NCData.groups['inversion'], NCData=NCData, verbose = verbose)
        if verbose:
            print('Successfully saved inversion class instance ' + 'inversion')


def serialize_nested_results_struct(parent_struct_name, address_of_struct, group, NCData, verbose = False):
    '''
        This function takes a results.solution class instance and saves the solutionstep instances from <solution>.steps to the NCData. 

        To do this, we get the number of dimensions (substructs) of the parent struct (list).
        Next, we iterate through each substruct and record the data. 
        For each substruct, we create a subgroup of the main struct.
        For each variable, we create dimensions that are assigned to each subgroup uniquely.
    '''
    if verbose:
        print("Beginning transfer of nested MATLAB struct to the NCData")
    
    # make a new subgroup to contain all the others:
    group = group.createGroup(str(parent_struct_name))

    # make sure other systems can flag the nested struct type
    serialize_string('this_is_a_nested', 'struct', group, list=False, NCData=NCData, verbose = verbose)

    # other systems know the name of the parent struct because it's covered by the results/qmu functions above
    no_of_dims = len(address_of_struct)
    for substruct in range(0, no_of_dims):
        # we start by making subgroups with nice names like "1x4"
        name_of_subgroup = '1x' + str(substruct)
        subgroup = group.createGroup(str(name_of_subgroup))

        # do some housekeeping to keep track of the current layer
        current_substruct = address_of_struct[substruct]
        substruct_fields = current_substruct.__dict__.keys()

        # now we need to iterate over each variable of the nested struct and save it to this new subgroup
        for variable in substruct_fields:
            address_of_child = current_substruct.__dict__[variable]
            serialize_var(variable, address_of_child, subgroup, NCData, verbose = verbose)
    
    if verbose:
        print(f'Successfully transferred struct {parent_struct_name} to the NCData\n')




def serialize_array_of_objects(list_name, address_of_child, group, NCData, verbose):
    if verbose: 
        print(f"Serializing array of objects.")
    
    # Get the dimensions of the cell array
    if len(np.shape(address_of_child)) > 1: 
        rows, cols = np.shape(address_of_child)
    else: rows, cols = 1, np.shape(address_of_child)[0]

    # Make subgroup to represent the array
    name_of_subgroup = f"{str(rows)}x{str(cols)}_cell_array_of_objects"
    subgroup = group.createGroup(name_of_subgroup)

    # Save the name of the cell array
    serialize_string('name_of_cell_array', list_name, subgroup, NCData, verbose)

    # Save the dimensions of the cell array
    rowsID = subgroup.createVariable('rows', int, ('int',))
    colsID = subgroup.createVariable('cols', int, ('int',))
    rowsID[:] = rows
    colsID[:] = cols


    # If this is a multidimensional cell array, iterate over rows here and cols in serialize_objects
    if rows > 1:
        for row in range(rows):
            # Make a subgroup for each row
            name_of_subgroup = f"Row_{row+1}_of_{rows}"
            subgroup = group.createGroup(name_of_subgroup)
            serialize_objects(address_of_child, subgroup, NCData, cols, verbose)
    else:
        serialize_objects(address_of_child, subgroup, NCData, cols, verbose)
        
    if verbose:
        print(f"Successfully serialized array of objects: {list_name}")


def serialize_objects(address_of_child, group, NCData, cols, verbose):
    for col in range(cols):
        # Make subgroup to contain each col of array
        name_of_subgroup = f'Col_{col+1}_of_{cols}'
        subgroup = group.createGroup(name_of_subgroup)

        # index the current item
        variable = address_of_child[col]

        # Get the kind of object we're working with:
        # see if it's a solution instance
        if isinstance(variable, solution) and len(variable.steps) != 0:
            pass
            # this needs more work...
        
        # see if it's a general class -- assume ISSM classes all have __dict__
        elif hasattr(variable, '__dict__'):
            # Handle class instances
            serialize_class_instance(variable, subgroup, NCData, verbose)
        else:
            print('ERROR: Cell arrays of mixed types are not yet supported in read_NCData!')
            print('Deserialization will not be able to complete!')
            # Handle regular data structures that are already supported
            serialize_var(variable_name, variable, subgroup, NCData, verbose)


def serialize_class_instance(instance, group, NCData, verbose):
    # get parent class name:
    name = instance.__class__.__name__

    # save the name of the class
    serialize_string(variable_name='class_is_a', address_of_child=name, group=group, NCData=NCData, verbose = verbose)

    # make subgroup to contain attributes
    name_of_subgroup = 'Properties_of_' + name
    subgroup = group.createGroup(name_of_subgroup)

    # get attributes
    keys = instance.__dict__.keys()

    for name in keys:
        serialize_var(name, instance.__dict__[name], subgroup, NCData, verbose)
    


        
def serialize_var(variable_name, address_of_child, group, NCData, verbose = False):
    # There are lots of different variable types that we need to handle from the md class
    
    # This first conditional statement will catch numpy arrays of any dimension and save them
    if isinstance(address_of_child, np.ndarray):
        serialize_numpy_array(variable_name, address_of_child, group, NCData, verbose=verbose)
    
    # check if it's an int
    elif isinstance(address_of_child, int) or isinstance(address_of_child, np.integer):
        variable = group.createVariable(variable_name, int, ('int',))
        variable[:] = address_of_child
    
    # or a float
    elif isinstance(address_of_child, float) or isinstance(address_of_child, np.floating):
        variable = group.createVariable(variable_name, float, ('float',))
        variable[:] = address_of_child

    # or a string
    elif isinstance(address_of_child, str):
        serialize_string(variable_name, address_of_child, group, NCData, verbose=verbose)

    #or a bool
    elif isinstance(address_of_child, bool) or isinstance(address_of_child, np.bool_):
        # NetCDF can't handle bool types like True/False so we convert all to int 1/0 and add an attribute named units with value 'bool'
        variable = group.createVariable(variable_name, int, ('int',))
        variable[:] = int(address_of_child)
        variable.units = "bool"
        
    # or an empty list
    elif isinstance(address_of_child, list) and len(address_of_child)==0:
        variable = group.createVariable(variable_name, int, ('int',))

    # or a list of strings -- this needs work as it can only handle a list of 1 string
    elif isinstance(address_of_child,list) and isinstance(address_of_child[0],str):
        for string in address_of_child:
            serialize_string(variable_name, string, group, list=True, NCData=NCData, verbose=verbose)

    # or a regular list
    elif isinstance(address_of_child, list):
        variable = group.createVariable(variable_name, type(address_of_child[0]), ('Unlim',))
        variable[:] = address_of_child

    # anything else... (will likely need to add more cases; ie helpers.OrderedStruct)
    else:
        try:
            variable = group.createVariable(variable_name, type(address_of_child), ('Unlim',))
            variable[:] = address_of_child
            print(f'Unrecognized variable was saved {variable_name}')
        except TypeError: pass # this would mean that we have an object, so we just let this continue to feed thru the recursive function above
        except Exception as e:
            print(f'There was error with {variable_name} in {group}')
            print("The error message is:")
            print(e)
            print('Datatype given: ' + str(type(address_of_child)))
    
    if verbose:
        print(f'Successfully transferred data from {variable_name} to the NCData')
    

def serialize_string(variable_name, address_of_child, group, list=False, NCData=None, verbose = False):
    # NCData and strings dont get along.. we have to do it 'custom':
    # if we hand it an address we need to do it this way:
    if list:
        """    
        Convert a list of strings to a numpy.char_array with utf-8 encoded elements
        and size rows x cols with each row the same # of cols and save to NCData
        as char array.
        """
        try:
            strings = address_of_child
            # get dims of array to save
            rows = len(strings)
            cols = len(max(strings, key = len))
    
            # Define dimensions for the strings
            rows_name = 'rows' + str(rows)
            cols_name = 'cols' + str(cols)
            try:
                group.createDimension(rows_name, rows)
            except: pass

            try:
                group.createDimension(cols_name, cols)
            except: pass
                
            # Create a variable to store the strings
            string_var = group.createVariable(str(variable_name), 'S1', (rows_name, cols_name))
    
            # break the list into a list of lists of words with the same length as the longest word:
            # make words same sizes by adding spaces 
            modded_strings = [word + ' ' * (len(max(strings, key=len)) - len(word)) for word in strings]
            # encoded words into list of encoded lists
            new_list = [[s.encode('utf-8') for s in word] for word in modded_strings]
    
            # make numpy char array with dims rows x cols
            arr = np.chararray((rows, cols))
    
            # fill array with list of encoded lists
            for i in range(len(new_list)):
                arr[i] = new_list[i]
    
            # save array to NCData file
            string_var[:] = arr

            if verbose:
                print(f'Saved {len(modded_strings)} strings to {variable_name}')
    
        except Exception as e:
            print(f'Error: {e}')
        
    else:
        the_string_to_save = address_of_child
        length_of_the_string = len(the_string_to_save)
        numpy_datatype = 'S' + str(length_of_the_string)
        str_out = netCDF4.stringtochar(np.array([the_string_to_save], dtype=numpy_datatype))        
    
        # we'll need to make a new dimension for the string if it doesn't already exist
        name_of_dimension = 'char' + str(length_of_the_string)
        try: 
            group.createDimension(name_of_dimension, length_of_the_string)
        except: pass
        # this is another band-aid to the results sub classes...
        try:
            # now we can make a variable in this dimension:
            string = group.createVariable(variable_name, 'S1', (name_of_dimension))
            #finally we can write the variable:
            string[:] = str_out
        #except RuntimeError: pass
        except Exception as e:
            print(f'There was an error saving a string from {variable_name}')
            print(e)


def serialize_numpy_array(variable_name, address_of_child, group, NCData, verbose = False):
    # to make a nested array in NCData, we have to get the dimensions of the array,
    # create corresponding dimensions in the NCData file, then we can make a variable
    # in the NCData with dimensions identical to those in the original array
    
    # start by getting the data type at the lowest level in the array:
    typeis = address_of_child.dtype

    # catch boolean arrays here
    if typeis == bool:
        # sometimes an array has just 1 element in it, we account for those cases here:
        if len(address_of_child) == 1:
            variable = group.createVariable(variable_name, int, ('int',))
            variable[:] = int(address_of_child)
            variable.units = "bool"
        else:
            # make the dimensions
            dimensions = []
            for dimension in np.shape(address_of_child):
                dimensions.append(str('dim' + str(dimension)))
                # if the dimension already exists we can't have a duplicate
                try:
                    group.createDimension(str('dim' + str(dimension)), dimension)
                except: pass # this would mean that the dimension already exists
    
            # create the variable:
            variable = group.createVariable(variable_name, int, tuple(dimensions))
            # write the variable:
            variable[:] = address_of_child.astype(int)
            variable.units = "bool"

    # handle all other datatypes here
    else:
        # sometimes an array has just 1 element in it, we account for those cases here:
        if len(address_of_child) == 1:
            if typeis is np.dtype('float64'):
                variable = group.createVariable(variable_name, typeis, ('float',))
                variable[:] = address_of_child[0]
            elif typeis is np.dtype('int64'):
                variable = group.createVariable(variable_name, typeis, ('int',))
                variable[:] = address_of_child[0]
            else:
                print(f'Encountered single datatype from {variable_name} that was not float64 or int64, saving under unlimited dimension, may cause errors.')
                variable = group.createVariable(variable_name, typeis, ('Unlim',))
                variable[:] = address_of_child[0]
    
        # This catches all arrays/lists:
        else:
            # make the dimensions
            dimensions = []
            for dimension in np.shape(address_of_child):
                dimensions.append(str('dim' + str(dimension)))
                # if the dimension already exists we can't have a duplicate
                try:
                    group.createDimension(str('dim' + str(dimension)), dimension)
                except: pass # this would mean that the dimension already exists
    
            # create the variable:
            variable = group.createVariable(variable_name, typeis, tuple(dimensions))
    
            # write the variable:
            variable[:] = address_of_child

            