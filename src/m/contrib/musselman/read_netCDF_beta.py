# imports
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from os import path, remove
from model import *
import re


'''
Given a NetCDF4 file, this set of functions will perform the following:
    1. Enter each group of the file.
    2. For each variable in each group, update an empty model with the variable's data
'''


# make a model framework to fill that is in the scope of this file
model_copy = model()


def read_netCDF(filename):
    # check if path exists
    if path.exists(filename):
        print('Opening {} for reading'.format(filename))

        # open the given netCDF4 file
        global NCData   
        NCData = Dataset(filename, 'r')
        # remove masks from numpy arrays for easy conversion
        NCData.set_auto_mask(False)
    

    # read the contents of the groups

    '''
    this function navigates like: 

    filename.groups.keys() -> filename.groups['group1'] -> 
    filename.groups['group1'].groups.keys() -> filename.groups['group1'].groups['group1.1'] ->
    filename.groups['group1'].groups['group1.1'].groups.keys() ->
    filename.groups['group1'].groups['group1.1'].groups['group1.1.1'] etc. etc.
    '''
    # walk through each group looking for subgroups and variables
    for group in NCData.groups.keys():
        print('walking ' + str(group))
        # have to send a custom name to this function: filename.groups['group']
        name = "NCData.groups['" + str(group) + "']"
        print('name sent to walker is: ' + name)
        walk_nested_groups(name)
    
    return model_copy


def walk_nested_groups(group_location_in_file):
    # first, we enter the group by: filename.groups['group_name']
    # second we search the current level for variables: filename.groups['group_name'].variables.keys()
    # third we get nested group keys by: filename.groups['group_name'].groups.keys()
    # if the variables exist, copy the data to the model framework by calling a custom function
    # if the nested groups exist, repeat. 

    for variable in eval(group_location_in_file + '.variables.keys()'):
        print('got a variable: ' + str(variable))
        location_of_variable_in_file = group_location_in_file + ".variables['" + str(variable) + "']"
        # group_location_in_file is like filename.groups['group1'].groups['group1.1'].groups['group1.1.1']
        # Define the regex pattern to match the groups within brackets
        pattern = r"\['(.*?)'\]"
        # Use regex to find all matches and return something like 'group1.group1.1.group1.1.1 ...' where the last value is the name of the variable
        matches = re.findall(pattern, location_of_variable_in_file)
        variable_name = matches[-1]
        location_of_variable_in_model = '.'.join(matches[:-1])
        copy_variable_data_to_new_model(location_of_variable_in_file, location_of_variable_in_model, variable_name)

    for nested_group in eval(group_location_in_file + '.groups.keys()'):
        print('got a nested group: ' + nested_group)
        new_nested_group = group_location_in_file + ".groups['" + str(nested_group) + "']"
        print('the location of this nested group in the file is: ' + new_nested_group)
        walk_nested_groups(new_nested_group)



def copy_variable_data_to_new_model(location_of_variable_in_file, location_of_variable_in_model, variable_name):
    # this should be as simple as navigating to the location_of_variable_in_model and setting it equal to the location_of_variable_in_file
    print('adress in file: ' + location_of_variable_in_file)
    print('adress in model: ' + location_of_variable_in_model)
    print('the value of the variable is: ')
    print(eval(location_of_variable_in_file + '[:]'))
    print('the name of the varialbe is: ' + variable_name)
    print('the type of the variable is: ' + str(type(eval(location_of_variable_in_file + '[:]'))))
    setattr(eval('model_copy.' + location_of_variable_in_model), variable_name, eval(location_of_variable_in_file + '[:]'))
    print('successfully saved var to model')