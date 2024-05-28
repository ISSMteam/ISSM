#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
import sys
import re
import os
import shutil

#get names of all directories to process
ISSM_DIR = os.getenv('ISSM_DIR')
if(not ISSM_DIR):
    raise NameError('ISSM_DIR undefined')
newclassesdir = ISSM_DIR + '/src/m/classes/'
oldclassesdir = ISSM_DIR + '/src/m/oldclasses/'

#make new directory
if(os.path.exists(oldclassesdir)):
    shutil.rmtree(oldclassesdir)
os.mkdir(oldclassesdir)

#prepare subsref and subsasgn
#{{{
subsasgntext = r'''
function obj = subsasgn(obj, index, val)
obj = builtin('subsasgn', obj, index, val)
'''
subsreftext = r'''
function obj = subsref(obj, index)
obj = builtin('subsref', obj, index)
'''
#}}}

#copy all files from new classes
files = os.listdir(newclassesdir)
for filename in files:
    newpath = newclassesdir + filename
    oldpath = oldclassesdir + filename
    if(filename == ".svn"):
        continue
    if(filename.endswith(".m")):
        shutil.copy(newpath, oldpath)
    if(filename.startswith("@")):
        shutil.copytree(newpath, oldpath)
    if(filename == "clusters"):
        files2 = os.listdir(newpath)
        for filename2 in files2:
            if(filename2 == ".svn"):
                continue
            newpath = newclassesdir + filename + '/' + filename2
            oldpath = oldclassesdir + filename2
            shutil.copy(newpath, oldpath)
    if(filename == "model"):
        shutil.copy(newpath + '/model.m', oldclassesdir + 'model.m')
    if(filename == "qmu"):
        shutil.copytree(newpath + '/@dakota_method', oldclassesdir + '/@dakota_method')
        files2 = os.listdir(newpath)
        for filename2 in files2:
            if(filename2 == ".svn"):
                continue
            if(filename2 == "@dakota_method"):
                continue
            newpath = newclassesdir + filename + '/' + filename2
            oldpath = oldclassesdir + filename2
            shutil.copy(newpath, oldpath)

files = os.listdir(oldclassesdir)
#prepare properties
#{{{
propertiesfile = open(oldclassesdir + '/properties.m', 'w')
propertiesfile.write("function out = properties(classname)\n")
#}}}
for file in files:
    if(not file.endswith(".m")):
        continue
    print("converting " + file + " from new to old Matlab class definition...")
    infile = open(oldclassesdir + file, 'r')
    classname = (re.compile(r"\.m")).sub("", file)
    dirname = oldclassesdir + '/@' + classname
    step = 0
    properties = ""

    #create directory
    if(not os.path.exists(dirname)):
        #print "Directory " + dirname + " does not exist, creating..."
        os.mkdir(dirname)

    #Process file
    file_text = infile.readlines()
    for i in range(len(file_text) - 2):
        mystr = file_text[i]

        if("properties" in mystr and step == 0):
            step = 1
            continue
        if("methods" in mystr and step == 1):
            step = 2
            continue
        if("function " in mystr and step == 2):
            step += 1

        if(step == 1):
            if("end" in mystr):
                continue
            property = mystr.lstrip()
            property = re.sub(r"%. * $", "", property)
            property = re.sub(r"\n", "", property)
            if(len(property)):
                properties = properties + 'OBJ' + property + ";\n"

        if("function " in mystr):

            #close previous file
            if(step > 3):
                outfile.close()

            #get function name
            mystr2 = (re.compile("=")).sub(" ", mystr)  #replaces equal signs by blank space
            mystr2 = (re.compile(r"\(")).sub(" (", mystr2)  #add blank spaces before and after (
            list = mystr2.split()
            for j in range(len(list)):
                word = list[j]
                if(word == '('):
                    break
            objectname = list[1]
            functionname = list[j - 1]
            objectname = re.sub(r"\[", "", objectname)
            objectname = re.sub(r"\]", "", objectname)
            if(functionname == "disp"):
                functionname = "display"
            outfile = open(dirname + '/' + functionname + '.m', 'w')

            #deal with constructor
            if(functionname == classname):

                properties2 = re.sub("OBJ", objectname + '.', properties)
                #write function declaration
                outfile.write(mystr)
                #write properties
                outfile.write(properties2)
                #write set class
                outfile.write(objectname + "=class(" + objectname + ", '" + classname + "');\n")

                #update properties list
                properties2 = properties2.split('\n')
                propertiesfile2 = open(dirname + '/properties.m', 'w')
                propertiesfile2.write("function out = properties(obj), \n")
                propertiesfile2.write('\tout = cell(' + str(len(properties2) - 1) + ', 1);\n')
                propertiesfile.write("if strcmp(classname, '" + classname + "'), \n")
                propertiesfile.write('\tout = cell(' + str(len(properties2) - 1) + ', 1);\n')
                for j in range(len(properties2) - 1):
                    property = re.sub(r"=. * $", "", properties2[j])
                    property = property.strip()
                    property = re.sub(objectname + '.', "", property)
                    propertiesfile.write('\tout{' + str(j + 1) + "} = '" + property + "';\n")
                    propertiesfile2.write('\tout{' + str(j + 1) + "} = '" + property + "';\n")
                propertiesfile.write('end\n')
                continue

        #write file
        if(step > 2):
            outfile.write(mystr)

    #close all files and delete m file
    if(step > 3):
        outfile.close()
    infile.close()
    os.remove(oldclassesdir + file)

    #Add subsref and subsasgn
    outfile = open(dirname + '/subsasgn.m', 'w')
    outfile.write(subsasgntext)
    outfile.close()
    outfile = open(dirname + '/subsref.m', 'w')
    outfile.write(subsreftext)
    outfile.close()


#close all files
propertiesfile.close()
#shutil.rmtree(newclassesdir)
