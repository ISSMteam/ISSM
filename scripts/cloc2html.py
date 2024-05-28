#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-
#inspired from http://qwiki.stanford.edu/images/d/df/Latex2qwiki.txt
import re
import os

ISSM_DIR = os.getenv('ISSM_DIR')
if(not ISSM_DIR):
    raise NameError('ISSM_DIR undefined')

infile = open('temp', 'r')
outfile = open('temp.html', 'w')
file_text = infile.readlines()

#write header
outfile.write('<table width="600px" rules=none border=0 bordercolor="#000000" cellpadding="3" align="center" style="border-collapse:collapse;">\n')
style_r = 'style="text-align:right;"'
style_c = 'style="text-align:center;"'
style_l = 'style="text-align:left;"'
color = ' bgcolor=#7AA9DD '  #dark blue
color1 = ' bgcolor=#C6E2FF '  #light blue
color2 = ' bgcolor=#FFFFFF '  #white

count = 0
toggle = 0
for i in range(len(file_text)):

    #Get current lines except if first line
    if(i == 0):
        continue
    line = file_text[i]

    pattern = r"----------------"
    if (re.search(pattern, line)):
        count += 1
        continue

    if(count == 1):
        mystr = ' <tr> \n'
        column = 1
        for i in line.split():
            if(column == 1):
                mystr += '<th ' + color + style_l + '>' + i + '</th>'
                column += 1
            else:
                mystr += '<th ' + color + style_r + '>' + i + '</th>'
        mystr += '<th ' + color + style_r + '>Total</th>\n</th>\n'
    elif(count == 2):
        total = 0
        column = 1
        if(toggle):
            mystr = '<tr>\n<th ' + color1 + style_l + '>'
        else:
            mystr = '<tr>\n<th ' + color2 + style_l + '>'
        for i in line.split():
            if(not i.isdigit() or (i.isdigit and int(i) == 77) or (i.isdigit and int(i) == 90)):
                mystr += ' ' + i + ' '
            else:
                if(column == 1):
                    mystr += '</th>'
                if(column >= 2):
                    total += int(i)
                if(toggle):
                    mystr += '<td ' + color1 + style_r + '>' + i + '</td>'
                else:
                    mystr += '<td ' + color2 + style_r + '>' + i + '</td>'
                column += 1
        if(toggle):
            mystr += '<td ' + color1 + style_r + '>' + str(total) + '</td>\n</tr>\n'
        else:
            mystr += '<td ' + color2 + style_r + '>' + str(total) + '</td>\n</tr>\n'
        toggle = 1 - toggle
    elif(count == 3):
        total = 0
        column = 1
        if(toggle):
            mystr = '<tr>\n<th ' + color1 + style_l + '>'
        else:
            mystr = '<tr>\n<th ' + color2 + style_l + '>'
        for i in line.split():
            if(not i.isdigit()):
                mystr += ' ' + i + ' '
            else:
                if(column == 1):
                    mystr += '</th>'
                if(column >= 2):
                    total += int(i)
                if(toggle):
                    mystr += '<td ' + color1 + style_r + '>' + i + '</td>'
                else:
                    mystr += '<td ' + color2 + style_r + '>' + i + '</td>'
                column += 1
        if(toggle):
            mystr += '<td ' + color1 + style_r + '>' + str(total) + '</td>\n</tr>\n'
        else:
            mystr += '<td ' + color2 + style_r + '>' + str(total) + '</td>\n</tr>\n'
    else:
        continue

    outfile.write(mystr)

#write header
outfile.write(' </table>\n')

#close all files
infile.close()
outfile.close()
