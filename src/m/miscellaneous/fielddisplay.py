import numpy as np

import MatlabFuncs as m


def fielddisplay(md, name, comment):
    '''
    FIELDDISPLAY - display model field

       Usage:
          fielddisplay(md, name, comment)
    '''

    #get field
    field = getattr(md, name)

    #disp corresponding line as a function of field type (offset set as 9 spaces)
    return parsedisplay("         ", name, field, comment)


def parsedisplay(offset, name, field, comment):  # {{{
    #string
    if isinstance(field, str):
        if len(field) > 30:
            string = displayunit(offset, name, "not displayed", comment)
        else:
            string = displayunit(offset, name, "'%s'" % field, comment)

    #numeric
    elif isinstance(field, (int, float)):
        string = displayunit(offset, name, str(field), comment)

    #matrix
    elif isinstance(field, np.ndarray):
        string = displayunit(offset, name, str(field.shape), comment)

    #logical
    elif isinstance(field, bool):
        if field:
            string = displayunit(offset, name, "True", comment)
        else:
            string = displayunit(offset, name, "False", comment)

    #dictionary
    elif isinstance(field, dict):
        string = dict_display(offset, name, field, comment)

    #list or tuple
    elif isinstance(field, (list, tuple)):
        string = list_display(offset, name, field, comment)

    #None
    elif field is None:
        string = displayunit(offset, name, "None", comment)

    else:
        string = displayunit(offset, name, "not displayed", comment)

    return string
    # }}}


def dict_display(offset, name, field, comment):  # {{{
    if field:
        string = displayunit(offset, name, '{dictionary}', comment) + '\n'
        offset += '   '

        for structure_field, sfield in list(field.items()):
            string += parsedisplay(offset, str(structure_field), sfield, '') + '\n'

        if string and string[-1] == '\n':
            string = string[:-1]

    else:
        string = displayunit(offset, name, 'N/A', comment)

    return string
    # }}}


def list_display(offset, name, field, comment):  # {{{
    #initialization
    if isinstance(field, list):
        sbeg = '['
        send = ']'
    elif isinstance(field, tuple):
        sbeg = '('
        send = ')'
    string = sbeg

    #go through the cell and fill string
    if len(field) < 5:
        for fieldi in field:
            if isinstance(fieldi, str):
                string += "'%s', " % fieldi
            elif isinstance(fieldi, (bool, int, float)):
                string += "%s, " % str(fieldi)
            else:
                string = sbeg
                break

    if m.strcmp(string, sbeg):
        string = "%s%dx1%s" % (sbeg, len(field), send)
    else:
        string = string[:-1] + send

    #call displayunit
    return displayunit(offset, name, string, comment)
    # }}}


def displayunit(offset, name, characterization, comment):  #{{{
    #take care of name
    if len(name) > 23:
        name = "%s..." % name[:20]

    #take care of characterization
    if characterization in ["''", '""', 'nan', np.nan, 'NaN', "[0x1]"]:
        characterization = "N/A"

    if len(characterization) > 15:
        characterization = "%s..." % characterization[:12]

    #print
    if not comment:
        string = "%s% - 23s: % - 15s" % (offset, name, characterization)
    else:
        if isinstance(comment, str):
            string = "%s% - 23s: % - 15s -- %s" % (offset, name, characterization, comment)
        elif isinstance(comment, list):
            string = "%s% - 23s: % - 15s -- %s" % (offset, name, characterization, comment[0])
            for commenti in comment:
                string += "\n%s% - 23s  % - 15s    %s" % (offset, '', '', commenti)
        else:
            raise RuntimeError("fielddisplay error message: format for comment not supported yet")

    return string
    # }}}
