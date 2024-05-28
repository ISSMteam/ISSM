import numpy as np
from writejsfield import writejsfield


def writejsfile(filename, model, keyname):
    #WRITEJSFILE - write model file to javascript database
    #
    #   Usage:
    #      writejsfile(filename, model, keyname)
    #

    nods = len(model.x)
    nel = len(model.index)
    nx = len(model.contourx1)
    print(filename)
    fid = open(filename, 'w', 0)

    fid.write('model = {};\n')
    fid.write('model["title"] = "{0}";\n'.format(model.title))
    fid.write('model["initialZoomFactor"]={0};\n'.format(model.initialZoomFactor))
    #write index:
    fid.write('<!-- model["index"]{{{-->\n')
    fid.write('model["index"]=[')
    for i in range(0, nel - 1):
        fid.write('[{0},{1},{2}], '.format(model.index[i][0], model.index[i][1], model.index[i][2]))
    fid.write('[{0},{1},{2}]];\n'.format(model.index[-1][0], model.index[-1][1], model.index[-1][2]))
    fid.write('<!--}}}-->\n')
    print('writing model coordinates')
    writejsfield(fid, 'model["x"]', model.x, nods)
    writejsfield(fid, 'model["y"]', model.y, nods)
    writejsfield(fid, 'model["z"]', model.z, nods)
    writejsfield(fid, 'model["surface"]', model.surface, nods)
    writejsfield(fid, 'model["contourx1"]', model.contourx1, nx)
    writejsfield(fid, 'model["contoury1"]', model.contoury1, nx)
    writejsfield(fid, 'model["contourz1"]', model.contourz1, nx)
    writejsfield(fid, 'model["contourx2"]', model.contourx2, nx)
    writejsfield(fid, 'model["contoury2"]', model.contoury2, nx)
    writejsfield(fid, 'model["contourz2"]', model.contourz2, nx)

    print('writing results')
    results = model.results
    fid.write('results={};\n')

    for i in range(0, len(results)):
        fid.write('result={};\n')
        writejsfield(fid, 'result["data"]', results[i].data, nods)
        fid.write('<!--{{{-->\n')
        fid.write('result["caxis"]=[{0},{1}];\n'.format(results[i].caxis[0], results[i].caxis[1]))
        fid.write('result["label"]="{0}";\n'.format(results[i].label))
        fid.write('result["shortlabel"]="{0}";\n'.format(results[i].shortlabel))
        fid.write('result["unit"]="{0}";\n'.format(results[i].unit))
        if type(results[i].data) == np.float64:
            fid.write('result["time_range"]=[{0},{1}];\n'.format(results[i].time_range[0], results[i].time_range[1]))
        fid.write('results["{0}"]=result;\n'.format(i))
        fid.write('<!--}}}-->\n')
    fid.write('model.results=results;\n')
    fid.write('models["{0}"]=model;\n'.format(keyname))

    fid.close()
