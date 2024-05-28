from helpers import *
from QmuSetupResponses import *


def expandresponses(md, responses):
    
    fnames = fieldnames(responses)

    # maintain order attributes were added
    dresp = OrderedStruct()

    for key in fnames:
        value = getattr(responses, key)
        exec('dresp.{} = type(value)()'.format(key))
        for j in range(len(value)):
            #call setupdesign
            exec('dresp.{}=QmuSetupResponses(md, dresp.{}, value[j])'.format(key, key))

    return dresp
