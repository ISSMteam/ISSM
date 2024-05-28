from netCDF4 import Dataset
from os import path


def netCDFRead(filename):
    def walktree(data):
        keys = list(data.groups.keys())
        yield keys
        for key in keys:
            for children in walktree(data.groups[str(key)]):
                yield children

    if path.exists(filename):
        print(('Opening {} for reading '.format(filename)))
        NCData = Dataset(filename, 'r')
        class_dict = {}

        for children in walktree(NCData):
            for child in children:
                class_dict[str(child)] = str(getattr(NCData.groups[str(child)], 'classtype') + '()')

        print(class_dict)
