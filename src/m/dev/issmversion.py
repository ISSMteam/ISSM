from IssmConfig import IssmConfig


def issmversion():
    """
    ISSMVERSION - display ISSM version

        Usage:
            issmversion()
    """


print(' ')
print((IssmConfig('PACKAGE_NAME')[0] + ' Version ' + IssmConfig('PACKAGE_VERSION')[0]))
print(('(website: ' + IssmConfig('PACKAGE_URL')[0] + ' contact: ' + IssmConfig('PACKAGE_BUGREPORT')[0] + ')'))
print(' ')
print(('Build date: ' + IssmConfig('PACKAGE_BUILD_DATE')[0]))
print('Copyright (c) 2009-2024 California Institute of Technology')
print(' ')
print('    to get started type: issmdoc')
print(' ')
