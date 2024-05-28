from IssmConfig_python import IssmConfig_python


def IssmConfig(string):
    """
    ISSMCONFIG

        Usage:
            value = IssmConfig('string')
    """

    # Call mex module
    value = IssmConfig_python(string)
    # Return
    return value
