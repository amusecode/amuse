
class ErrorCode(object):
    
    def __init__(self, majorcode, minorcode, description, formatstring = None):
        self.majorcode = majorcode
        self.minorcode = minorcode
        self.description = description
        self.formatstring = formatstring
    
    def __str__(self):
        return "E{0}.{1}".format(self.majorcode, self.minorcode)


class AmuseException(Exception):
    errorcode = ErrorCode(-1,-1, "unknown error")
    
    def __init__(self, *args):
        Exception.__init__(self, *args)
    
    def __str__(self):
        if self.errorcode.formatstring:
            return self.errorcode.formatstring.format(*self.args)
        return self.args[0]


class AmuseWarning(Warning):
    pass


class CoreException(AmuseException):
    majorcode = 0
    errorcode = ErrorCode(majorcode,-1, "core error")


class LegacyException(AmuseException):
    majorcode = 1
    errorcode = ErrorCode(majorcode,-1, "legacy code error")

