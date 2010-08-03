
class ErrorCode(object):
    
    def __init__(self, majorcode, minorcode, description, formatstring = None):
        self.majorcode = majorcode
        self.minorcode = minorcode
        self.description = description
        self.formatstring = formatstring
    
    def __str__(self):
        return "E{0}.{1}".format(self.majorcode, self.minorcode)


class AmuseException(Exception):
    formatstring = "{0}"
    
    def __init__(self, *arguments):
        Exception.__init__(self)
        self.arguments = arguments
        
    
    def __str__(self):
        return self.formatstring.format(*self.arguments)


class AmuseWarning(Warning):
    pass


class CoreException(AmuseException):
    majorcode = 0
    errorcode = ErrorCode(majorcode,-1, "core error")


class LegacyException(AmuseException):
    majorcode = 1
    errorcode = ErrorCode(majorcode,-1, "legacy code error")

