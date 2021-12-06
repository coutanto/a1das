#
# definition of a few exceptions
#
# TODO: mettre à jour les infos de chaque classe

__doc__="This module contains error messages, not for end users"
class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class FileFormatError(Error):
    """
    Exception raised for errors in the file format

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class DataTypeError(Error):
    """
    Exception raised for errors in theData type

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class ReadDataError(Error):
    """
    Exception raised for errors in theData type

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class WrongValueError(Error):
    """
    Exception raised for errors in theData type

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class VirtualServerError(Error):
    """
    Exception raised for errors in theData type

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message

class ZMQSocketError(Error):
    """
    Exception raised for errors in theData type

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message