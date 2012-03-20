.. _options-label:

Options
=======

.. automodule:: amuse.support.options

    Introduction
    ------------
    The AMUSE framework provides a generic way to handle optional
    attributes for any class. When an optional attribute is requested
    the class can read the value from an configuration file or
    return a default value when nothing is found in the file. 
    The values for all options are stored in one configuration file. 
    
    Configuration file location
    ----------------------------
    The AMUSE framework will search for the configuration file in 
    following locations:

    1. First, the framework tries to find a file named **amuserc** 
       in the current working directory
    2. Second, the framework tries to find a hidden file named 
       **.amuserc** in the home directory.
    3. Last, the framework tries to find a file named **amuserc**
       in the amuse installation directory
    
    When no file is found all options will refer to the default-values.
    
    
    Configuration file format
    -------------------------
    The configuration file is formatted similar to Microsoft Windows 
    INI files. The configuration file consists of sections, led by a 
    [section] header and followed by name: value entries.
        
    For example::
        
        [channel]
        redirection=none
        debugger=xterm
    
    
    Option values
    -------------
    Values for optional attributes are determined in four
    different ways:
    
    * the attribute is set on an instance of the object::
        
        channel.debugger = "ddd"
            
    * given when creating an object of the class::
    
        channel = MessageChannnel(debugger = "ddd")
     
    * given in the ini file::
    
        [channel]
        debugger = ddd
    
    * set as default value::
    
        @option
        def debugger(self):
            return "ddd"
    
    
    Sections lookup
    ---------------
    Options can be set in different sections. To determine
    the value of an option the framework first searches
    the sections of a class and then the sections of the option.
    
    
    Use
    ---
    
    .. autoclass:: option
        :members:
        
    .. autoclass:: OptionalAttributes
        :members:
    
    
    
