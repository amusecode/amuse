from amuse.legacy.interface import create_definition

class LegacyDocStringProperty(object):            
    """
    Return a docstring generated from a legacy
    function specification
    """
    def __get__(self, instance, owner):
        if instance is None:
            if hasattr(owner, "__init__"):
                return owner.__init__.__doc__
            else:
                return self
        
        usecase = create_definition.CreateDescriptionOfALegacyFunctionDefinition()
        usecase.specification = instance.specification
        usecase.start()
        return usecase.out.string
