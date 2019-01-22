from amuse.support.core import late
from amuse.support import exceptions
from amuse.units import nbody_system

from amuse.rfi.async_request import DependentASyncRequest, AbstractASyncRequest, FakeASyncRequest

import inspect

class AbstractCodeMethodWrapper(object):
    
    def __init__(self, method):
        self.method = method
      
    @late
    def method_is_legacy(self):
        return hasattr(self.method, 'specification')
        
    @late
    def method_is_code(self):
        return hasattr(self.method, 'method_input_argument_names')
    
    @late
    def is_async_supported(self):
        if hasattr(self.method, 'is_async_supported'):
            return self.method.is_async_supported
        elif self.method_is_legacy:
            return True
        else:
            return False
            
    @late
    def legacy_specification(self):
        if self.method_is_code:
            return self.method.legacy_specification
        elif self.method_is_legacy:
            return self.method.specification
        else:
            return None
            
    @late
    def method_input_argument_names(self):
        if self.method_is_code:
            return self.method.method_input_argument_names
        elif self.method_is_legacy:
            return map(lambda x : x.name, self.method.specification.input_parameters)
        else:
            args = inspect.getargspec(self.method).args
            if args:
                if args[0] == 'self' or args[0] == 'cls':
                    return args[1:]
            return args
        
    @late
    def optional_method_input_argument_names(self):
        if self.method_is_code:
            return self.method.optional_method_input_argument_names
        elif self.method_is_legacy:
            return map(lambda x : x.name, self.method.specification.iter_optional_input_parameters())
        else:
            argspec = inspect.getargspec(self.method)
            defaults = argspec.defaults
            if defaults is None or len(defaults) == 0:
                return []
            else:
                return argspec.args[-len(defaults):]
             
      
    @late
    def method_output_argument_names(self):
        if self.method_is_code:
            return self.method.method_output_argument_names
        elif self.method_is_legacy:
            return map(lambda x : x.name, self.method.specification.output_parameters)
        else:
            return ()
           
    @late
    def index_input_attributes(self):
        if self.method_is_code:
            return self.method.index_input_attributes
        else:
            return None
            
    @late
    def nbody_input_attributes(self):
        if self.method_is_code:
            return self.method.nbody_input_attributes
        else:
            return [False] * len(self.method_input_argument_names)
    
    @late
    def index_output_attributes(self):
        if self.method_is_code:
            return self.method.index_output_attributes
        else:
            return None
            

class CodeMethodWrapper(AbstractCodeMethodWrapper):
    
    def __init__(self, method, definition):
        self.method = method
        self.definition = definition
        self.definition.check_wrapped_method(self)
    
    def __call__(self, *list_arguments, **keyword_arguments):
        async_dependency=keyword_arguments.pop("async_dependency", None)
        return_request=keyword_arguments.pop("return_request", False)
        
        if any(isinstance(x, AbstractASyncRequest) for x in list_arguments) or \
           any(isinstance(x, AbstractASyncRequest) for x in keyword_arguments):
            list_arguments_=[]
            keyword_arguments_=dict()
            for arg in list_arguments:
                if isinstance(arg, AbstractASyncRequest):
                        async_dependency=arg.join(async_dependency)              
            for key,arg in keyword_arguments.items():
                if isinstance(arg, AbstractASyncRequest):
                        async_dependency=arg.join(async_dependency)              

            def dummy_factory():
                return FakeASyncRequest()

            # need this step in between to make sure results are available
            request=DependentASyncRequest(async_dependency, dummy_factory)

            def factory():
                list_arguments_=[]
                keyword_arguments_=dict()
                for arg in list_arguments:
                    if isinstance(arg, AbstractASyncRequest):
                        list_arguments_.append(arg.result())
                    else:
                        list_arguments_.append(arg)                    
                for key,arg in keyword_arguments.items():
                    if isinstance(arg, AbstractASyncRequest):
                        keyword_arguments_[key]=arg.result()
                    else:
                        keyword_arguments_[key]=arg
                  
                return self.asynchronous(*list_arguments_, **keyword_arguments_)

            request=DependentASyncRequest(request, factory)
            if return_request:
                request._result_index=self.convert_result_index()
                return request
            else:
                return request.result()

        if async_dependency is not None:
            def factory():
                return self.asynchronous(*list_arguments, **keyword_arguments)
            request = DependentASyncRequest(async_dependency, factory)
            if return_request:
                request._result_index=self.convert_result_index()
                return request
            else:
                return request.result()

        if return_request:
            request = self.asynchronous(*list_arguments, **keyword_arguments)
            request._result_index=self.convert_result_index()
            return request

        object = self.precall()
        list_arguments, keyword_arguments = self.convert_arguments(list_arguments, keyword_arguments)
        result = self.method(*list_arguments, **keyword_arguments)
        
        result = self.convert_result(result)
        
        self.postcall(object)
        
        return result
    
    def asynchronous(self, *list_arguments, **keyword_arguments):
        if not self.is_async_supported:
            raise exceptions.AmuseException("asynchronous call is not supported for this method")
        
        
        object = self.precall()
        
        list_arguments, keyword_arguments = self.convert_arguments(list_arguments, keyword_arguments)
        
        request = self.method.asynchronous(*list_arguments, **keyword_arguments)
        
        def handle_result(function):
            
            result = function()

            result = self.convert_result(result)
        
            self.postcall(object)
                        
            return result
        
        request.add_result_handler(handle_result)
        
        return request
        
        
    def convert_arguments(self, list_arguments, keyword_arguments):
        return self.definition.convert_arguments(self, list_arguments, keyword_arguments)
    
    def convert_result(self, result):
        return self.definition.convert_result(self, result)

    def convert_result_index(self):
        return self.definition.convert_result_index(self)
    
    def precall(self):
        return self.definition.precall(self)
        
    def postcall(self, object):
        self.definition.postcall(self, object)
        

    def __str__(self):
        return 'wrapped<{0}>'.format(self.method)
    
    
class CodeMethodWrapperDefinition(object):
    
    def check_wrapped_method(self, method):
        pass
        
    def precall(self, method):
        return None
        
    def postcall(self, method, object):
        pass
        
    def convert_arguments(self, method, list_arguments, keyword_arguments):
        return (list_arguments, keyword_arguments)
    
    def convert_result(self, method, result):
        return result
    
    
        
class ProxyingMethodWrapper(AbstractCodeMethodWrapper):
    
    def __init__(self, code_interface, attribute_name):
        self.code_interface = code_interface
        self.attribute_name = attribute_name
        self.method = getattr(code_interface, attribute_name)
        
    def __getstate__(self):
        return {
            "code_interface": self.code_interface,
            "attribute_name": self.attribute_name
        }
    
    def __setstate__(self, state):
        self.__dict__.update(state)
        self.method = getattr(self.code_interface, self.attribute_name)
        
    
    def __call__(self, *list_arguments, **keyword_arguments):
        return self.method(*list_arguments, **keyword_arguments)
    
    def asynchronous(self, *list_arguments, **keyword_arguments):
        return self.method.asynchronous(*list_arguments, **keyword_arguments)

    def __str__(self):
        return 'wrapped<{0}>'.format(self.method)

class IncorrectWrappedMethodException(exceptions.AmuseException):
    formatstring = "{0}"

