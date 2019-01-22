import select
import operator

import channel

class AbstractASyncRequest(object):
    def __nonzero__(self):
        return not self.is_finished

    def waitone(self):
        return self.wait()

    def waitall(self):
        while not self.is_finished:
            self.wait()

    def wait(self):
        raise Exception("not implemented")

    @property
    def is_finished(self):
        return self._is_finished

    @property
    def is_result_set(self):
        return self._is_result_set

    def is_result_available(self):
        raise Exception("not implemented")

    def is_failed(self):
        if not self.is_finished:
            return False
        return not self.is_result_set

    def result(self):
        raise Exception("not implemented")

    @property
    def results(self):
        return [self.result()]

    def add_result_handler(self, function, args = ()): 
        self.result_handlers.append([function,args])

    def is_mpi_request(self):
        return False

    def is_socket_request(self):
        return False
 
    def is_other(self):
        return not self.is_mpi_request() and not self.is_socket_request()
        
    def get_mpi_request(self):
        raise Exception("not implemented")

    def get_socket(self):
        raise Exception("not implemented")
        
    #~ def is_pool(self):
        #~ return False
        
    def join(self, other):
        if other is None:
            return self
        elif isinstance(other, AbstractASyncRequest):
            pool = AsyncRequestsPool()
            pool.add_request(self, lambda x: x.result())
            pool.add_request(other, lambda x: x.result())
        elif isinstance(other, AsyncRequestsPool):
            return other.join(self)
        else:
            raise Exception("error: join only possible with ASyncRequest or Pool")
        return pool

    def waits_for(self):
        if self.is_finished:
            return None
        return self

    def __getitem__(self, index):
        return IndexedASyncRequest(self,index)

    #~ def __getattr__(self, name):
        #~ print name, "<<"
        
    def __add__(self, other):
        return baseOperatorASyncRequest(self,other, operator.add)
    def __radd__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.add(y,x))
    def __sub__(self, other):
        return baseOperatorASyncRequest(self,other, operator.sub)
    def __rsub__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.sub(y,x))
    def __mul__(self, other):
        return baseOperatorASyncRequest(self,other, operator.__mul__)
    def __rmul__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.mul(y,x))
    def __truediv__(self, other):
        return baseOperatorASyncRequest(self,other, operator.truediv)
    def __rtruediv__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.truediv(y,x))
    def __floordiv__(self, other):
        return baseOperatorASyncRequest(self,other, operator.floordiv)
    def __rfloordiv__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.floordiv(y,x))
    def __div__(self, other):
        return baseOperatorASyncRequest(self,other, operator.div)
    def __rdiv__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.div(y,x))
    def __pow__(self, other):
        return baseOperatorASyncRequest(self,other, operator.pow)
    def __rpow__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.pow(y,x))
    def __mod__(self, other):
        return baseOperatorASyncRequest(self,other, operator.mod)
    def __rmod__(self, other):
        return baseOperatorASyncRequest(self,other, lambda x,y: operator.mod(y,x))
    def __neg__(self):
        return baseOperatorASyncRequest(self, None, operator.neg)
    
    def __iter__(self):
        if self._result_index:
            for i in self._result_index:
                yield self[i]
        else:
            yield self
       
    #~ def __call__(self):
        #~ return self.result()
        
class DependentASyncRequest(AbstractASyncRequest):
    def __init__(self, parent, request_factory):
        
        self._result_index=None
        self.request=None
        self.parent=parent
        if isinstance(parent, AsyncRequestsPool):
            self.parent=PoolDependentASyncRequest(parent)
        
        def handler(arg):
            result=arg()
            self.request=request_factory()
            for h in self.result_handlers:
                self.request.add_result_handler(*h)
            return result

        self.parent.add_result_handler(handler)
        
        self.result_handlers = []

    @property
    def is_result_set(self):
        if self.request is None:
            return False
        return self.request.is_result_set

    @property
    def is_finished(self):
        if self.request is None:
            if self.parent.is_finished:
                return True  
            else:
                return False
        
        return self.request.is_finished

    def wait(self):    
        try:
            self.parent.waitall()
        except Exception, ex:
            message=str(ex)
            if not message.startswith("Error in dependent call: "):
                message="Error in dependent call: "+str(ex)
            raise type(ex)(message)
        if self.request is None:
            raise Exception("something went wrong (exception of parent?)")
        self.request.wait()

    def is_result_available(self):
        if self.is_finished:
            if self.request is None:
                return False

        #~ if not self.parent.is_finished:
            #~ return False

        if self.request is None:
            return False
            #~ raise Exception("something went wrong (exception of parent?)")
            
        return self.request.is_result_available()

    def result(self):
        self.wait()

        if not self.request.is_result_set:
            raise Exception("result unexpectedly not available")

        return self.request.result()

    @property
    def results(self):
        return self.parent.results+[self.result()]        

    def add_result_handler(self, function, args = ()):
        if self.request is None:
            self.result_handlers.append([function,args])
        else:
            self.request.add_result_handler(function,args)

    def is_mpi_request(self):
        if self.request is None:
            return self.parent.is_mpi_request()
        else:
            return self.request.is_mpi_request()
    
    def is_socket_request(self):
        if self.request is None:
            return self.parent.is_socket_request()
        else:
            return self.request.is_socket_request()

    def waits_for(self):
        if self.is_finished:
            return None
        if self.request is not None:
            return self.request
        else:
            return self.parent.waits_for()

class PoolDependentASyncRequest(DependentASyncRequest):
    def __init__(self, parent):
        self.parent=parent
        self.request=FakeASyncRequest()        
        self.result_handlers = []

class IndexedASyncRequest(DependentASyncRequest):
    def __init__(self, parent, index):
        self.parent=parent
        self.index=index
        self.request=FakeASyncRequest()        
        self.result_handlers = []
        try:
            self._result_index=parent._result_index[index]
        except:
            self._result_index=None

    def result(self):
        self.wait()
        return self.parent.result().__getitem__(self.index)

class baseOperatorASyncRequest(DependentASyncRequest):
    def __init__(self, first, second, operator):
        self._first=first
        self._second=second
        self._operator=operator
        if isinstance( second, AbstractASyncRequest):
            pool=AsyncRequestsPool(first,second)
            self.parent=PoolDependentASyncRequest(pool)
        else:
            self.parent=first
        self.request=FakeASyncRequest()        
        self.result_handlers = []
        try:
            self._result_index=first._result_index[index]
        except:
            self._result_index=None
                    
    def result(self):
        self.wait()
        first=self._first.result()
        second=self._second.result() if isinstance( self._second, AbstractASyncRequest) else self._second
        if second is None: 
            return self._operator(first)
        return self._operator(first,second)

class ASyncRequest(AbstractASyncRequest):
        
    def __init__(self, request, message, comm, header):
        self.request = request
        self.message = message
        self.header = header
        self.comm = comm
        self._is_finished = False
        self._is_result_set = False
        self._called_set_result = False
        self._result = None
        self.result_handlers = []
        self._result_index=None

    def wait(self):
        if self.is_finished:
            return

        self._is_finished = True
                
        self.request.Wait()
        self._set_result()
    
    def is_result_available(self):
        if self.is_finished:
            return self._is_result_set
        return self.request.Test()

    def get_message(self):
        return self.message
        
    def _set_result(self):
        if self._called_set_result:
            return
        self._called_set_result=True
      
        class CallingChain(object):
            def __init__(self, outer, args,  inner):
                self.outer = outer
                self.inner = inner
                self.args = args
                
            def __call__(self):
                return self.outer(self.inner, *self.args)
                
        self.message.receive_content(self.comm, self.header)
        
        current = self.get_message
        for x, args in self.result_handlers:
            current = CallingChain(x, args, current)
        
        self._result = current()
        
        self._is_result_set = True
        
    def result(self):
        self.wait()

        if not self._is_result_set:
            raise Exception("result unexpectedly not available")

        return self._result

    def is_mpi_request(self):
        return True

class ASyncSocketRequest(AbstractASyncRequest):
        
    def __init__(self, message, socket):
        self.message = message
        self.socket = socket
        
        self._is_finished = False
        self._is_result_set = False
        self._called_set_result = False
        self._result = None
        self.result_handlers = []
        self._result_index=None

    def wait(self):
        if self.is_finished:
            return

        self._is_finished = True

        while True:
            readables, _r, _x = select.select([self.socket], [], [])
            if len(readables) == 1:
                break        
        self._set_result()

    def is_result_available(self):
        if self.is_finished:
            return self._is_result_set
            
        readables, _r, _x = select.select([self.socket], [], [], 0.001)
        
        return len(readables) == 1
            
    def get_message(self):
        return self.message
        
    def _set_result(self):
        if self._called_set_result:
            return
        self._called_set_result=True

        class CallingChain(object):
            def __init__(self, outer, args, inner):
                self.outer = outer
                self.inner = inner
                self.args=args
                
            def __call__(self):
                return self.outer(self.inner, *self.args)
                
        self.message.receive(self.socket)
        
        current = self.get_message
        for x,args in self.result_handlers:
            current = CallingChain(x, args, current)
        
        self._result = current()
        
        self._is_result_set = True
        
    def result(self):
        self.wait()

        if not self._is_result_set:
            raise Exception("result unexpectedly not available")

        return self._result

    def is_socket_request(self):
        return True

class FakeASyncRequest(AbstractASyncRequest):
        
    def __init__(self, result=None):
        self._is_finished = False
        self._is_result_set = False
        self._called_set_result = False
        self._result = None
        self.__result = result
        self.result_handlers = []
        self._result_index=None

    def wait(self):
        if self.is_finished:
            return
        self._is_finished = True                
        self._set_result()
        
    def is_result_available(self):
        return True
            
    def _set_result(self):
        if self._called_set_result:
            return
        self._called_set_result=True

        class CallingChain(object):
            def __init__(self, outer, args,  inner):
                self.outer = outer
                self.inner = inner
                self.args = args
                
            def __call__(self):
                return self.outer(self.inner, *self.args)
                        
        current = lambda : self.__result

        for x, args in self.result_handlers:
            current = CallingChain(x, args, current)
        
        self._result = current()
        
        self._is_result_set = True
        
    def result(self):
        self.wait()

        if not self._is_result_set:
            raise Exception("result unexpectedly not available")

        return self._result

class ASyncRequestSequence(AbstractASyncRequest):
        
    def __init__(self, create_next_request, args = ()):
        self.create_next_request = create_next_request
        self.args = args
        self.index = 0
        self.current_async_request = self.create_next_request(self.index, *self.args)
        self.request = self.current_async_request.request
        self._is_finished = False
        self._is_result_set = False
        self._called_set_result = False
        self._result = None
        self.result_handlers = []
        self._results = []
        self._result_index=None

    @property
    def is_finished(self):
        self._next_request()
        return self.current_async_request is None

    def wait(self):
        if self.is_finished:
            return
            
        self._is_finished=True

        while self.current_async_request is not None:
            self.current_async_request.wait()
        
            self._next_request()

        self._set_result()

    def waitone(self):
        if self.is_finished:
            return

        self.current_async_request.wait()
        
        self._next_request()
        
        if self.current_async_request is None:
            self._is_finished=True
            self._set_result()
        

    def _next_request(self):
        if self.current_async_request is not None and \
           self.current_async_request.is_result_available():
            self._results.append(self.current_async_request.result())
            self.index += 1
            self.current_async_request = self.create_next_request(self.index, *self.args)
            if not self.current_async_request is None:
                self.request = self.current_async_request.request
            else:
                self._set_result()
      

    @property
    def results(self):
        return self._results
            
    def is_result_available(self):
        if self.is_finished:
            return True
        
        self._next_request()
            
        return self.current_async_request is None
        
    def add_result_handler(self, function, args = ()):
        self.result_handlers.append([function,args])
    
    def get_message(self):
        return self._results
        
    def _set_result(self):
        if self._called_set_result:
            return
        self._called_set_result=True

        class CallingChain(object):
            def __init__(self, outer, args,  inner):
                self.outer = outer
                self.inner = inner
                self.args = args
                
            def __call__(self):
                return self.outer(self.inner, *self.args)
                
        current = self.get_message
        for x, args in self.result_handlers:
            current = CallingChain(x, args, current)
        
        self._result = current()
        
        self._is_result_set = True
        
    def result(self):
        self.wait()
        
        if not self._is_result_set:
            raise Exception("result unexpectedly not available")

        return self._result

    def is_mpi_request(self):
        return self.current_async_request.is_mpi_request()

    def is_socket_request(self):
        return self.current_async_request.is_socket_request()

    def waits_for(self):
        return self.current_async_request

class AsyncRequestWithHandler(object):
    
    def __init__(self, pool, async_request, result_handler, args=(), kwargs={}):
        self.async_request = async_request
        if result_handler is None:
            def empty(request):
                return request.result()
            result_handler = empty
        self.result_handler = result_handler
        self.args = args
        self.kwargs = kwargs
        self.pool = pool
    

    def run(self):
        self.result_handler(self.async_request, *self.args, **self.kwargs)
        
class AsyncRequestsPool(object):
    
    def __init__(self, *requests):
        self.requests_and_handlers = []
        self.registered_requests = set([])
        for x in requests:
            self.add_request(x)
        
    def add_request(self, async_request, result_handler = None, args=(), kwargs={}):
        if async_request is None:
            return
        if async_request in self.registered_requests:
            return
            #~ raise Exception("Request is already registered, cannot register a request more than once")
            
        self.registered_requests.add(async_request)
        
        self.requests_and_handlers.append(
            AsyncRequestWithHandler(
                self,
                async_request,
                result_handler,
                args,
                kwargs
            )
        )
    

    def waitall(self):
        while len(self) > 0:
            self.wait()

    def waitone(self):
        return self.wait()
        
    def wait(self):
        
        # TODO need to cleanup this code
        #
        while len(self.requests_and_handlers) > 0:
            requests = [x.async_request.waits_for() for x in self.requests_and_handlers if x.async_request.is_other()]
            indices = [i for i, x in enumerate(self.requests_and_handlers) if x.async_request.is_other()]
            if len(requests) > 0:
                for index, x in zip(indices, requests):
                    x.waits_for().waitone()

                    request_and_handler = self.requests_and_handlers[index]
                    if request_and_handler.async_request.is_result_available():
                        self.registered_requests.remove(request_and_handler.async_request)
                    
                        self.requests_and_handlers.pop(index)
                    
                        request_and_handler.run()
                        break

            requests_ = [x.async_request.waits_for().request for x in self.requests_and_handlers if x.async_request.is_mpi_request()]
            indices_ = [i for i, x in enumerate(self.requests_and_handlers) if x.async_request.is_mpi_request()]
            
            requests=[]
            indices=[]
            for r,i in zip(requests_, indices_):
                if r not in requests:
                    requests.append(r)
                    indices.append(i)
                        
            if len(requests) > 0:
                index = channel.MPI.Request.Waitany(requests)
                
                index = indices[index]
                
                request_and_handler = self.requests_and_handlers[index]
                
                request_and_handler.async_request.waits_for().waitone()  # will set the finished flag
                
                if request_and_handler.async_request.is_result_available():
                    self.registered_requests.remove(request_and_handler.async_request)
                    
                    self.requests_and_handlers.pop(index)
                    
                    request_and_handler.run()
                    break
                            
            sockets_ = [x.async_request.waits_for().socket for x in self.requests_and_handlers if x.async_request.is_socket_request()]
            indices_ = [i for i, x in enumerate(self.requests_and_handlers) if x.async_request.is_socket_request()]

            sockets=[]
            indices=[]
            for r,i in zip(sockets_, indices_):
                if r not in sockets:
                    sockets.append(r)
                    indices.append(i)

            if len(sockets) > 0:
                readable, _, _ = select.select(sockets, [], [])
                indices_to_delete = []
                for read_socket in readable:
                    
                    index = sockets.index(read_socket)
                    
                    index = indices[index]
                    
                    request_and_handler = self.requests_and_handlers[index]

                    request_and_handler.async_request.waits_for().waitone()  # will set the finished flag

                    if request_and_handler.async_request.is_result_available():
  
                        self.registered_requests.remove(request_and_handler.async_request)
                    
                        indices_to_delete.append(index)
                    
                        request_and_handler.run()
                    
                for x in reversed(list(sorted(indices_to_delete))):
                    
                    self.requests_and_handlers.pop(x)
                
                if len(indices_to_delete) > 0:
                    break

    def join(self, other):
        if other is None:
            return self
        elif isinstance(other, AbstractASyncRequest):
            self.add_request(other, lambda x: x.result())
        elif isinstance(other, AsyncRequestsPool):
            for x in other.requests_and_handlers:
                self.add_request(
                    x.async_request,
                    x.result_handler,
                    args = x.args,
                    kwargs = x.kwargs
                    )
        else:
            raise Exception("can only join request or pool")
        return self
            
    def __len__(self):
        return len(self.requests_and_handlers)
        
    def __nonzero__(self):
        return len(self)==0

    def waits_for(self):
        raise Exception("pool has no waits for, should never be called")
