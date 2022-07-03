.. _asynchronous:

Asynchronous calls
==================

Asynchronous function calls and requests is a mechanism in AMUSE to
let multiple models or code instances do work concurrently.

Normally when a function call to a code is made (be it explicitly or 
implicitly by e.g. referencing a particle attribute) the framework 
script waites on the call to finish. Any method call on an AMUSE code
can be made asynchronous, by appending ``.asynchronous`` to the method's 
name or (preferred) adding a keyword argument ``return=request``.
When such asynchronous call is made to a code, the Python script
continues to run while the call is being made.  The Python script can
then perform computations or calls to other codes in AMUSE.
Subsequent calls to the same code instance can be made, but they
will be performed in sequence.

An asynchronous AMUSE function call returns  immediately,
i.e. it does not wait for the worker code to finish the function.
Instead of the normal return value, the asynchronous function returns a request
object, which can be used to access the result of the function call later.
The request can also be added to a request pool, which manages
several concurrent requests.
::

   request1 = gravity.evolve_model(target_time, return_request=True)
   # ... do something else ...
   print(request1.result()) 

Example
-------

Running two models simultaneously::

    from amuse.community.huayno.interface import Huayno
    from amuse.units import nbody_system
    from amuse.ic.plummer import new_plummer_model
    
    h1=Huayno()
    h2=Huayno()
    
    p1=new_plummer_model(100)
    p2=new_plummer_model(100)
    
    h1.particles.add_particles(p1)
    h2.particles.add_particles(p2)
        
    request1 = h1.evolve_model(1| nbody_system.time, return_request=True)
    request2 = h2.evolve_model(1| nbody_system.time, return_request=True)
    
    # join the request to make request pool that can be waited on
    pool=request1.join(request2)
    
    # wait for the requests to finish
    pool.waitall()
  
Requests can be joined to form a pool. An empty pool can explicitly
be constructed by importing ``from amuse.rfi.async_request import AsyncRequestsPool``.

Asynchronous variable access::
  
    # setting attribute
    # normal synchronous call
    h1.particles[0:10].mass = 1 | nbody_system.mass
    
    # asynchronous call
    h2.particles[0:10].request.mass = 1 | nbody_system.mass
    
    # getting attribute
    # synchronous call, returns an array
    xpos = h1.particles[0:10].x
    
    # asynchronous call, returns a request.
    request = h2.particles.request.x
    xpos = request.result() # retrieve result. Implicit wait for the request to finish
  
Requests can be used in arithmetic operations, and such operation returns a new 
request that can be waited on (with the corresponding result)::

    x1=h1.particles.request.x
    x2=h2.particles.request.x
    
    relx=x2-x1
    
    relx=relx.result()
  

Caveats
-------

* Mixing asynchronous and synchronous calls produces correct results,
  but has consequences for the sequencing of the work: Performing a
  normal, synchronous call to a code after one or more asynchronous
  calls, causes an implicit wait for the asynchronous calls to complete.

* Accessing ``result()`` on a request causes a wait for the request to
  complete.

* When making several calls to the same code instance, the first call
  begins immediately. The second call begins only when the code is waited on,
  not automatically when the first call completes.

* asynchronous access does not work (yet) with compound attributes like ``velocity``.
