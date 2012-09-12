import inspect

from amuse.units import units
from amuse.datamodel import Particle, Particles
from amuse.support.exceptions import AmuseException


class CollisionHandler(object):
    """
    Generic class for handling stellar collisions.
    
    
    """

    def __init__(
            self, 
            collision_code, 
            collision_code_arguments = dict(),
            collision_code_parameters = dict(),
            gravity_code = None, 
            stellar_evolution_code = None, 
            verbose = False,
            **options
        ):
        
        if collision_code.stellar_evolution_code_required and stellar_evolution_code is None:
            raise AmuseException("{0} requires a stellar evolution code: "
                "CollisionHandler(..., stellar_evolution_code=x)".format(collision_code.__class__.__name__))
        if collision_code.gravity_code_required and gravity_code is None:
            raise AmuseException("{0} requires a gravity code: "
                "CollisionHandler(..., gravity_code=x)".format(collision_code.__class__.__name__))
        
        self.collision_code = collision_code
        self.collision_code_arguments = collision_code_arguments
        self.collision_code_parameters = collision_code_parameters
        
        self.gravity_code = gravity_code
        self.stellar_evolution_code = stellar_evolution_code
        self.verbose = verbose
        self.options = options
    
    def handle_collisions(self, primaries, secondaries):
        result = Particles()
        for primary, secondary in zip(primaries.as_set(), secondaries.as_set()):
            result.add_particle(self.handle_collision(primary, secondary))
        return result
    
    def handle_collision(self, primary, secondary):
        colliders = (primary + secondary).copy_to_memory()
        
        if self.verbose:
            print "Handling collision between stars with masses {0}.".format(colliders.mass)
        
        if inspect.isclass(self.collision_code):
            collision_code = self.collision_code(**self.collision_code_arguments)
            if hasattr(collision_code, "initialize_code"):
                collision_code.initialize_code()
            for par, value in self.collision_code_parameters.iteritems():
                setattr(collision_code.parameters, par, value)
            if hasattr(collision_code, "commit_parameters"):
                collision_code.commit_parameters()
        else:
            collision_code = self.collision_code
        
        handle_collision_args = self.options.copy()
        if collision_code.stellar_evolution_code_required:
            handle_collision_args["stellar_evolution_code"] = self.stellar_evolution_code
        if collision_code.gravity_code_required:
            handle_collision_args["gravity_code"] = self.gravity_code
        
        merge_product = collision_code.handle_collision(primary, secondary, **handle_collision_args)
        
        if self.verbose:
            print "{0} concluded with return value:\n{1}".format(self.collision_code.__class__.__name__, merge_product)
        
        if not self.stellar_evolution_code is None:
            if (hasattr(self.stellar_evolution_code, "new_particle_from_model") and 
                (hasattr(merge_product, "internal_structure"))):
                self.stellar_evolution_code.new_particle_from_model(
                    merge_product.internal_structure(), 
                    0.0 | units.Myr, 
                    key = merge_product.key
                )
            else:
                self.stellar_evolution_code.particles.add_particle(merge_product)
            self.stellar_evolution_code.particles.remove_particles(colliders)
            if self.verbose:
                print "Colliders have been replaced by merge product in {0}.".format(self.stellar_evolution_code.__class__.__name__)
        
        if not self.gravity_code is None:
            new_grav_particle = Particle(key=merge_product.key)
            new_grav_particle.mass = merge_product.mass
            
            if hasattr(merge_product, "radius"):
                new_grav_particle.radius = merge_product.radius
            elif hasattr(colliders, "radius"):
                new_grav_particle.radius = max(colliders.radius)
            
            if hasattr(merge_product, "x"):
                new_grav_particle.position = merge_product.position
            else:
                new_grav_particle.position = colliders.center_of_mass()
            
            if hasattr(merge_product, "vx"):
                new_grav_particle.velocity = merge_product.velocity
            else:
                new_grav_particle.velocity = colliders.center_of_mass_velocity()
            
            self.gravity_code.particles.add_particle(new_grav_particle)
            self.gravity_code.particles.remove_particles(colliders)
            if self.verbose:
                print "Colliders have been replaced by merge product in {0}.".format(self.gravity_code.__class__.__name__)
        
        return merge_product
    


