"""Classes with plotting methods for equilibrium objects.

Since the equilibrium objects are implemented in Rust, the only reasonable way to simulate
inheritance is to 'attach' these methods to the object right after they are imported from the
interpreter.

Although the plotter classes appear as parent classes in the stub files, this is only done for
building the documentation and has no actual effect on the object.

The methods only require what is guaranteed to be implemented by the corresponding evaluation
trait, so they can be called on all analytical and numerical equilibrium objects.

Some optional parameters that require non-essential attributes are ignored if the attributes
are not found.
"""

__all__ = []
