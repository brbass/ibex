----
ibex
----

The ibex code implements meshless methods for the solution of the neutron transport equation and heat conduction equation. The discretization methods include the strong-form collocation method and the weak-form meshless-local Petrov-Galerkin (MLPG) method.

The code is hosted at https://github.com/brbass/ibex.

For the theory behind the code, see https://github.com/brbass/umich_dissertation. 

------------
Dependencies
------------

As mentioned in the INSTALL file, the only dependency of the ibex code that is not included in the source files is Trilinos (https://trilinos.org/).

The external dependencies that are included in the source code include all the files in the "packages/external" folder:

- Eigen (http://eigen.tuxfamily.org), MLP2 license
- nanoflann (https://github.com/jlblancoc/nanoflann), BSD license
- pugixml (https://pugixml.org), MIT license
- quadrule (http://people.sc.fsu.edu/~jburkardt/cpp_src/quadrule/quadrule.html), LGPL license

The individual files contain the full license information.

See "packages/angular/LDFE_Quadrature.hh" for information on the LDFE quadrature, which was provided by Josh Jarrell. 
