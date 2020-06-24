.. raw:: latex

   \clearpage

.. _benchmarking:

Running a Benchmark
===================

Rayleigh has been programmed with internal testing suite so that its
results may be compared against benchmarks described in Christensen et al. (2001)
:cite:`CHRISTENSEN200125` and Jones et al. (2011)
:cite:`JONES2011120`

We recommend running a benchmark whenever running Rayleigh on a new
machine for the first time, or after recompiling the code. The
Christensen et al. (2001) :cite:`CHRISTENSEN200125` reference describes two Boussinesq tests that
Rayleigh’s results may be compared against. The Jones et al. (2011) :cite:`JONES2011120`
reference describes anelastic tests. Rayleigh has been tested
successfully against two benchmarks from each of these papers. Input
files for these different tests are enumerated in Table table_benchmark_
below. In addition to the
input files listed in Table table_benchmark_,
input examples appropriate for use as a template for new runs are
provided with the *\_input* suffix (as opposed to the *minimal* suffix.
These input files still have benchmark_mode active. Be sure to turn this
flag off if not running a benchmark.



**Important:** If you are not running a benchmark, but only wish to
modify an existing benchmark-input file, delete the line containing the
text “*benchmark_mode=X*.” When benchmark mode is active, custom inputs,
such as Rayleigh number, are overridden and reset to their
benchmark-appropriate values.

**We suggest using the c2001_case0_minimal input file for installation
verification**. Algorithmically, there is little difference between the
MHD, non-MHD, Boussinesq, and anelastic modes of Rayleigh. As a result,
when installing the code on a new machine, it is normally sufficient to
run the cheapest benchmark, case 0 from Christensen 2001 :cite:`CHRISTENSEN200125`.

To run this benchmark, create a directory from within which to run your
benchmark, and follow along with the commands below. Modify the
directory structure a each step as appropriate:

#. mkdir path_to_my_sim

#. cd path_to_my_sim

#. cp
   path_to_rayleigh/Rayleigh/input_examples/c2001_case0_minimal   main_input

#. cp path_to_rayleigh/Rayleigh/bin/rayleigh.opt   rayleigh.opt (or use
   *ln -s* in lieu of *cp*)

#. mpiexec -np **N** ./rayleigh.opt -nprow **X** -npcol **Y** -nr **R**
   -ntheta **T**

For the value **N**, select the number of cores you wish to run with.
For this short test, 32 cores is more than sufficient. Even with only
four cores, the lower-resolution test suggested below will only take
around half an hour. The values **X** and **Y** are integers that
describe the process grid. They should both be at least 2, and must
satisfy the expression

.. math:: N=X \times Y.

Some suggested combinations are {N,X,Y} = {32,4,8}, {16,4,4}, {8,2,4},
{4,2,2}. The values **R** and **T** denote the number of radial and
latitudinal collocation points respectively. Select either {R,T}={48,64}
or {R,T}={64,96}. The lower-resolution case takes about 3 minutes to run
on 32 Intel Haswell cores. The higher-resolution case takes about 12
minutes to run on 32 Intel Haswell cores.

Once your simulation has run, examine the file
path_to_my_sim/Benchmark_Reports/00025000. You should see output similar
to that presented in Tables table_benchmark_low_ or table_benchmark_high_ . Your numbers may differ
slightly, but all values should have a % difference of less than 1. If
this condition is satisfied, your installation is working correctly.

.. _table_benchmark:

.. centered:: **Table. Benchmark.**

Benchmark-input examples useful for verifying Rayleigh’s installation.
Those from Christensen et al. (2001) :cite:`CHRISTENSEN200125`
are Boussinesq. Those from Jones et al. (2011) :cite:`JONES2011120` are anelastic. Examples are found
in the directory: Rayleigh/input_examples/

+-----------------------+-----------------+--------------------------------+
| Paper                 | Benchmark       | Input File                     |
+=======================+=================+================================+
| Christensen et al.    | Case 0          | c2001_case0_minimal            |
+-----------------------+-----------------+--------------------------------+
| Christensen et al.    | Case 1(MHD)     | c2001_case1_minimal            |
+-----------------------+-----------------+--------------------------------+
| Jones et al. 2011     | Steady Hydro    | j2011_steady_hydro_minimal     |
+-----------------------+-----------------+--------------------------------+
| Jones et al. 2011     | Steady MHD      | j2011_steady_MHD_minimal       |
+-----------------------+-----------------+--------------------------------+


.. _table_benchmark_low:

.. centered:: **Table. Benchmark Low.**

Rayleigh benchmark report for Christensen
et al. (2001) :cite:`CHRISTENSEN200125` case 0 when run with nr=48 and ntheta=64. Run time was
approximately 3 minutes when run on 32 Intel Haswell cores.

Run command:

.. code-block::

 mpiexec -np 32 ./rayleigh.opt -nprow 4 -npcol 8 -nr 48 -ntheta 64

+-----------------+------------+------------+--------------+-----------+
| Observable      | Measured   | Suggested  | % Difference | Std. Dev. |
+=================+============+============+==============+===========+
| Kinetic Energy  | 58.347827  | 58.348000  | -0.000297    | 0.000000  |
+-----------------+------------+------------+--------------+-----------+
| Temperature     | 0.427416   | 0.428120   | -0.164525    | 0.000090  |
+-----------------+------------+------------+--------------+-----------+
| Vphi            | -10.118053 | -10.157100 | -0.384434    | 0.012386  |
+-----------------+------------+------------+--------------+-----------+
| Drift Frequency | 0.183272   | 0.182400   | 0.477962     | 0.007073  |
+-----------------+------------+------------+--------------+-----------+


.. _table_benchmark_high:


.. centered:: **Table. Benchmark High.**

Rayleigh benchmark report for Christensen
et al. (2001) :cite:`CHRISTENSEN200125` case 0 when run with nr=64 and ntheta=96. Run time was
approximately 12 minutes when run on 32 Intel Haswell cores.

Run command:

.. code-block::

  mpiexec -np 32 ./rayleigh.opt -nprow 4 -npcol 8 -nr 64 -ntheta 96

+-----------------+------------+------------+--------------+-----------+
| Observable      | Measured   | Suggested  | % Difference | Std. Dev. |
+=================+============+============+==============+===========+
| Kinetic Energy  | 58.347829  | 58.348000  | -0.000294    | 0.000000  |
+-----------------+------------+------------+--------------+-----------+
| Temperature     | 0.427786   | 0.428120   | -0.077927    | 0.000043  |
+-----------------+------------+------------+--------------+-----------+
| Vphi            | -10.140183 | -10.157100 | -0.166551    | 0.005891  |
+-----------------+------------+------------+--------------+-----------+
| Drift Frequency | 0.182276   | 0.182400   | -0.067994    | 0.004877  |
+-----------------+------------+------------+--------------+-----------+
