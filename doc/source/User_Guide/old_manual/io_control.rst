.. raw:: latex

   \clearpage

.. _io:

I/O Control
===============

Some aspects of Rayleigh's I/O can be controlled through variables found in the io_controls namelist.  This section covers:

`I/O Format Controls`_    

`I/O Redirection`_

I/O Format Controls
*********************


By default, integer output is reported with 8 digits and padded with leading zeros.  This includes integer iteration
numbers reported to stdout at each timestep and integer-number filenames created through diagnostics and checkpointing
output.  If desired, the number of digits may be controlled through the **integer_output_digits** variable.   When
reading in a Checkpoint created with a different number of digits, set the **integer_input_digits** variable to an appropriate
value.  

At several points in the code, floating-point output is sent to stdout.  This output is formatted using scienific notation, with
three digits to the right of the decimal place.     The number of digits after the decimal can be controlled through the
**decimal_places** variable.

As an example, the following combination of inputs
::
   &temporal_controls_namelist
   checkpoint_interval=10
   /
   &io_controls_namelist
   integer_output_digits=5
   integer_input_digits=3
   decimal_places=5
   /
   &initial_conditions_namelist
   init_type=-1
   restart_iter=10
   /

would restart from checkpoint files with the prefix formatted as:

::

   Checkpoints/010_grid_etc.

It would generate status line, shell_slice output, and checkpoints formatted as:

::

   Iteration:  00033   DeltaT:  1.00000E-04   Iter/sec:  2.68500E+00
   Shell_Slices/00020
   Checkpoints/00020_grid_etc.

**Developer's Note:**  The format codes generated through the values of these three variables are declared (with
descriptive comments) in Controls.F90.   For integer variables that may take on a negative value, additional format codes with one extra
digit (for the negative sign) are also provided.


I/O Redirection
*********************

Rayleigh writes all text output (e.g., error messages, iteration
counter, etc.) to stdout by default. Different computing centers handle
stdout in different ways, but typically one of two path is taken. On
some machines, a log file is created immediately and updated
continuously as the simulation runs. On other machines, stdout is
buffered on-node and written to disk only when the run has terminated.

There are situations where it can be advantageous to have a regularly
updated log file whose update frequency may be controlled. This feature
exists in Rayleigh and may be accessed by assigning values to
**stdout_flush_interval** and **stdout_file** in the io controls
namelist.

::

   &io_controls_namelist
   stdout_flush_interval = 1000
   stdout_file = 'routput'
   /

Set stdout_file to the name of a file that will contain Rayleigh’s text
output. In the example above, a file named *routput* will be appear in
the simulation directory and will be updated periodically throughout the
run. The variable stdout_flush_interval determines how many lines of
text are buffered before they are flushed to routput. Rayleigh prints
time-step information during each time step, and so setting this
variable to a relatively large number (e.g., 100+) prevents excessive
disk access from occurring throughout the run. In the example above, a
text buffer flush will occur once 1000 lines of text have been
accumulated.

Changes in the time-step size and self-termination of the run will also
force a text-buffer flush. Unexpected crashes and sudden termination by
the system job scheduler do not force a buffer flush. Note that the
default value of stdout_file is **‘nofile’**. If this value is
specified, output will directed to normal stdout.

To save on disk space for logs of very long runs, the number of status outputs
can be reduced by specifying **statusline_interval** in the
**io_controls_namelist**. This causes only every n-th status line to be
written.
