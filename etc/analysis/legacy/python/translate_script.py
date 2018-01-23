import numpy as np
import os
from checkpoint_reading import *

iteration = 2
#1.
#translates iteration 2 checkpiont files FROM the native endianness

translate_checkpoint(iteration, '','translation')

#2. 
#translates iteration 2 checkpiont files TO the native endianness
#translate_checkpoint(iteration, '','translation',tonative=True)

