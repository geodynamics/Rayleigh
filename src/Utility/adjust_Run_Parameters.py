"""
Adjust lines in Run_Parameters.F90 that need to be split over multiple
lines with proper Fortran line continuation characters

Usage:
    adjust_Run_Parameters.py <file> <outputname>

    <file> is the Run_Parameters.F90 source file to parse/adjust
    <outputname> is name of new file that holds the changes

"""
from __future__ import print_function

def main(filename, outputname):
    """
    fix/adjust the Run_Parameters.F90 file
    """
    with open(outputname, 'w') as new:
        with open(filename, 'r') as mf:
            for line in mf:
                if (len(line) > 85 and "character" in line.lower() and "::" in line):
                    if (line.endswith("\n")):   # strip trailing newline
                        line = line[:-1]
                    i = line.rfind("//Char(0)") # strip trailing "Char(0)"
                    line = line[:i]
                    first_split = 75            # make first continuation at this column
                    n = 50                      # split into sections of length n
                    ir = line.rfind('"')
                    while (first_split >= ir):  # make sure first split occurs between quotes
                        first_split -= 5
                    continuation = "\"// &\n"
                    prefix = "\t   \""
                    left  = line[:first_split]
                    right = line[first_split:]
                    sections = [right[i:i+n] for i in range(0, len(right), n)]
                    sections[0] = prefix + sections[0]
                    joint = continuation+prefix
                    newline = joint.join(sections)
                    newline = left + continuation + newline
                    new.write(newline+"//Char(0)\n") # add back the newline and Char(0)
                else:
                    new.write(line)

if __name__=="__main__":
    import sys
    main(sys.argv[1], sys.argv[2])

