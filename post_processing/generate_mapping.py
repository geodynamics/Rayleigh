"""
Generate the look-up-table for Rayleigh diagnostic output quantities.

The look-up-table maps the quantity code to the Fortran variable name associated
with that quantity code. The default behavior is to generate and store the mapping
to the file called lut_mapping.py. The main user interface is located in the file
called lut.py, which will import the lut_mapping.py data. Rayleigh does *NOT* need
to be compiled in order to build the mapping data.

If any custom diagnostic outputs will be defined using the "--with-custom=<CUSTOM ROOT>"
option in the configure script, then these quantities can be included in the mapping
by specifying the custom root location using the "--custom-dir" option. When a custom
directory is provided, the generated map will be stored in a file called
lut_mapping_custom.py and the lut.py user interface will import the custom map data.

Usage:
    generate_mapping.py [options] <Rayleigh_top_directory>

    <Rayleigh_top_directory> is where Rayleigh was installed. This location should
        contain the following files: configure, LICENCE, INSTALL, Makefile, etc.

Options:
    --custom-dir=<c>  Specify the <CUSTOM ROOT> directory to include custom diagnostics
    --overwrite       Overwrite existing lut mapping file [default: False]
"""
from __future__ import print_function
import os
import datetime
import subprocess

def run_command(cmd):
    """
    Run a command in the shell and return the output

    cmd should be a list
    """
    kw = {"stdout":subprocess.PIPE, "stderr":subprocess.PIPE}

    p = subprocess.Popen(cmd, **kw) # run command and store output/error
    out, err = p.communicate()

    try: # convert to standard string
        out = out.decode('utf-8')
        err = err.decode('utf-8')
    except:
        out = out
        err = err

    return out

def find_repo_attributes(directory):
    """Find the git commit/url/branch information"""
    cwd = os.getcwd()
    os.chdir(directory)

    # find some git characteristics
    commit = run_command(["git", "rev-parse", "HEAD"])
    url = run_command(["git", "config", "remote.origin.url"])
    branch = run_command(["git", "rev-parse", "--abbrev-ref", "HEAD"])

    os.chdir(cwd)

    return commit, url, branch

def substring_indices(line, substr):
    """
    Find all indices where the substring is found in the given line

    For example:
        line = "hello,  hello,   hello"
        ind = substring_indices(line, "he")

        should produce: ind = [0,8,17]

    Args
    ----
    line : str
        The line to search
    substr : str
        The substring to find

    Returns
    -------
    indices : list of int
        The indices of each match
    """
    indices = []
    n = len(line)
    for i in range(n):
        ind = line.find(substr, i) # find 'substr' in line[i:]

        # ind = -1 for no match
        if (ind >= 0 and (ind not in indices)): # avoid double counting substrings
            indices.append(ind)

    return indices

def _detexify(line):
    """
    Make the non-math-mode LaTeX compile by fixing special characters

    Args
    ----
    line : str
        The line that has special characters that need to be fixed

    Returns
    -------
    repaired : str
        The repaired line
    """
    line = line.replace("{", "\{") # fix old brackets first
    line = line.replace("}", "\}")
    line = line.replace("^", "\^{}") # now add in new brackets
    line = line.replace("_", "\_")
    line = line.replace("<", "$<$")
    line = line.replace(">", "$>$")

    return line

def _ensure_texable(line, sep=":tex:", verbose=True):
    """
    Verify that each line will compile with LaTeX (this is mostly not required anymore)

    Args
    ----
    line : str
        The line containing the contents ".... :tex: .... $formula$ ...."
    sep : str, optional
        The string that separates the LaTeX formula from everything else
    verbose : bool, optional
        Print status information

    Returns
    -------
    repaired : str
        The repaired line
    """
    entry = line.split(sep)[1].strip()

    if (entry.count("$") % 2 == 1):
        entry = entry + "$"
        if (verbose):
            print("...fixed a :tex: line (missing '$') = {}".format(entry))

    # line is "......$math$......" and there ar "_" or "{" in the 1st/2nd part
    if (not (entry.lstrip().startswith("$") and entry.rstrip().endswith("$"))):
        lind = entry.find("$")
        rind = entry.rfind("$")
        first = _detexify(entry[:lind])
        second = _detexify(entry[rind+1:])
        entry = first + entry[lind:rind+1] + second
        if (verbose):
            print("...fixed a :tex: line (special characters) = {}".format(entry))

    return entry

class Quantity:
    """
    A Rayleigh output quantity

    Attributes
    ----------
    code : int
        The quantity code
    name : str
        The variable name
    tex : str
        The LaTeX formula, including "$" symbols
    """

    def __init__(self, code, name, tex=None):
        """
        Args
        ----
        code : int
            The quantity code
        name : str
            The variable name
        tex : str, optional
            The LaTeX formula, including "$" symbols
        """
        self.code = code
        self.name = name.lower()
        self.tex = tex

class OutputQuantities:
    """
    A collection of Rayleigh output quantities found by parsing Diagnostic_Base.F90

    Attributes
    ----------
    quantities : list of Quantity objects
        The available quantities
    diagnostic_types : dict
        The available quantities sorted by diagnostic type. The keys are the available
        diagnostic types, e.g., "Velocity_Field", "Energies", etc. The value is a list of
        Quantity objects associated with this type.
    """

    def __init__(self, rayleigh_dir, default_location=True, custom_dir=None):
        """
        Args
        ----
        rayleigh_dir : str
            Path to the top of the Rayleigh source tree
        default_location : bool, optional
            If True, the search location will be rayleigh_dir/src/Diagnostics/. If False,
            the search location is just rayleigh_dir/
        custom_dir : str, optional
            Path to the CUSTOM ROOT directory with user defined Diagnostics files
        """
        rayleigh_dir = os.path.expanduser(rayleigh_dir) # expand "~" into $HOME

        self.rayleigh_dir = os.path.abspath(rayleigh_dir)
        if (default_location):
            self.diag_dir = os.path.join(self.rayleigh_dir, "src", "Diagnostics")
        else:
            self.diag_dir = self.rayleigh_dir

        # use "standard" location
        basefile = "Diagnostics_Base.F90"
        self.diag_base = os.path.join(self.diag_dir, basefile)

        # if custom basefile was modified, use it
        if (custom_dir is not None):
            custom_dir = os.path.expanduser(custom_dir) # expand "~" into $HOME
            self.custom_dir = os.path.abspath(custom_dir)

            cfiles = os.listdir(self.custom_dir) # list of filenames, no path info
            if (basefile in cfiles):
                self.diag_base = os.path.join(self.custom_dir, basefile) # add full path
        else:
            self.custom_dir = None

        # main storage structures
        self.quantities = [] # list of available Quantity objects

        # quantity codes organized by where they are defined
        # key = diagnostic_type, value = collection of Quantity objects
        self.diagnostic_types = {}

        # parsing tools
        self.offsets = {} # key=string name, value=integer value

        # parse the various elements
        self._parse_basefile()
        self._parse_diagnostic_files()

    def _parse_quantity_code(self, name, code):
        """evaluate right hand side of "name = offset + value" entries"""

        if (("off" in name) and ("off" not in code)): # case: "offset = number"
            code = int(code)
            if (name not in self.offsets.keys()):
                self.offsets[name] = code # this is an offset definition, save it

        elif (("off" not in name) and ("off" in code)): # case: name = offset + number"
            vals = code.split("+")
            off = self.offsets[vals[0].strip()] # get existing offset
            val = int(vals[1].strip())

            # compute code
            code = off + val

        elif (("off" in name) and ("off" in code)): # case: offset = offset + number"
            vals = code.split("+")
            off = self.offsets[vals[0].strip()] # get existing offset
            val = int(vals[1].strip())

            code = off + val
            if (name not in self.offsets.keys()):
                self.offsets[name] = code # new offset definition, save it

        else: # case: name = number
            code = int(code)

        return code

    def _parse_line(self, Line):
        """parse a line of the form: Integer, parameter :: variable_name = index (! comments)"""
        line = Line.lower()

        result = None

        # ignore empty lines and comments
        if (line.lstrip().startswith("!") or (line.strip() == "")): return result

        # valid lines include all three
        if (("integer" in line) and ("parameter" in line) and ("=" in line)):
            line = line.strip()
            Line = Line.strip() # to maintain case sensitivity of comments

            quant = line.split("::")[1] # everything to right of "::"
            inds = quant.split("!")     # split trailing comments, if any
            quantity = inds[0]          # "name = index" part
            if (len(inds) == 1):
                comment = ''
            else:
                _q = Line.split("::")[1] # restore/maintain case sensitivity
                comment = (_q.split("!")[1]).strip()

            q = quantity.split("=") # parse out the name and index/code
            name = q[0].strip()
            code = q[1].strip()

            code = self._parse_quantity_code(name, code) # convert code to integer

            # ensure LaTeX in the comment will compile; matching $, etc.
            if (("off" not in name) and (":tex:" in comment)):
                comment = _ensure_texable(comment, verbose=False)

            result = (name, code, comment)

        return result

    def _parse_basefile(self):
        """parse the Diagnostic_Base.F90 file for valid quantity codes"""

        with open(self.diag_base, "r") as f:
            base_lines = f.readlines()

        quantities = []

        for Line in base_lines: # loop over base file
            line = Line.lower()

            if (line.lstrip().startswith("include")): # parse the included file

                # included filename is between quotes, so remove those
                if ("'" in line):
                    inc_file = Line.split("'")[1]
                elif ('"' in line):
                    inc_file = Line.split('"')[1]
                else:
                    # this is probably true...
                    raise ValueError("This include line will not compile: {}".format(Line))

                fname = None

                # check "standard" location for the include file
                if (os.path.isfile(os.path.join(self.diag_dir, inc_file))):
                    fname = os.path.join(self.diag_dir, inc_file)

                # use custom location instead, if include file is there
                if (self.custom_dir is not None):
                    if (os.path.isfile(os.path.join(self.custom_dir, inc_file))):
                        fname = os.path.join(self.custom_dir, inc_file)

                if (fname is None):
                    raise ValueError("Could not find the included file: {}".format(inc_file))

                # parse the include file
                with open(fname, "r") as mf: # parse file
                    for l in mf:
                        Q = self._parse_line(l)
                        if (Q is not None): quantities.append(Q)

            else: # no fancy include, so parse the line
                Q = self._parse_line(line)
                if (Q is not None): quantities.append(Q)

        quantities.sort(key=lambda x: x[1]) # sort by quantity code

        for q in quantities: # store results as Quantity objects
            Q = Quantity(q[1], q[0], tex=q[2])
            self.quantities.append(Q)

    def _parse_diagnostic_files(self):
        """parse the Diagnostic_<type>.F90 files to find where each quantity is defined"""

        # find all Diagnostics_....F90 files, will not include the base directory
        files = os.listdir(self.diag_dir)
        files = [f for f in files if "diagnostics" in f.lower()]
        l_files = [f.lower() for f in files]

        # dont bother parsing some files: they are known not to contain definitions
        ignore_files = ["diagnostics_base.f90", "diagnostics_interface.f90",
                        "diagnostics_adotgradb.f90", "diagnostics_mean_correction.f90"]
        for x in ignore_files:
            if (x in l_files):
                ind = l_files.index(x)
                del files[ind]
                del l_files[ind] # required because we search based on l_files to delete from files

        # ensure only fortran files will be parsed
        good_ext = [".F90", ".f90", ".f", ".F", ".F03", ".f03", ".F08", ".f08"]
        files = [f for f in files if os.path.splitext(f)[1] in good_ext]

        # get list of custom files
        if (self.custom_dir is not None):
            cfiles = os.listdir(self.custom_dir)
        else:
            cfiles = []

        # add parent directory, either "standard" location or custom location
        for i,f in enumerate(files):
            if (f not in cfiles):
                parent_dir = self.diag_dir
            else:
                parent_dir = self.custom_dir

            files[i] = os.path.join(parent_dir, f)

        quantity_names = [q.name for q in self.quantities] # all available quantity names

        # loop over each file and determine what quantities are computed here
        for f in files:
            diag_quants = []
            with open(f, "r") as mf:
                for Line in mf:
                    line = Line.lower()
                    if (line.startswith("!") or (line.strip() == "")): continue

                    # determine if this line constitutes a quantity definition
                    quantities = self._find_quantities(line)
                    if (quantities is not None):
                        for q in quantities:
                            if (q not in diag_quants):
                                diag_quants.append(q)

            # ensure unique entries
            diag_quants = list(set(diag_quants))

            # diagnostic type from filename, strip off ".F90": /path/Diagnostics_<type>.F90
            diag_type = os.path.basename(f).split("_", 1)[1] # get <type>.F90
            diag_type = os.path.splitext(diag_type)[0]       # strip off extension

            # allocate space
            if (diag_type not in self.diagnostic_types.keys()):
                self.diagnostic_types[diag_type] = []

            # get names already associated with this diagnostic type
            qnames = [x.name for x in self.diagnostic_types[diag_type]]

            # ensure unique entries
            for q in diag_quants: # loop over found quantities
                if (q not in quantity_names):
                    continue

                if (q in qnames): continue # already added

                ind = quantity_names.index(q) # find associated Quantity object
                Q = self.quantities[ind]
                self.diagnostic_types[diag_type].append(Q)

        # remove any empty entries
        keys = list(self.diagnostic_types.keys())
        for k in keys:
            if (len(self.diagnostic_types[k]) == 0):
                del self.diagnostic_types[k]

        # sort entries by quantity code
        for k in self.diagnostic_types.keys():
              self.diagnostic_types[k].sort(key=lambda x: x.code)

    def _find_quantities(self, line):
        """find all instances of "compute_quantity(Q)" in the line"""
        func_name = "compute_quantity"
        if (func_name not in line): return None

        length = len(func_name)

        quants = []
        indices = substring_indices(line, func_name)
        for ind in indices:
            # line looks like "compute_quantity (  var_name   )", extract var_name
            start = ind + length
            open_paren = line[start:].find("(") # find opening parenthesis
            close_paren = line[start:].find(")") # find closing parenthesis

            var_name = line[start+open_paren+1:start+close_paren].strip()

            quants.append(var_name)

        return list(set(quants))

if __name__ == "__main__":
    from docopt import docopt
    args = docopt(__doc__)

    search_path = args['<Rayleigh_top_directory>']
    overwrite = args['--overwrite']
    custom_dir = args['--custom-dir']

    my_path = os.path.dirname(os.path.realpath(__file__))
    if (custom_dir is not None):
        output = os.path.join(my_path, "lut_mapping_custom.py")
    else:
        output = os.path.join(my_path, "lut_mapping.py")

    if (os.path.isfile(output) and (not overwrite)):
        print("\nOutput file already exists: {}".format(output))
        print("\n\tTo overwrite this file use the \"--overwrite\" argument\n")

        import sys; sys.exit()

    # build a header with a doc string and some python executable code
    header = "\"\"\"\n" +\
        "DO NOT EDIT THIS FILE, unless you really know what you are doing...\n\n" +\
        "This file is automatically generated by generate_mapping.py. To make\n" +\
        "changes to this file, do so by re-running that code on the\n" +\
        "appropriate directories in the Rayleigh source code.\n\n" +\
        "There are three dictionaries defined here:\n\n" +\
        "    name_given_code  --- key is integer quantity code, value is string name\n" +\
        "    code_given_name  --- key is string name, value is integer quantity code\n" +\
        "    tex_given_code   --- key is integer quantity code, value is string LaTeX\n" +\
        "@@CUSTOM@@" +\
        "\"\"\"\n" +\
        "from collections import OrderedDict\n\n" +\
        "name_given_code = OrderedDict()\n" +\
        "code_given_name = OrderedDict()\n" +\
        "tex_given_code  = OrderedDict()\n\n"

    footer = "\n"

    # add more information to the custom header
    cheader = "\nAutomatically generated on @@DATE@@\n" +\
              "\nRayleigh information:\n\n" +\
              "\tcommit : @@COMMIT@@\n" +\
              "\t  url  : @@URL@@\n" +\
              "\tbranch : @@BRANCH@@\n" +\
              "\nCustom information:\n\n" +\
              "\t  dir  : @@DIR@@\n\n"

    # parse for the output quantities
    outputQ = OutputQuantities(search_path, default_location=True, custom_dir=custom_dir)
    offsets = list(outputQ.offsets.keys())

    # get Rayleigh information
    commit, url, branch = find_repo_attributes(search_path)

    # build header
    if (custom_dir is not None):
        today = datetime.date.today()
        date = today.strftime("%b-%d-%Y") # format as "May-8-2021"

        cheader = cheader.replace("@@DATE@@", date)
        cheader = cheader.replace("@@COMMIT@@", commit)
        cheader = cheader.replace("@@URL@@", url)
        cheader = cheader.replace("@@BRANCH@@", branch)
        cheader = cheader.replace("@@DIR@@", os.path.abspath(os.path.expanduser(custom_dir)))

        header = header.replace("@@CUSTOM@@", cheader)
    else:
        header = header.replace("@@CUSTOM@@", "\n")

    if (overwrite):
        print("\nOverwriting quantity code mapping ...")
    else:
        print("\nWriting quantity code mapping ...")
    with open(output, "w") as f:
        f.write(header)
        for Q in outputQ.quantities:
            if (Q.name in offsets): continue

            # some LaTeX has single quotes, so use double here
            line1 = "name_given_code[{}] = \"{}\"\n"
            line2 = "code_given_name[\"{}\"] = {}\n"
            line3 = "tex_given_code[{}] = r\"{}\"\n\n"

            n = Q.name; c = Q.code; t = Q.tex
            f.write(line1.format(c, n))
            f.write(line2.format(n, c))
            f.write(line3.format(c, t))

        f.write(footer)

    print("\nSaved mapping to: {}\n".format(output))

