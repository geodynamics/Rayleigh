"""
Module to help manage mapping between quantity names and quantity codes
"""
from __future__ import print_function

# if custom mapping exists, use it
try:
    from lut_mapping_custom import code_given_name, name_given_code, tex_given_code
except ImportError:
    from lut_mapping import code_given_name, name_given_code, tex_given_code

from lut_shortcuts import shortcuts

def lookup(quantity, __check_shortcuts=True):
    """
    Search the look up table for the quantity. If the quantity code was
    provided, then return the name. If the name was given, return the code.

    __check_shortcuts is an internally used keyword and should always be set to True
    """
    code_given = False
    name_given = False

    try: # assume integer-like
        code = int(quantity)
        code_given = True
    except:
        name = quantity
        name_given = True

    if (code_given):
        if (code in name_given_code.keys()): # code given, grab name
            return name_given_code[code]
        else:
            return None

    if (name_given):
        if (__check_shortcuts): # check shortcuts first
            ind = shortcut_lookup(name)

            if (ind is not None): # found in shortcuts
                return ind

        name = name.lower()
        if (name in code_given_name.keys()): # name give, grab code
            return code_given_name[name]
        else:
            return None

    return None

def shortcut_lookup(name):
    """
    Search the user defined shortcuts for the given quantity name, return the code
    """
    if (name not in shortcuts.keys()):
        return None

    qname = shortcuts[name]
    ind = lookup(qname, __check_shortcuts=False)

    return ind

def latex_formula(qcode):
    """
    Return the LaTeX formula for the given quantity code
    """
    if (qcode in tex_given_code.keys()):
        return tex_given_code[qcode]
    return None

def parse_quantity(quantity):
    """
    Given a quantity, return the index and the name

    Examples:
        parse_quantity("v_r") returns (1, "v_r")
        parse_quantity(  1  ) returns (1, "v_r")
        parse_quantity( "1" ) returns (1, "v_r")
    """
    try:
        code = int(quantity)
        name = lookup(code)
    except:
        name = quantity
        code = lookup(name)
    return code, name

def parse_quantities(quantities):
    """
    Given a list of quantities, return the indices and the names

    Examples:
        parse_quantity(["v_r", 2, "3"]) returns [1,2,3], ["v_r", "v_theta", "v_phi"]
    """
    codes = []; names = []

    for q in quantities:
        c, n = parse_quantity(q)
        codes.append(c)
        names.append(n)

    return codes, names

def find_possible(search_string):
    """
    Search the quantity names for the given string. Return possible codes & names
    """
    codes = []; names = []
    search_string = search_string.lower()
    for c,n in name_given_code.items():

        if (search_string in n):
            codes.append(c)
            names.append(n)

    return codes, names

def quantity_available(quantity):
    """
    Check if the quantity is available, returns True/False
    """
    available = False

    try:
        code = int(quantity)
        if (code in name_given_code.keys()):
            available = True
    except:
        name = quantity.lower()
        if (name in code_given_name.keys()):
            available = True

    return available

def quantities_available(quantities):
    """
    Check if multiple quantities are available, returns list of True/False
    """
    available = []
    for q in quantities:
        available.append(quantity_available(q))
    return available

