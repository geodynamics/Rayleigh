#
#  Copyright (C) 2018 by the authors of the RAYLEIGH code.
#
#  This file is part of RAYLEIGH.
#
#  RAYLEIGH is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  RAYLEIGH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RAYLEIGH; see the file LICENSE.  If not see
#  <http://www.gnu.org/licenses/>.
#

nqmax = 4000  # The maximum possible quantity code.  Should be consistent with the number in SphericalIO.F90
glob_echecks = (nqmax+1)*[0]  # used to check for duplicates across all quantity codes
def classify_line(codeline):
    """Description:  Classifies a line of code as offset-code, quantity-code, or other
       Inputs:  
                codeline -- string variable containing a line of code
       Outputs:
                result -- integer {2 = offset code line, 1 = quantity code line, 0 = other} """
    test=codeline.split('OFFSET CODE')
    if (len(test) > 1):
        return 2
    test=codeline.split('off')
    if (len(test) > 1):
        return 1
    return 0

def parse_qcode_line(codeline):

    tmp=codeline.split('=')
    t1 = tmp[0].split('::')
    varname = t1[1]

    t2 = tmp[1].split('!')
    t3 = t2[0].split('+')
    quantity_code=int(t3[1])

    tex_code = '--'
    tmp=codeline.split(':tex:')
    if (len(tmp) > 1):
        #tex_code=tmp[1]
        tmp2=tmp[1].split('$')
        #tmp2 = tmp.split('\n')
        tex_code='$'+tmp2[1]+'$'
    return (varname, quantity_code,tex_code)

def parse_offset_line(codeline):
    test=codeline.split('OFFSET CODE')
    tmp = test[0]
    tmp2 = tmp.split('!')
    tmp3 = tmp2[0].split('::')
    #print(tmp)
    #print(tmp2)
    #print(tmp3)
    splits = tmp3[1].split('=')
    #print(splits)
    check = splits[1].split('+')
    this_vname=''.join(splits[0].split(' ')) # offset code name

    
    if (len(check) == 1):
        # This is (hopfully) the first offset (the velocity offset of 0)
        oval=int(splits[1])
        rel_offset=0
        rel_vname='initial'
        #mydict[this_vname]=oval
        #print('vname: ',this_vname)
    else:
        rel_vname=''.join(check[0].split(' '))
        
        #print('\n')
        #print('vname: ',this_vname)
        #print('relative to: ', rel_vname)
        #print('\n')
        #offset=mydict[rel_vname]+int(check[1])
        rel_offset=int(check[1])
        mydict[this_vname]=offset
    #Check to see if we have a
    return (this_vname,rel_offset,rel_vname)

indir='../../../src/Diagnostics'
inprefs=[]
inprefs.append(['velocity_field'     ,2])
inprefs.append(['mass_flux'          ,2])
inprefs.append(['vorticity_field'    ,2])
inprefs.append(['kinetic_energy'     ,2])
inprefs.append(['thermal_field'      ,2])
inprefs.append(['thermal_energy'     ,2])
inprefs.append(['magnetic_field'     ,2])
inprefs.append(['current_density'    ,2])
inprefs.append(['magnetic_energy'    ,2])
inprefs.append(['momentum_equation'  ,2])
inprefs.append(['thermal_equation'   ,2])
inprefs.append(['induction_equation' ,2])
inprefs.append(['amom_equation'      ,2])
inprefs.append(['ke_equation'      ,2])
inprefs.append(['me_equation'      ,2])

#['mass_flux', 'vorticity_field', 'kinetic_energy', 'thermal_field', 'thermal_energy', 'magnetic_field', 'current_density', 'magnetic_energy', 'momentum_equation', 'thermal_equation']
#infiles=['velocity_field_codes.F'] #, 'thermal_field_codes.F', 'vorticity_field_codes.F']


mydict={}
mydict['initial']=0

maxsets=3
#ncol=3
for inpref in inprefs:
    lines=[]
    infile=inpref[0]+'_codes.F'
    filehandle=open(indir+'/'+infile,"r")
    while True:                            # Keep reading forever
        theline = filehandle.readline()   # Try to read next line
        if len(theline) == 0:              # If there are no more lines
            break                          #     leave the loop
        lines.append(theline)
    filehandle.close()

    
    varnames = []
    qcodes = []
    texcodes=[]

    loc_echecks = (nqmax+1)*[0]  # Used to check for duplicated indices within each .F file 
    for aline in lines:
        ltype=classify_line(aline)
        if (ltype == 1):
            #print('aline: ', aline)
            pieces=parse_qcode_line(aline)
            vname = pieces[0]
            qcode = pieces[1]
            tcode = pieces[2]
            varnames.append(vname)
            qcodes.append(qcode+offset)
            texcodes.append(tcode)
            loc_echecks[int(qcode+offset)]+=1
            glob_echecks[int(qcode+offset)]+=1
        if (ltype == 2):
            pieces=parse_offset_line(aline)
            vname=pieces[0]
            roffset=pieces[1]
            relvname=pieces[2]
            mydict[vname]=mydict[relvname]+roffset
            offset=mydict[vname]

    for ecvi, ecv in enumerate(loc_echecks):
        if (ecv > 1):
            print("Error in ", infile, ' qcode duplicated: ', ecvi) 
    
    tfcnt=0
    ofile=inpref[0]+'_table'+str(tfcnt)+'.tex'
    myfile= open(ofile,"w")
    nrows = 0
    ncol=inpref[1]
    maxsets=inpref[1]
    print('ncol = ', ncol)
    absmaxrows=20
    nvars = len(varnames)
    maxrows=nvars//ncol
    remain = nvars%ncol
    if (remain != 0):
        maxrows=maxrows+1
    if (maxrows > absmaxrows):
        maxrows=absmaxrows
    # I need to think about logic to keep the columns compact here
    # Will revisit this.
    #LATEST NOTE:  I think this should be split into "chunks" of length absmaxrows and then
    #              the logic above should be carried out on the last chunk
        
 
    outlines=[]
    estring=' \\\\'+'[10pt] \n'
    for i,v in enumerate(varnames):
        vname='\\_'.join(v.split('_'))
        outline = ' '+texcodes[i]+' & '+ str(qcodes[i]) +' & ' + vname
        if (nrows < maxrows):
            outlines.append(outline)
        else:
            outlines[nrows % maxrows] = outlines[i % maxrows]+' & '+outline


        nrows+=1
        if (nrows == maxrows*maxsets) or (i == (len(varnames)-1)):
            for outline in outlines:
                myfile.write(outline+estring)
            outlines=[]
            nrows=0
            tfcnt+=1
            ofile=inpref[0]+'_table'+str(tfcnt)+'.tex'
            myfile.close()
            if (i != (len(varnames)-1)):
                myfile= open(ofile,"w")
    #print(len(varnames))
    myfile.close()

#Finally, perform a consistency check across all quantity codes parsed
for ecvi, ecv in enumerate(glob_echecks):
    if (ecv > 1):
        print('Global Error -- qcode duplicated: ', ecvi) 
