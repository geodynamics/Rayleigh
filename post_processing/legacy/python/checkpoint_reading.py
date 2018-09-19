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

import numpy as np
import os

#Routines for reading checkpoint grid files and for translating endianness of checkpoint files

class Grid_Info():
    def __init__(self, filename='nothing',byteswap = False):
        self.tpar = 0
        self.gpars = 0
        self.radius = 0
        self.simtime = 0
        if (filename != 'nothing'):
            self.read_grid(filename,bswap=byteswap)
    def read_grid(self,gfile,bswap=False):
        #Grid files are currently written in fortran 77 unformatted style.
        #That means there are several "spacer" 4-byte integers preceeding and following each record
        fd = open(gfile,'rb')
        nr_grid_lmax = []
        grid = Grid_Info()
        # First read dimensions
        for i in range(3):
            sp = swapread(fd,dtype='int32',count=1,swap=bswap)
            val = swapread(fd,dtype='int32',count=1,swap=bswap)
            sp = swapread(fd,dtype='int32',count=1,swap=bswap)
            nr_grid_lmax.append(val)


        #Next get dt, and old dt
        dt_old_dt = []
        sp = swapread(fd,dtype='int32',count=1,swap=bswap)
        val = swapread(fd,dtype='float64',count=1,swap=bswap)
        sp = swapread(fd,dtype='int32',count=1,swap=bswap)
        dt_old_dt.append(val)

        sp = swapread(fd,dtype='int32',count=1,swap=bswap)
        val = swapread(fd,dtype='float64',count=1,swap=bswap)
        sp = swapread(fd,dtype='int32',count=1,swap=bswap)
        dt_old_dt.append(val)



        #Now get the radius array
        nr = nr_grid_lmax[0]
        sp = swapread(fd,dtype='int32',count=1,swap=bswap)
        radius = swapread(fd,dtype='float64',count=nr,swap=bswap, verbose = True)
        sp = swapread(fd,dtype='int32',count=1,swap=bswap)


        self.gpars = nr_grid_lmax
        self.tpars = dt_old_dt
        self.radius = radius

        sp = swapread(fd,dtype='int32',count=1,swap=bswap)
        simtime = swapread(fd,dtype='float64',count=1,swap=bswap)
        sp = swapread(fd,dtype='int32',count=1,swap=bswap)

        self.simtime = simtime
        fd.close()


    def write_grid(self,gfile,bswap=False):
        
        fd = open(gfile,'wb')


        # First write dimensions
        four = np.ndarray(shape=(1),dtype='int32')
        four[0] = 4
        eight = np.ndarray(shape=(1),dtype='int32')
        eight[0] = 8
        for i in range(3):
            val = self.gpars[i]
            swapwrite(four,fd,swap=bswap)
            swapwrite(val,fd,swap=bswap)
            swapwrite(four,fd,swap=bswap)



        #Next write dt, and old dt
        val = self.tpars[0]
        swapwrite(eight , fd,swap=bswap)
        swapwrite(val, fd,swap=bswap)
        swapwrite(eight , fd,swap=bswap)

        val = self.tpars[1]
        swapwrite(eight,fd,swap=bswap)
        swapwrite(val,fd,swap=bswap)
        swapwrite(eight,fd,swap=bswap)


        #Now write the radius array
        nr = self.gpars[0]
        eightnr = np.ndarray(shape=(1),dtype='int32')
        eightnr[0] = 8*nr

        radius = self.radius
        swapwrite(eightnr,fd,swap=bswap)
        swapwrite(radius,fd,swap=bswap,verbose = True, array = True)
        swapwrite(eightnr,fd,swap=bswap)


        #Finally write the simulation time


        swapwrite(eight,fd,swap=bswap)
        swapwrite(self.simtime,fd,swap=bswap)
        swapwrite(eight,fd,swap=bswap)

        fd.close()

def swapread(fd,dtype='float64',count=1,swap=False, verbose = False):
        #simple wrapper to numpy.fromfile that allows byteswapping based on Boolean swap
        # Set swap to true to read bytes in different endianness than current machine
        if (swap):
                if(verbose):
                    print "Swapping on read."
                val = np.fromfile(fd,dtype=dtype,count=count).byteswap()
        else:
                val = np.fromfile(fd,dtype=dtype,count=count)
        if (len(val) == 1):
                val = val[0]
        return val

def swapwrite(val,fd,swap=False,verbose=False, array = False):
        #simple wrapper to numpy.tofile that allows byteswapping based on Boolean swap
        #set swap to true to write bytes in different endianness than current machine

        if (swap):
                if (verbose):
                    print "Swapping on write."
                if (array):
                    if (verbose):
                        print "Swapping entire array of bytes"
                    val2 = val.byteswap().newbyteorder()                
                else:    
                    val2 = val.newbyteorder()
                val2.tofile(fd)
        else:
                val.tofile(fd)


def swap_translate(infile,ofile,val_count,dtype='float64',bswap=False):
    if (bswap):
        #reads infile in opposite endianess
        #Writes contents to ofile in native endianess
        fd = open(infile,'rb')
        val = np.fromfile(fd,dtype=dtype,count=val_count).byteswap()
        fd.close()
        fd = open(ofile,'wb')
        val.tofile(fd)
        fd.close()

    else:

        #Reads infile in native endianess.
        #Writes contents to ofile in opposite endianess
        fd = open(infile,'rb')
        val = np.fromfile(fd,dtype=dtype,count=val_count)
        fd.close()
        fd = open(ofile,'wb')
        val2 = val.byteswap().newbyteorder()
        val2.tofile(fd)
        fd.close()




def translate_checkpoint(iteration,indir,outdir,magnetism=False, tonative= False):
    if (indir == outdir):
        print "Error:  input directory is same as output directory"
        return
    istring = str(iteration)
    #Convert iteration format to i8.8
    for i in range(8):
        if (iteration < 10**i):
            istring = '0'+istring
    

    gridfile = istring+'_grid_etc'
    if (indir != ''):
        last_char = indir[len(indir)-1]
        if (last_char != '/'):
            indir = indir+'/'
        gridfile = indir+gridfile

    ogridfile = istring+'_grid_etc'
    if (outdir != ''):
        last_char = outdir[len(outdir)-1]
        if (last_char != '/'):
            outdir = outdir+'/'
        ogridfile = outdir+ogridfile

    #///////////////////////////////
    # Translate the grid file
    print 'Translating '+gridfile+' to '+ogridfile
    thisgrid = Grid_Info(gridfile,byteswap=tonative)
    print thisgrid.gpars
    thisgrid.write_grid(ogridfile,bswap= not tonative)

    nr = thisgrid.gpars[0]
    lmax = thisgrid.gpars[2]
    nell = lmax+1
    item_count = nell*(nell+1)*nr  # Number of 8-byte reals in each checkpoint file
    fsize = item_count*8   # file size (in bytes)


    #//////////////////////////////////////////
    #Translate the individual checkpoint files

    field_strings = ['W', 'P', 'T', 'Z'] 
    if (magnetism):
        field_strings.append('A')
        field_strings.append('C')
    adds = ['', 'AB']
    numfield = len(field_strings)
    for j in range(numfield):
        for i in range(2):
            suffix = istring+'_'+field_strings[j]+adds[i]
            ifile = indir+suffix
            ofile = outdir+suffix
            print 'Translating '+ifile+' to '+ofile
            swap_translate(ifile,ofile,item_count,dtype='float64',bswap = not tonative)



