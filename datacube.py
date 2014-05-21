#! /usr/bin/env python
import numpy as np
import pyfits
import sys
import glob

##########
### DataObject
##########
class DataObject:
    # Construct DataObject from filename
    def __init__(self,data,header,normalize=False):
        #self.filename=filename
        #self.data,self.header = pyfits.getdata(filename,0,header=True)
        self.data = data
        self.header = header

        if normalize:
            # Normalize by exptime
            self.data /= np.float(self.header['EXPTIME'])/1000.0 # ms

    def __call__(self,*args):
        return self

    @classmethod
    def fromfilename(cls,filename,dext=0,hext=0,normalize=False):
        f = pyfits.open(filename)
        #data,header = pyfits.getdata(filename,0,header=True)
        data = f[dext].data
        header = f[hext].header
        obj = cls(data,header,normalize=normalize)
        obj.filename = filename
        return obj
        

##########
### DataCube
##########
class DataCube:
    @staticmethod
    def get_filenum(filename,path,source_name,ext):
        num = int(filename[len(''.join([path,source_name])):-len(ext)])
        return num
        

    @staticmethod
    def get_filelist(path,source_name,first_frame,last_frame,ext):
        if path[-1] != '/':
            path += '/'

        if ext[0] != '.':
            ext = ''.join(['.',ext])
            
        base = ''.join([path,source_name,'*',ext])

        # Filter if file number between first and last frame
        filelist = []
        for f in glob.glob(base):
            try:
                filenum = DataCube.get_filenum(f,path,source_name,ext)
            except:
                continue
            if (filenum <= last_frame) and (filenum >= first_frame):
                filelist.append(f)
        filelist.sort()
        return filelist
        
    @classmethod
    def fromfilelist(cls,filelist,dext=0,hext=0,normalize=False):
        datalist = []
        for pfile in filelist:
            datalist.append(DataObject.fromfilename(pfile,dext=dext,hext=hext,normalize=normalize))

        cube = cls(datalist)
        cube.filelist = filelist
        return cube
        
        
    # Construct DataCube from list of files
    def __init__(self,objectlist):#path,source_name,first_frame,last_frame,ext):
        self.datalist = objectlist
        
    # Return DataObject
    def __getitem__(self,val):
        return self.datalist[val]

    def __len__(self):
        return len(self.datalist)
    
    def __call__(self,*args):
        return self

    def combine(self,method="median"):
        c = None

        if method == "average":
            c = np.mean([d.data for d in self],axis=0)
        elif method == "median":
            c = np.median([d.data for d in self],axis=0)
        elif method == "sum":
            c = np.sum([d.data for d in self],axis=0)

        return c

    
#############
### Main
############
def main():

    im_path = 'PATH'
    base_name = 'FILENAME'
    ext = 'EXT'
    start = '0'
    end = '10'

    print 'Reading from %s/%s*%s...' % (im_path,base_name,ext)
    SOURCE_list = DataCube.get_filelist(im_path,base_name,start,end,ext)
    SOURCE_dc = DataCube.fromfilelist(SOURCE_list,normalize=True)
    
if __name__ == "__main__":
    sys.exit(main())
