#! /usr/bin/env python

# sdm.py -- P. Demorest, 2015/03

import os.path
from lxml import etree, objectify

from .scan import Scan
from .mime import MIMEPart, MIMEHeader

_install_dir = os.path.abspath(os.path.dirname(__file__))
_xsd_dir = os.path.join(_install_dir, 'xsd')
_sdm_xsd = os.path.join(_xsd_dir, 'sdm_all.xsd')
# Might want to make the schema location settable somehow?
_sdm_parser = objectify.makeparser(schema=etree.XMLSchema(file=_sdm_xsd))

class SDM(object):
    """
    Top-level class to represent an SDM.

    Init arguments:
      path = path to SDM directory
      bdfdir = different directory to search for bdfs (optional, for pre-archive SDMs)

    Attributes:
      tables = list of tables
      path   = full path to SDM directory

    SDM['TableName'] returns the relevant SDMTable object.
    """
    def __init__(self,path='.',use_xsd=True, bdfdir=''):
        if use_xsd: parser = _sdm_parser
        else: parser = None
        self._tables = {}
        self._schemaVersion = {}
        self.path = os.path.abspath(path)
        self.bdfdir = bdfdir
        self._asdmtree = objectify.parse(path+'/ASDM.xml',parser)
        self.asdm = self._asdmtree.getroot()
        for tab in self.asdm.Table:
            self._schemaVersion[tab.Name] = tab.Entity.attrib['schemaVersion']
            # TODO, compare schema versions, relax parsing if they don't match
            self._tables[tab.Name] = sdmtable(tab.Name,path,use_xsd=use_xsd)

    @property
    def tables(self):
        """Return the list of table names"""
        return self._tables.keys()

    def __getitem__(self,key):
        return self._tables[key]

    def scan(self,idx):
        """Return a Scan object for the given scan number."""
        return Scan(self,str(idx))

    def scans(self, hasbdf=False):
        """Iterate over scans.  Set hasbdf=True to only return scans for which BDFs exist."""
        # List of SDM scan numbers:
        if hasbdf:
            scanidx = [s.scanNumber for s in self['Scan']
                       if self.scan(s.scanNumber).bdf.exists]
        else:
            scanidx = [s.scanNumber for s in self['Scan']]

        for idx in scanidx:
            yield self.scan(idx)

    def _update_ASDM(self):
        """Updates the ASDM table with the current number of rows."""
        # TODO could check UIDs as well
        for tab in self.asdm.Table:
            try:
                nrow = len(self[tab.Name])
                tab.NumberRows = nrow
            except TypeError:
                # This is a workaround since len is not yet implemented for
                # binary tables.
                pass

    def write(self,newpath):
        """Write the SDM out to the new path location.  Currently does not
        copy the BDFs or anything else under ASDMBinary."""
        self._update_ASDM()
        if not os.path.exists(newpath): os.mkdir(newpath)
        # Write ASDM.xml
        objectify.deannotate(self._asdmtree,cleanup_namespaces=True)
        self._asdmtree.write(newpath+'/ASDM.xml',
                encoding='utf-8', pretty_print=True, standalone=True)
        # Call each table's write method for the rest
        for tab in self.tables: self[tab].write(newpath)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

def sdmtable(name,path,*args,**kwargs):
    """
    Return the correct type of SDM table object (binary or XML).
    """
    fnamebase = os.path.join(path,str(name))
    if os.path.exists(fnamebase + '.xml'):
        return SDMTable(name,path,*args,**kwargs)
    elif os.path.exists(fnamebase + '.bin'):
        return SDMBinaryTable(name,path,*args,**kwargs)
    return None

def decap(s):
    return s[:1].lower() + s[1:] if s else ''

class SDMTable(object):
    """
    Class for an individual SDM table.

    Generally this should not be used directly, but as part of a full
    SDM via the SDM class.  

    Init arguments:
      name = Name of table (not including .xml extension)
      path = Path to SDM directory
      idtag = Name of ID tag in table (defaults to nameId)

    SDMTable[i] returns the i-th row as an lxml objectify object
    """

    # Any non-standard Id tag names can be listed here.
    # Note that for some tables a unique key is not a single tag
    # but is defined as a combination of several tags.  We are
    # just going to ignore this for now and do something that is
    # convenient and seems to work for most EVLA data sets.
    _idtags = {
            'Main': 'scanNumber',
            'Scan': 'scanNumber',
            'Subscan': 'scanNumber',
            }

    def __init__(self,name,path='.',idtag=None,use_xsd=True):
        self.name = name
        if idtag is None:
            try:
                self.idtag = self._idtags[name]
            except KeyError:
                self.idtag = decap(str(name)) + 'Id'
        else:
            self.idtag = idtag
        if use_xsd: parser = _sdm_parser
        else: parser = None
        self._tree = objectify.parse(path+'/'+name+'.xml',parser)
        self._table = self._tree.getroot()

    @property
    def entityId(self):
        """Shortcut to entityId"""
        return self._table.Entity.get('entityId')

    @property
    def containerId(self):
        """Shortcut to ContainerEntity entityId"""
        return self._table.ContainerEntity.get('entityId')

    def __getitem__(self,key):
        if type(key)==int:
            return self._table.row[key]
        else:
            # Search through the table and find the first with
            # the matching id tag:
            for r in self._table.row:
                try:
                    if str(getattr(r,self.idtag)) == key:
                        return r
                except AttributeError:
                    pass
            # No matching rows:
            raise KeyError(key)

    def __len__(self):
        try:
            return len(self._table.row)
        except AttributeError:
            return 0

    def write(self,newpath,fname=None):
        """
        Write the updated XML file to the specified path.  Will be named 
        TableName.xml unless overridden via the fname argument.
        """
        objectify.deannotate(self._tree, cleanup_namespaces=True)
        if fname is None:
            outf = os.path.join(newpath,self.name+'.xml')
        else:
            outf = os.path.join(newpath,fname)
        self._tree.write(outf, encoding='utf-8', pretty_print=True,
                standalone=True)

#class SDMTableRow(objectify.ObjectifiedElement):
#    """
#    In case we want to add any extra functionality to table rows, 
#    could try to get this working.  Will require some lxml-fu, and is
#    not currently implemented.
#    """
#
#    def __str__(self):
#        return str(self.__dict__)
#
#    @property
#    def keys(self):
#        return self.__dict__.keys()

class SDMBinaryTable(object):
    """
    Represents an SDM binary table.  Not really implemented yet, but will
    read the data and write it back out when asked to do so by the main
    SDM class.
    """
    # Notes: 
    #
    # Binary tables use MIME multipart format.  Should only have
    # two parts.  First has content-id "<header.xml>" and is an XML
    # description of the table (column names only).  Second has content-id
    # "<content.bin>" and is the table in binary format.
    # 
    # All binary numbers are big endian.
    #
    # In the binary table, strings are encoded as an int giving string length
    # followed by the string data.  Same with vectors of floats, etc.
    #
    # Entities are sets of 5 strings: entityId, entityIdEncrypted, 
    # entityTypeName, schemaVersion, documentVersion
    #
    # The binary table has the following contents:
    #  - table entity
    #  - container entity
    #  - number of rows (int)
    #  - N_row times row data
    #
    # Optional entries within a row are preceded by a 1-byte boolean
    # that is 1 if the entry exists, 0 if not.
    #
    # I have not yet found a way of automatically determining what 
    # types of data are contained within the row.  This may not be 
    # actually be possible!
    #
    def __init__(self,name,path,use_xsd=None):
        self.name = name
        fp = open(path+'/'+name+'.bin','r')
        self._data = fp.read()
        fp.seek(0)
        try:
            mimetmp = MIMEPart(fp,recurse=True)
            ## Assume part 0 is header; TODO do more checks
            self.header = objectify.fromstring(mimetmp.body[0].body)
            self._doffs = mimetmp.body[1].body
            self._dsize = mimetmp.body[1].size
        except RuntimeError:
            # Don't die on a truncated file.  
            # Probably some better way to handle this..
            pass
        fp.close()

    def write(self,newpath,fname=None):
        if fname is None:
            outf = os.path.join(newpath,self.name+'.bin')
        else:
            outf = os.path.join(newpath,fname)
        open(outf,'w').write(self._data)

    def get_bytes(self,offs,nbytes):
        return self._data[self._doffs+offs:self._doffs+offs+nbytes]

