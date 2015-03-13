#! /usr/bin/env python

# sdm.py -- P. Demorest, 2015/03

try:
    from lxml import etree
except ImportError:
    from xml.etree import ElementTree as etree

class SDM(object):
    """
    Top-level class to represent an SDM.

    Init arguments:
      path = path to SDM directory

    Attributes:
      tables = list of tables with non-zero number of rows

    SDM['TableName'] returns the relevant SDMTable object.
    """
    def __init__(self,path='.'):
        self._tables = {}
        asdm = etree.parse(path+'/ASDM.xml').getroot()
        for tab in asdm.iter('Table'):
            name = tab.find('Name').text
            nrows = int(tab.find('NumberRows').text)
            if nrows>0:
                try:
                    self._tables[name] = SDMTable(name,path)
                except IOError:
                    pass

    @property
    def tables(self):
        """Return the list of non-empty table names"""
        return self._tables.keys()

    def __getitem__(self,key):
        return self._tables[key]


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

    SDMTable[i] returns the i-th row as a SDMTableRow object
    """
    def __init__(self,name,path='.',idtag=None):
        self.name = name
        if idtag is None:
            self.idtag = decap(name) + 'Id'
        else:
            self.idtag = idtag
        table = etree.parse(path + '/' + name + '.xml').getroot()
        entity = table.find('Entity')
        self.entityId = entity.attrib['entityId']
        self.rows = []
        self._ids = {}
        for row in table.iter('row'):
            newrow = SDMTableRow(row)
            self.rows.append(newrow)
            # TODO : ids can map to more than one row..
            if hasattr(newrow,self.idtag):
                self._ids[getattr(newrow,self.idtag)] = len(self.rows)-1

    def __getitem__(self,key):
        if type(key)==int:
            return self.rows[key]
        else:
            return self.rows[self._ids[key]]

    def __len__(self):
        return len(self.rows)

    @property
    def ids(self):
        return self._ids.keys()

class SDMTableRow(object):
    """
    Represents an individual row in an SDM Table.

    Generally this should not be used directly, but as part of a full
    SDM via the SDM and SDMTable classes.  

    Init arguments:
      element = etree XML element for row

    Attributes:
      keys = list of field names in row

    Values accessed as SDMTableRow.keyname.  All values are kept as strings
    for now.
    """
    def __init__(self,element):
        for c in element.getchildren():
            ent_ref = c.find('EntityRef')
            if ent_ref is None:
                setattr(self,c.tag,c.text)
            else:
                # Could read more of the entity reference stuff..
                setattr(self,c.tag,ent_ref.attrib['entityId'])

    def __str__(self):
        return str(self.__dict__)

    @property
    def keys(self):
        return self.__dict__.keys()

