#! /usr/bin/env python

try:
    from lxml import etree
except ImportError:
    from xml.etree import ElementTree as etree

class SDM(object):
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
        return self._tables.keys()

    def __getitem__(self,key):
        return self._tables[key]


def decap(s):
    return s[:1].lower() + s[1:] if s else ''

class SDMTable(object):
    def __init__(self,name,path='.',idtag=None):
        self.name = name
        if idtag is None:
            self.idtag = decap(name) + 'Id'
        else:
            self.idtag = idtag
        table = etree.parse(path + '/' + name + '.xml').getroot()
        self.rows = []
        self._ids = {} # Maybe ordereddict?
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
    def __init__(self,element):
        for c in element.getchildren():
            setattr(self,c.tag,c.text)

    def __str__(self):
        return str(self.__dict__)

    @property
    def keys(self):
        return self.__dict__.keys()
