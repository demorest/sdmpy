#! /usr/bin/env python

# bdf.py -- PBD 2015/03/24

# This class provides access to EVLA/ALMA Binary Data Format (BDF)
# files.  The approach used here is based somewhat on 'bdfparse.py' 
# originally by Peter Williams.

try:
    from lxml import etree
except ImportError:
    from xml.etree import ElementTree as etree

import sys
import string
import re
import mmap
import numpy
from collections import namedtuple

class MIMEPart(namedtuple('MIMEPart','hdr body')):
    """Simple class for representing one part of a MIME message.
    Has two member variable:

      hdr  = Dict of MIME header key/value pairs
      body = Body of message.  In our usage, can be a file offset in bytes
             (for binary parts), a string (for text) or a list of MIMEPart
             objects (for multipart).

    The loc property is a shortcut for the Content-Location header
    parameter.

    The type property is a shortcut for Content-Type
    """

    @property
    def loc(self):
        try:
            return self.hdr['Content-Location'][0]
        except KeyError:
            return None

    @property
    def type(self):
        try:
            return self.hdr['Content-Type'][0]
        except KeyError:
            return None

class BDF(object):

    def __init__(self, fname):
        self.fname = fname
        self.fp = open(fname, 'r')
        self.mmdata = mmap.mmap(self.fp.fileno(), 0, mmap.MAP_PRIVATE,
                mmap.PROT_READ)
        # TODO: more stuff

    @staticmethod
    def split_mime(line):
        idx = line.index(':')
        key = line[:idx]
        vals = map(string.strip, line[idx+1:].split(';'))
        return (key, vals)

    @staticmethod
    def mime_boundary(mime_hdr):
        if mime_hdr['Content-Type'][0].startswith('multipart/'):
            for v in mime_hdr['Content-Type']:
                if v.startswith('boundary='):
                    return v[v.index('=')+1:]
        return None

    def read_mime_part(self,boundary=None,recurse=False):
        """
        Read a MIME content part starting at the current file location.
        Return value is (hdr, body) tuple:

            hdr    dict of MIME header key / value pairs

            body   string if Content-Type was 'text/xml', offset into
                   the file if 'application/octet-stream', or None
                   for a 'multipart/*'.

        File pointer will be left at the start of the next part.
        TODO Do we really want to do this recursive for multiparts?
        """
        hdr = {}
        body = None
        in_hdr = True
        binary_type = False
        multipart_type = False
        # Note, need to use readline() rather than iterating over file
        # because we need to recover file positions and seek ahead.
        # The "for line in file" type loop reads ahead so is not compatible
        # with this approach.
        while True:

            line = self.fp.readline()

            # hit EOF
            if line=='':
                return MIMEPart(hdr, body)

            # Check for multipart boundary marker
            if boundary is not None:
                if in_hdr:
                    # If we are starting, ignore a 'start' marker,
                    # quit on a 'done' marker
                    if line=='--'+boundary+'\n':
                        continue
                    elif line=='--'+boundary+'--\n':
                        return MIMEPart({}, None)
                else:
                    # This marks the end of a part, rewind so that the 
                    # next part can be parsed, and return results
                    if line.startswith('--' + boundary):
                        self.fp.seek(-len(line),1)
                        return MIMEPart(hdr, body)

            if line=='\n':
                # Got blank line, the next part will be body.  We
                # want to skip it if this is a binary part, otherwise
                # read and return the body.
                in_hdr = False
                if binary_type:
                    # Note the location within the file and skip
                    # ahead by the correct amount.
                    bin_name = hdr['Content-Location'][0].split('/')[-1]
                    bin_name = re.sub('\.bin$','',bin_name)
                    body = self.fp.tell()
                    # Need to add one extra byte for the newline
                    self.fp.seek(self.bin_size[bin_name] 
                            * self.bin_dtype_size[bin_name] + 1, 1)
                elif multipart_type:
                    if recurse:
                        # Parse the parts and add to a list
                        while True:
                            #print "recur b='%s'" % boundary 
                            pmime = self.read_mime_part(boundary=boundary,
                                        recurse=True)
                            #print phdr
                            #print pbody
                            if pmime.hdr == {}:
                                return MIMEPart(hdr,body)
                            else:
                                body.append(pmime)
                continue

            if in_hdr:
                # Still in the header, parse the line as MIME key/val
                (key, vals) = self.split_mime(line)
                hdr[key] = vals
                if key=='Content-Type':
                    if vals[0].startswith('multipart/'):
                        multipart_type = True
                        boundary = self.mime_boundary(hdr)
                        body = []
                    elif vals[0] == 'application/octet-stream':
                        binary_type = True
            else:
                if not binary_type:
                    # In body part of a non-binary type
                    if body is None: body = line
                    else: body += line
                else:
                    # Should not really get here, means size calculation
                    # failed or file is otherwise messed up... what to do?
                    raise RuntimeError('BDF MIME parsing failure')

    @staticmethod
    def stripns(tag):
        return re.sub('{.+}','',tag)

    # Size in bytes of each data element. In principle, the crossData 
    # type needs to be read from the headers while all others have 
    # pre-set values...
    bin_dtype_size = {
            'flags':           4, # INT32
            'actualTimes':     8, # INT64
            'actualDurations': 8, # INT64
            'zeroLags':        4, # FLOAT32
            'autoData':        4, # FLOAT32
            'crossData':       4, # FLOAT32 but should be determined from file
            }

    def read_top(self):
        self.fp.seek(0,0) # Go back to start
        if not self.fp.readline().startswith('MIME-Version:'):
            raise RuntimeError('Invalid BDF: missing MIME-Version')
        # Read top-level MIME; do we need to save any of this stuff?
        mime_hdr = self.read_mime_part().hdr
        # Read main xml header
        sdmDataMime = self.read_mime_part(boundary=self.mime_boundary(mime_hdr))
        if sdmDataMime.loc != 'sdmDataHeader.xml':
            raise RuntimeError('Invalid BDF: missing sdmDataHeader.xml')
        self.sdmDataHeader = etree.fromstring(sdmDataMime.body)
        # Find the sizes of all binary parts
        self.bin_size = {}
        for e in self.sdmDataHeader.iter():
            if 'size' in e.attrib.keys() and 'axes' in e.attrib.keys():
                self.bin_size[self.stripns(e.tag)] = int(e.attrib['size'])
        self.fp.seek(0,0) # reset to start again

