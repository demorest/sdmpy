#! /usr/bin/env python

# mime.py -- PBD 2016/04

# This class provdes MIME-parsing functionality suitable for reading
# EVLA/ALMA Binary Data Format (BDF) files, and SDM binary tables.  It
# is not meant to be a totally general MIME reader.  In particular,
# contents of binary type are not returned directed, rather an offset
# into the file is given.

import string
from collections import OrderedDict

class MIMEHeader(OrderedDict):
    # MIMEHeader is a dict with keys equal to the mime header keywords
    # and values equal to lists of entries.  Helper functions for parsing
    # or generating the MIME-format headers are collected here.

    @property
    def boundary(self):
        if self['Content-Type'][0].startswith('multipart/'):
            for v in self['Content-Type']:
                if v.startswith('boundary='):
                    return v[v.index('=')+1:]
        return None

    def addline(self,line):
        """
        Given a single line from a mime header, split it into key/val and
        add it to the dict.
        """
        idx = line.index(':')
        key = line[:idx]
        vals = map(string.strip, line[idx+1:].split(';'))
        self[key] = vals

    @staticmethod
    def _asline(key,val):
        """Convert given key and value list to MIME header line."""
        return key + ': ' + string.join(val, '; ') + '\n'

    def tostring(self,key=None):
        """
        Return contents as in MIME-header format.  If key is given, only
        the line corresponding to the requested key will be returned,
        otherwise the full header will be returned.
        """
        if key is not None:
            return self._asline(key,self[key])
        else:
            out = ''
            for k in self.keys():
                out += self._asline(k,self[k])
            return out

    def __str__(self):
        return self.tostring()

# TODO make a utils.py to hold stuff like this
import os
def basename_noext(path):
    return os.path.basename(os.path.splitext(path)[0])

class MIMEPart(object):
    """
    Class for representing one part of a MIME message.
    Has two member variable:

      hdr  = Dict of MIME header key/value pairs
      body = Body of message.  In our usage, can be a file offset in bytes
             (for binary parts), a string (for text) or a list of MIMEPart
             objects (for multipart).

    The loc property is a shortcut for the Content-Location header
    parameter.

    The type property is a shortcut for Content-Type
    """

    def __init__(self,fp,boundary=None,recurse=False,binary_size=None):
        """
        Read a MIME content part starting at the current file location.
        Return value is a MIMEPart object, which has elements:

            hdr    dict of MIME header key / value pairs

            body   string if Content-Type was 'text/xml', offset into
                   the file if 'application/octet-stream', or list of
                   other MIMEParts for a 'multipart/*'.

        If recurse is True, will read/return the contents of a multipart
        (and any multiparts found at lower levels).  Otherwise will read
        one header/body unit and pointer will be left at the start of 
        the next one (or first sub-part for multiparts).

        binary_size is a dict of sizes of binary components by type.
        For each binary part found, if Content-Location agrees with type,
        the binary data will be skipped over rather than read (this is
        for reading BDF files).  If binary_size is not given, or if an
        unknown type is found, the data must be read to determine its 
        size, however this last part is not implemented yet.
        """
        self.hdr = MIMEHeader({})
        self.body = None
        in_hdr = True
        binary_type = False
        multipart_type = False
        # Note, need to use readline() rather than iterating over file
        # because we need to recover file positions and seek ahead.
        # The "for line in file" type loop reads ahead so is not compatible
        # with this approach.
        while True:

            line = fp.readline()

            # hit EOF
            if line=='':
                return

            # Check for multipart boundary marker
            if boundary is not None:
                if in_hdr:
                    # If we are starting, ignore a 'start' marker,
                    # quit on a 'done' marker
                    if line=='--'+boundary+'\n':
                        continue
                    elif line=='--'+boundary+'--\n':
                        self.hdr = MIMEHeader({})
                        self.body = None
                        return
                else:
                    # This marks the end of a part, rewind so that the 
                    # next part can be parsed, and return results
                    if line.startswith('--' + boundary):
                        fp.seek(-len(line),1)
                        return

            if line=='\n':
                # Got blank line, the next part will be body.  We
                # want to skip it if this is a binary part, otherwise
                # read and return the body.
                in_hdr = False
                if binary_type:
                    # Note the location within the file and skip
                    # ahead by the correct amount.
                    bin_name = basename_noext(self.hdr['Content-Location'][0])
                    self.body = fp.tell()
                    # NOTE if the list of binary sizes is not given, or
                    # if the bin_name is unknown, we should read the data
                    # until the boundary marker is found to determine the 
                    # size.  THIS NEEDS TO BE IMPLEMENTED SOMETIME.  For
                    # now this case will raise an error.
                    if ((binary_size is None) 
                            or (bin_name not in binary_size.keys())):
                        raise RuntimeError("Unknown binary type '%s' found"
                                % bin_name)
                    # Need to add one extra byte for the newline
                    fp.seek(binary_size[bin_name]+1, 1)
                elif multipart_type:
                    if recurse:
                        # Parse the parts and add to a list
                        while True:
                            pmime = MIMEPart(fp, boundary=boundary,
                                        recurse=True,
                                        binary_size=binary_size)
                            if pmime.hdr == {}:
                                return
                            else:
                                self.body.append(pmime)
                continue

            if in_hdr:
                # Still in the header, parse the line as MIME key/val
                self.hdr.addline(line)
                if 'Content-Type' in line:
                    vals = self.hdr['Content-Type']
                    if vals[0].startswith('multipart/'):
                        multipart_type = True
                        boundary = self.hdr.boundary
                        self.body = []
                    elif vals[0] == 'application/octet-stream':
                        binary_type = True
            else:
                if not binary_type:
                    # In body part of a non-binary type
                    if self.body is None: self.body = line
                    else: self.body += line
                else:
                    # Should not really get here, means size calculation
                    # failed or file is otherwise messed up... what to do?
                    raise RuntimeError('MIME parsing failure')

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

