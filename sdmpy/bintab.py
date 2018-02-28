from __future__ import print_function, division, absolute_import, unicode_literals # not casa compatible
from builtins import bytes, dict, object, range, map, input#, str # not casa compatible
from future.utils import itervalues, viewitems, iteritems, listvalues, listitems
from io import open

import struct
import numpy

# Contains code for unpacking data from SDM binary tables.


def unpacker(sdmtable):
    """Return an appropriate unpacker given an SDMBinaryTable object."""
    # Could do something more fancy with a registry, etc but it's probably
    # not worth it given the small number of binary tables.
    if sdmtable.name == 'SysPower':
        return SysPowerUnpacker(sdmtable)
    else:
        raise RuntimeError("Unknown binary table type: '%s'" % sdmtable.name)


class BinaryTableUnpacker(object):
    """Base class for SDM binary table unpackers.  This should not be called
    directly, because it needs a column data structure definition in order
    to be useful.  This is implemented in derived classes, e.g.
    SysPowerUnpacker for the SysPower table."""

    # Unpackers for various binary data types
    _unpack_int = struct.Struct('>i')
    _unpack_long = struct.Struct('>q')
    _unpack_float = struct.Struct('>f')

    # Column definitions go in derived classes, see SysPowerUnpacker
    # below for an example.
    columns = []

    def __init__(self, sdmtable):
        self._sdmtable = sdmtable  # The original SDMBinaryTable object
        self._pos = 0  # Pointer to current position within the table
        # Read the entity stuff at the start of the table, then
        # record the start of row data.  Just use tuples for the
        # entity things now, until building a better entity data
        # structure seems worthwhile.
        self.table_entity = tuple(self._get_val('S') for i in range(5))
        self.container_entity = tuple(self._get_val('S') for i in range(5))
        self.nrows = self._get_val('i4')  # -1 in all data I've looked at...
        self._pos0 = self._pos
        self.row = []

    def _readtab(self, nbytes):
        # Read/return nbytes from the table, and increment the
        # position marker.
        result = self._sdmtable.get_bytes(self._pos, nbytes)
        if len(result) != nbytes:
            return None
        self._pos += nbytes
        return result

    def _get_val(self, dtype, optional=False, array=False):
        # Unpack one value of the specified type from the table,
        # and advance the position pointer appropriately.  dtype
        # uses values compatible with numpy, ie i4, f4, etc.
        string_val = False
        if dtype == 'i4':
            unpack = self._unpack_int
        elif dtype == 'i8':
            unpack = self._unpack_long
        elif dtype == 'f4':
            unpack = self._unpack_float
        elif dtype[0] == 'S':
            unpack = None
            string_val = True
        else:
            raise RuntimeError("Unknown data dtype (dtype='%s')" % dtype)
        try:
            if optional:
                has_data = self._readtab(1) != '\0'
                if not has_data:
                    # Return some default values..
                    if string_val:
                        return ''
                    else:
                        return 0
            if array or string_val:
                nval = self._unpack_int.unpack_from(self._readtab(4))[0]
            if string_val:
                return self._readtab(nval)
            dsize = unpack.size
            if array:
                return [unpack.unpack_from(self._readtab(dsize))[0]
                        for i in range(nval)]
            else:
                return unpack.unpack_from(self._readtab(dsize))[0]
        except struct.error:
            # Not enough bytes to unpack probably means we are at the end of
            # the table.
            return None

    @property
    def record_dtype(self):
        return [(col[0], col[2], col[3]) for col in self.columns]

    def _blank_row(self, nrows=1):
        return numpy.zeros(nrows, dtype=self.record_dtype)

    def _unpack_row(self):
        # Read the next row starting at current position and fill
        # into a single-row record array.
        row = self._blank_row()
        for col in self.columns:
            # Maybe use a better data struct for these...
            colname = col[0]
            isoptional = col[1]
            coldtype = col[2]
            isarray = col[3] != ()
            val = self._get_val(coldtype, isoptional, isarray)
            if val is None:
                return None
            if isarray:
                row[0][colname][:len(val)] = val
            else:
                row[0][colname] = val
        return row

    def iterrows(self):
        # Iterator over rows.  May or may not make sense to expose this
        # since everything is probably faster on the full unpacked array..
        self._pos = self._pos0  # reset to start
        more_data = True
        while more_data:
            row = self._unpack_row()
            if row is None:
                more_data = False
            else:
                yield row

    def unpack(self):
        self.row = numpy.recarray(1000, dtype=self.record_dtype)
        irow = 0
        for row in self.iterrows():
            if irow == len(self.row):
                self.row = numpy.rec.array(
                        numpy.resize(self.row, int(len(self.row)*1.2))
                        )
            self.row[irow] = row
            irow += 1
        self.row = numpy.rec.array(numpy.resize(self.row, irow))


class SysPowerUnpacker(BinaryTableUnpacker):
    """Unpacker for the SysPower binary table.  Currently makes some
    assumptions like there never being more than 2 measurements in
    the array-valued columns.  This is true for now for the VLA.
    """

    # Define columns using tuples like:
    #   (colname, optional, dtype, array_dims)
    # The column names do not have to exactly match the SDM definition
    # but seems like a reasonable idea to stay close.  The data types
    # should be recognized by numpy.recarray dtype argument.  This means
    # we need to specify an array dimension here.  array_dims should be set
    # to () for scalar columns.

    columns = [
            ("antennaId",        False, 'S32', ()),
            ("spectralWindowId", False, 'S32', ()),
            ("feedId",           False, 'i4', ()),
            ("timeMid",          False, 'i8', ()),
            ("interval",         False, 'i8', ()),
            ("numReceptor",      False, 'i4', ()),
            ("switchedPowerDifference", True, 'f4', (2,)),
            ("switchedPowerSum",        True, 'f4', (2,)),
            ("requantizerGain",         True, 'f4', (2,)),
            ]
