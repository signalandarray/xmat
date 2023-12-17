import numpy as np
import io
from functools import reduce  # for numel()
import operator


DTYPE_INFO_MAP = {}


def _add_type(typename_numpy, typename_format=None):
    if not typename_format:
        typename_format = typename_numpy
    DTYPE_INFO_MAP[typename_format] = np.dtype(typename_numpy)


# https://numpy.org/devdocs/user/basics.types.html
_add_type('S', 'char')
_add_type('int8')
_add_type('uint8')
_add_type('int32')
_add_type('uint32')
_add_type('int64')
_add_type('uint64')
_add_type('float32')
_add_type('float64')
_add_type('complex64')
_add_type('complex128')


def print_type_info():
    for name, tb in DTYPE_INFO_MAP.items():
        print(f'typename: {name: >12}, char: {tb.char}, itemsize: {tb.itemsize},'
              f'isaligned: {tb.isalignedstruct}, descr: {tb.descr[0][1]}')


class Header:

    FORMAT_SIGNATURE_SIZE = 4
    FORMAT_SIGNATURE = 'XYZ'
    FORMAT_FOOTER = 'end'
    UFIX_SIZE = 4
    MAX_BLOCK_NAME_LEN = 64
    MAX_TYPE_NAME_LEN = 32
    MAX_NDIM = 8

    def __init__(self,
                 signature=FORMAT_SIGNATURE,
                 ufix_name=UFIX_SIZE,
                 max_block_name_len=MAX_BLOCK_NAME_LEN,
                 max_type_name_len=MAX_TYPE_NAME_LEN,
                 max_ndim=MAX_NDIM):
        self.format_signature: str = signature
        self.ufix_size: int = ufix_name
        self.max_block_name_length: int = max_block_name_len
        self.max_type_name_length: int = max_type_name_len
        self.max_ndim: int = max_ndim

    @staticmethod
    def read(fid, order='little'):
        signature = fid.read(Header.FORMAT_SIGNATURE_SIZE).decode('ascii').split('\0')[0]
        if signature != Header.FORMAT_SIGNATURE:
            raise RuntimeError("wrong format signature header begin")
        ufix_size = int.from_bytes(fid.read(1), order)
        maxblocknamelen = int.from_bytes(fid.read(ufix_size), order)
        maxtypenamelen = int.from_bytes(fid.read(ufix_size), order)
        maxndim = int.from_bytes(fid.read(ufix_size), order)
        signarute_end = fid.read(Header.FORMAT_SIGNATURE_SIZE).decode('ascii').split('\0')[0]
        if signarute_end != signature:
            raise RuntimeError("wrong format signature header end")
        return Header(signature, ufix_size,
                      maxblocknamelen, maxtypenamelen, maxndim)

    def write(self, os):
        if os.tell():
            raise RuntimeError("file cursor isn't in the beginning if the file")
        _write_string(os, self.format_signature.ljust(Header.FORMAT_SIGNATURE_SIZE, '\0'))
        _write_ufix(os, self.ufix_size, 1)
        _write_ufix(os, self.max_block_name_length, self.ufix_size)
        _write_ufix(os, self.max_type_name_length, self.ufix_size)
        _write_ufix(os, self.max_ndim, self.ufix_size)
        _write_string(os, self.format_signature.ljust(Header.FORMAT_SIGNATURE_SIZE, '\0'))

    def byte_size(self):
        return 2*Header.FORMAT_SIGNATURE_SIZE + 1 + 3 * self.ufix_size

    def __str__(self):
        os = io.StringIO()
        for key, val in self.__dict__.items():
            print(key, ": ", val, file=os, end=', ')
        content = os.getvalue()
        os.close()
        return content


class BlockDescriptor:

    SIGNATURE_BEGIN = '<#>'
    SIGNATURE_END = '>#<'

    def __init__(self, name='', type='', typesize=0, ndim=0, shape=(), pos=-1):
        self.name = name
        self.type = type
        self.typesize = typesize
        self.ndim = ndim
        self.shape = shape
        self.numel = calc_numel(self.shape)
        self.pos = pos   # for support

    @staticmethod
    def make(obj):
        if isinstance(obj, str):
            pass
        elif isinstance(obj, np.ndarray):
            pass
        else:
            raise RuntimeError("wrong type")

    @staticmethod
    def read(fid, h: Header, order='little'):
        signature_begin = _read_string(fid, Header.FORMAT_SIGNATURE_SIZE)
        if signature_begin != BlockDescriptor.SIGNATURE_BEGIN:
            raise RuntimeError("wrong blockDescriptor signature begin")
        name = _read_string(fid, h.max_block_name_length)
        type = _read_string(fid, h.max_type_name_length)
        typesize = _read_int(fid, h.ufix_size)
        ndim = _read_int(fid, h.ufix_size)
        shape = [_read_int(fid, h.ufix_size) for _ in range(h.max_ndim)]
        shape = shape[:ndim]
        numel = _read_int(fid, h.ufix_size)

        signature_end = _read_string(fid, Header.FORMAT_SIGNATURE_SIZE)
        if signature_end != BlockDescriptor.SIGNATURE_END:
            raise RuntimeError("wrong blockDescriptor signature end")
        return BlockDescriptor(name, type, typesize, ndim, shape, numel)

    def write(self, os, h: Header):
        _write_string(os, self.SIGNATURE_BEGIN.ljust(Header.FORMAT_SIGNATURE_SIZE, '\0'))
        _write_string(os, self.name.ljust(h.max_block_name_length, '\0'))
        _write_string(os, self.type.ljust(h.max_type_name_length, '\0'))
        _write_ufix(os, self.typesize, h.ufix_size)
        _write_ufix(os, self.ndim, h.ufix_size)
        for i in self.shape:
            _write_ufix(os, i, h.ufix_size)
        for i in range(h.max_ndim - self.ndim):
            _write_ufix(os, 0, h.ufix_size)
        _write_ufix(os, self.numel, h.ufix_size)
        _write_string(os, self.SIGNATURE_END.ljust(Header.FORMAT_SIGNATURE_SIZE, '\0'))

    def buffersize(self):
        return self.numel*self.typesize

    def __str__(self):
        os = io.StringIO()
        for key, val in self.__dict__.items():
            print(key, ": ", val, file=os, end=', ')
        content = os.getvalue()
        os.close()
        return content

    def __repr__(self):
        return self.__str__()


# simple read/write functions
# ---------------------------
def _read_string(fid, length):
    return fid.read(length).split(b'\0')[0].decode('ascii')


def _write_string(fid, s: str):
    fid.write(s.encode('ascii'))


def _read_int(fid, len, order='little'):
    return int.from_bytes(fid.read(len), order)


def _write_ufix(fos, num: int, length, order='little'):
    fos.write((num).to_bytes(length, order))


def _read_block_buffer(buffer, block):
    dt = DTYPE_INFO_MAP[block.type]
    if dt.itemsize != block.typesize:
        raise RuntimeError('wrong itemsize')
    x = np.frombuffer(buffer, dt, block.numel)
    x = x.reshape(block.shape)
    return x


def calc_numel(shape):
    return reduce(operator.mul, shape, 1)


# User Classes
# ------------
class Writer:
    def __init__(self, filename, h: Header=None):
        self.os = open(filename, 'wb')
        if not h:
            self.h = Header()
        self.h.write(self.os)

        self.block = None
        self.content = {}

    def close(self):
        self.os.close()

    def save(self, name: str, obj):
        if isinstance(obj, str):
            dtname = 'int8'
            dt = DTYPE_INFO_MAP[dtname]
            obj += '\0'
            bd = BlockDescriptor(name, dtname, dt.itemsize, ndim=1, shape=[len(obj)])
            buffer = obj.encode('ascii')
        elif isinstance(obj, np.ndarray):
            dtname = obj.dtype.name
            if dtname not in DTYPE_INFO_MAP:
                raise RuntimeError(f'unsupported dtype')
            dt = DTYPE_INFO_MAP[dtname]
            bd = BlockDescriptor(name, dtname, dt.itemsize, ndim=obj.ndim, shape=obj.shape)
            buffer = obj.tobytes('C')
        else:
            raise RuntimeError("wrong type. only str, nympy.ndarray are supported")

        if bd.buffersize() != len(buffer):
            raise RuntimeError("wrong buffersize")
        bd.pos = self.os.tell()
        bd.write(self.os, self.h)
        n_ = self.os.write(buffer)
        if n_ != bd.buffersize():
            raise RuntimeError('not the all data has been written')
        self.block = bd
        return self.os.tell()

    def _write_footer(self):
        pass


class Reader:
    def __init__(self, filename):
        self.is_ = open(filename, 'rb')
        self.filelen = self.is_.seek(0, 2)
        self.is_.seek(0, 0)

        self.h: Header = Header.read(self.is_)
        self.block = None
        self.content = {}
        self._init_content_map()

    def close(self):
        self.is_.close()

    def reset(self):
        self.is_.seek(self.h.byte_size())

    def load(self, fieldname):
        if fieldname not in self.content:
            raise RuntimeError(f'wrong fieldname')
        block = self.content[fieldname]
        self.is_.seek(block.pos)
        self._read_descriptor()
        buffer = self.is_.read(self.block.buffersize())
        obj = _read_block_buffer(buffer, self.block)
        return obj

    def next(self):
        self._read_descriptor()
        if not self.block:
            return self.block
        self.is_.seek(self.block.buffersize(), 1)
        return self.block

    def _read_descriptor(self):
        signature = _read_string(self.is_, Header.FORMAT_SIGNATURE_SIZE)
        if signature != BlockDescriptor.SIGNATURE_BEGIN:
            # if signature != Header.FORMAT_FOOTER:
            #     raise RuntimeError(f'wrong end-signature: {signature}')
            self.block = None
            return self.block
        pos = self.is_.seek(-Header.FORMAT_SIGNATURE_SIZE, 1)
        self.block = BlockDescriptor.read(self.is_, self.h)
        self.block.pos = pos

    def _init_content_map(self):
        self.reset()
        while self.next():
            self.content[self.block.name] = self.block
        self.reset()

    def __repr__(self):
        os = io.StringIO()
        for block in self.content.values():
            print(block, file=os)
        content = os.getvalue()
        os.close()
        return content
