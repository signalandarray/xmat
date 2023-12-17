#pragma once

#include <cstdint>
#include <cstring>
#include <cassert>

#include <iostream>
#include <fstream>
#include <array>
#include <complex>
#include <exception>
#include <algorithm>

#include <vector>
#include <string>
#include <map>

#include "types.hpp"


namespace cml {
// constants
// ---------
using ufix = std::uint32_t;

const char kUfixSize = sizeof(ufix);
const ufix kMaxBlockNameLength = 64;
const ufix kMaxTypeNameLength = 32;
const ufix kMaxNDim = 8;

const ufix kFormatSignatureSize = 4;
const char kFormatSignature[kFormatSignatureSize] = "XYZ";
const char kDescriptorBeginSignature[kFormatSignatureSize] = "<#>";
const char kDescriptorEndSignature[kFormatSignatureSize] = ">#<";
const char kFormatFooterSignature[kFormatSignatureSize] = "end";

const ufix kSpaceSizeHeader = 2*kFormatSignatureSize + 1 + 3 * kUfixSize;
const ufix kSpaceSizeDescriptor = 
  kMaxBlockNameLength + kMaxTypeNameLength + 
  (1 + 1 + kMaxNDim + 1) * kUfixSize +
  2 * kFormatSignatureSize;
const ufix kPosSize = kSpaceSizeDescriptor - kFormatSignatureSize;


// supported types registration
// ----------------------------
template<typename T>
struct TypeInfo {
  static_assert(alwaysFalse<T>(), "using not registered data_type");
  static const bool registered = false;
  static const ufix size = 0;
  static constexpr char name[kMaxTypeNameLength] = "";
};

#define CML_ADDTYPE(T, Name) \
template<> struct TypeInfo<T> {                            \
  static const bool registered = true;                     \
  static const ufix size = sizeof(T);                      \
  static constexpr char name[kMaxTypeNameLength] = Name;   \
};                                                         \


CML_ADDTYPE(std::int8_t, "int8");
CML_ADDTYPE(std::uint8_t, "uint8");
CML_ADDTYPE(std::int32_t, "int32");
CML_ADDTYPE(std::uint32_t, "uint32");
CML_ADDTYPE(std::int64_t, "int64");
CML_ADDTYPE(std::uint64_t, "uint64");
CML_ADDTYPE(char, "char");
CML_ADDTYPE(float, "float32");
CML_ADDTYPE(double, "float64");
CML_ADDTYPE(std::complex<float>, "complex64")
CML_ADDTYPE(std::complex<double>, "complex128")


// file header operations
// ----------------------
inline
void writeHeader(std::ostream& os) {
  os.write(kFormatSignature, sizeof(kFormatSignature));
  os.write(reinterpret_cast<const char*>(&kUfixSize), 1);
  os.write(reinterpret_cast<const char*>(&kMaxBlockNameLength), kUfixSize);
  os.write(reinterpret_cast<const char*>(&kMaxTypeNameLength), kUfixSize);
  os.write(reinterpret_cast<const char*>(&kMaxNDim), kUfixSize);
  os.write(kFormatSignature, sizeof(kFormatSignature));
}


//
inline
int readHeader(std::istream& is) {
  char signature[kFormatSignatureSize];
  is.read(signature, sizeof(kFormatSignatureSize));
  if (std::strcmp(signature, kFormatSignature)) throw std::runtime_error("wrong beg-signature");

  char uintsize = -1;
  is.read(reinterpret_cast<char*>(&uintsize), 1);
  if (uintsize != kUfixSize) throw std::runtime_error("wrong size of int used");

  ufix maxnamelen = -1;
  is.read(reinterpret_cast<char*>(&maxnamelen), kUfixSize);
  if (maxnamelen != kMaxBlockNameLength) throw std::runtime_error("wrong max name len");

  ufix maxtypelen = -1;
  is.read(reinterpret_cast<char*>(&maxtypelen), kUfixSize);
  if (maxtypelen != kMaxTypeNameLength) throw std::runtime_error("wrong max type len");

  ufix maxndim = -1;
  is.read(reinterpret_cast<char*>(&maxndim), kUfixSize);
  if (maxndim != kMaxNDim) throw std::runtime_error("wrong max ndim");

  is.read(signature, sizeof(kFormatSignatureSize));
  if (std::strcmp(signature, kFormatSignature)) throw std::runtime_error("wrong end-signature");
  return 0;
}

// BlockDescriptor operations
// --------------------------
struct BlockDescriptor {
  // errors
  enum class State { undef, ok, err_0, err_1 };

  char name[kMaxBlockNameLength] = "undef";
  char type[kMaxTypeNameLength] = "undef";
  ufix typesize = 0;    // size if a single element in bytes
  ufix ndim = 0;
  ufix shape[kMaxNDim] = {0};
  ufix numel = 0;
  ufix blocksize() const { return numel * typesize; }  // size of the whole block in bytes

  template<typename T>
  void initType() {
    strcpy_s(type, sizeof(type), TypeInfo<T>::name);
    typesize = TypeInfo<T>::size;
  }

  template<typename T>
  void checkType() const {
    if (!TypeInfo<T>::registered) throw std::runtime_error("type isn't registrated");
    if (std::strcmp(type, TypeInfo<T>::name)) throw std::runtime_error("wrong type.name");
    if (typesize != TypeInfo<T>::size) throw std::runtime_error("wrong type.size");
    return;
  }

  // support attributes
  State state = State::undef;
  uint pos = 0;
};

//
template<typename T1, typename T2>
T1 numel(const T1* shape, T2 ndim) {
  assert(ndim);
  if (!ndim) return 0;
  T1 num = shape[0];
  for (T2 i = 1; i < ndim; ++i) num *= shape[i];
  return num;
}


//
inline
void writeBlockDescriptor(std::ostream& os, const BlockDescriptor& bd) {
  os.write(kDescriptorBeginSignature, kFormatSignatureSize);
  os.write(bd.name, kMaxBlockNameLength);
  os.write(bd.type, kMaxTypeNameLength);
  os.write(reinterpret_cast<const char*>(&bd.typesize), kUfixSize);
  os.write(reinterpret_cast<const char*>(&bd.ndim), kUfixSize);
  os.write(reinterpret_cast<const char*>(bd.shape), kUfixSize*kMaxNDim);
  os.write(reinterpret_cast<const char*>(&bd.numel), kUfixSize);
  os.write(kDescriptorEndSignature, kFormatSignatureSize);
}


//
inline
BlockDescriptor::State readBlockDescriptor(std::istream& is, BlockDescriptor& bd) {
  bd.pos = is.tellg();
  char descriptor_signature[kFormatSignatureSize] = "";
  is.read(descriptor_signature, kFormatSignatureSize);
  if (std::strcmp(descriptor_signature, kDescriptorBeginSignature)) {
    return BlockDescriptor::State::err_0;
  }
  is.read(bd.name, kMaxBlockNameLength);
  is.read(bd.type, kMaxTypeNameLength);
  is.read(reinterpret_cast<char*>(&bd.typesize), kUfixSize);
  is.read(reinterpret_cast<char*>(&bd.ndim), kUfixSize);
  is.read(reinterpret_cast<char*>(bd.shape), kUfixSize * kMaxNDim);
  is.read(reinterpret_cast<char*>(&bd.numel), kUfixSize);
  is.read(descriptor_signature, kFormatSignatureSize);
  if (std::strcmp(descriptor_signature, kDescriptorEndSignature)) {
    return BlockDescriptor::State::err_1;
  }
  bd.state = BlockDescriptor::State::ok;
  return BlockDescriptor::State::ok;
}


//
class WritingDescriptor {
 public:
  WritingDescriptor(const char* filename) {
    _os.open(filename, std::ios::out | std::ios::binary);
    if (!_os.is_open()) throw std::runtime_error("can't open file");
    writeHeader(_os);
  }
  
  WritingDescriptor(const std::string& filename) : WritingDescriptor(filename.c_str()) {}
  virtual ~WritingDescriptor() = default;

  // methods
  // -------
  void close() {
    _os.close();
  }
  const BlockDescriptor& block() const { return _bd; }
  
  //
  template<typename T>
  void save(const char* name, const T& x) {
    static_assert(TypeInfo<T>::registered, "type isn't registered");
    _bd.initType<T>();
    strcpy_s(_bd.name, sizeof(_bd.name), name);

    _bd.ndim = 0;
    _bd.shape[0] = 0;
    _bd.numel = 1;

    writeBlockDescriptor(_os, _bd);
    _os.write(reinterpret_cast<const char*>(&x), _bd.blocksize());
  }
 
 protected:
  std::ofstream _os;
  BlockDescriptor _bd;

 private:
   void writeFooter() {
     _os.write(kFormatFooterSignature, kFormatSignatureSize);
   }
}; // WritingDescriptor


//
class ReadingDescriptor {
 public:
  ReadingDescriptor(const char* filename) {
    _is.open(filename, std::ios::in | std::ios::binary);
    if (!_is.is_open()) throw std::runtime_error("can't open file");

    _is.seekg(0, _is.end);
    _filelen = _is.tellg();
    _is.seekg(0, _is.beg);
  }

  ReadingDescriptor(const std::string& filename) : ReadingDescriptor(filename.c_str()) {}
  virtual ~ReadingDescriptor() = default;

  // methods
  // -------
  void close() { _is.close(); }
  const BlockDescriptor& block() const { return _bd; }

  template<typename T>
  T load(const char* fieldname) {
    impl_find<T>(fieldname);
    if (_bd.numel != 1) throw std::runtime_error("wrong numel");

    T x;
    _is.read(reinterpret_cast<char*>(&x), TypeInfo<T>::size);
    return x;
  }

  void reset() { _is.seekg(kSpaceSizeHeader, _is.beg); }

  bool find(const std::string& fieldname) { return find(fieldname.c_str()); }

  //
  bool find(const char* fieldname) {
    reset();
    while (next()) {
      if (!std::strcmp(_bd.name, fieldname)) {
        _is.seekg(-int(_bd.blocksize()), _is.cur);
        return true;
      }
    }
    return false;
  }

  //
  /// \raise runtime_exception if errors
  bool validate() {
    _is.seekg(0, _is.beg);
    readHeader(_is);
    while (next()) {
      if (_bd.state != BlockDescriptor::State::ok) {
        return false;
      }
    }
    reset();
    return true;
  }


  // return map<fieldname: string, BlockDescriptor>
  std::map<std::string, BlockDescriptor> content() {
    std::map<std::string, BlockDescriptor> m;
    reset();
    while (next()) { m[_bd.name] = _bd; }
    reset();
    return m;
  }


  //
  bool next() {
    const uint toend = _filelen - _is.tellg();
    if (toend < kSpaceSizeDescriptor) return false;
    
    auto errcode = readBlockDescriptor(_is, _bd);
    if (errcode != BlockDescriptor::State::ok) return false;
    
    _is.seekg(_bd.blocksize(), _is.cur);
    return true;
  }

 protected:
  // use this int load implementations
  // find block with name, read descriptor, match descriptor with Type
  template<typename T>
  bool impl_find(const char* fieldname) {
    static_assert(TypeInfo<T>::registered, "type isn't registered"); // if T is even registered
    if (!find(fieldname)) throw std::runtime_error("wrong fieldname");
    _bd.checkType<T>();
    return true;
  }

  std::ifstream _is;
  BlockDescriptor _bd;
  bool _flag_valid = false;
  uint _filelen;
}; // ReadingDescriptor


// specific data types handling
// ----------------------------

// interface for specific data savind
class Saver : public WritingDescriptor {
 public:
  using WritingDescriptor::WritingDescriptor;

  // std containerss
  template<typename Arr>
  void save1d(const char* fieldname, const Arr& a) { save1d(fieldname, std::begin(a), std::end(a)); }

  //
  template<typename InputIt>
  void save1d(const char* fieldname, InputIt first, InputIt last) {
    using T = typename InputIt::value_type;
    static_assert(TypeInfo<T>::registered, "type isn't registered");
    _bd.initType<T>();
    strcpy_s(_bd.name, sizeof(_bd.name), fieldname);
    _bd.ndim = 1;
    _bd.shape[0] = last - first;
    _bd.numel = numel(_bd.shape, _bd.ndim);
    writeBlockDescriptor(_os, _bd);
    
    for (; first != last; ++first) {
      _os.write(reinterpret_cast<const char*>(&(*first)), TypeInfo<T>::size);
    }
  }
  
  //
  template<typename A>
  void saveNd(const char* fieldname, const A& a) {
    using T = typename A::value_type;
    static_assert(TypeInfo<T>::registered, "type isn't registered");
    _bd.initType<T>();
    strcpy_s(_bd.name, sizeof(_bd.name), fieldname);
    _bd.ndim = a.ndim();
    for (uint i = 0; i < a.ndim(); ++i) _bd.shape[i] = a.shape()[i];
    _bd.numel = numel(_bd.shape, _bd.ndim);
    writeBlockDescriptor(_os, _bd);

    for (ufix i = 0; i < _bd.numel; ++i) {
      _os.write(reinterpret_cast<const char*>(&(a[i])), TypeInfo<T>::size);
    }
  }
  
  //
  template<typename Mx>
  void saveMx(const char* fieldname, const Mx& mx) {
    using T = typename Mx::Scalar;
    static_assert(TypeInfo<T>::registered, "type isn't registered");
    _bd.initType<T>();
    strcpy_s(_bd.name, sizeof(_bd.name), fieldname);
    _bd.ndim = 2;
    _bd.shape[0] = mx.rows();
    _bd.shape[1] = mx.cols();
    _bd.numel = numel(_bd.shape, _bd.ndim);
    writeBlockDescriptor(_os, _bd);

    for (uint i0 = 0, N0 = mx.rows(); i0 < N0; ++i0) {
      for (uint i1 = 0, N1 = mx.cols(); i1 < N1; ++i1) {
        _os.write(reinterpret_cast<const char*>(&mx(i0, i1)), TypeInfo<T>::size);
      }
    }
  }
  
  // may be not used
  template<typename T, typename A, typename UInt>
  void implSaveNd(const char* fieldname, const A& a, const ufix ndim, const UInt* shape) {
    static_assert(TypeInfo<T>::registered, "type isn't registered");
    _bd.initType<T>();
    strcpy_s(_bd.name, sizeof(_bd.name), fieldname);
    _bd.ndim = ndim;
    for (uint i = 0; i < ndim; ++i) _bd.shape[i] = shape[i];
    _bd.numel = numel(_bd.shape, _bd.ndim);
    writeBlockDescriptor(_os, _bd);
    
    for (ufix i = 0; i < _bd.numel; ++i) {
      _os.write(reinterpret_cast<const char*>(&(a[i])), _bd.typesize);
    }
  }
}; // Saver -------------------


// interface for specific data loading
class Loader : public ReadingDescriptor {
public:
  using ReadingDescriptor::ReadingDescriptor;
  
  // std containers
  template<typename Arr>
  Arr load1d(const char* fieldname) {
    using T = typename Arr::value_type;
    impl_find<T>(fieldname);
    
    Arr a;
    a.resize(_bd.numel);
    load1d(fieldname, a.begin(), a.end(), true);
    return a;
  }

  //
  template<typename OutIt>
  void load1d(const char* fieldname, OutIt first, OutIt last, bool stable=false) {
    using T = typename OutIt::value_type;
    if (!stable) {
      impl_find<T>(fieldname);
    }
    ufix cursor = 0;
    for (; first != last && cursor < _bd.numel; ++first, ++cursor) {
      _is.read(reinterpret_cast<char*>(&(*first)), _bd.typesize);
    }
    _is.seekg(_bd.blocksize() - cursor, _is.cur);
  }

  // arrayNd
  template<typename A>
  A loadNd(const char* fieldname) {
    using T = typename A::value_type;
    impl_find<T>(fieldname);
    A a(_bd.shape, _bd.ndim);
    loadNd(fieldname, a, true);
    return a;
  }

  //
  template<typename A>
  void loadNd(const char* fieldname, A& a, bool stable=false) {
    using T = typename A::value_type;
    if (!stable) {
      impl_find<T>(fieldname);
    }
    assert(a.size() == _bd.numel);
    for (uint i = 0; i < _bd.numel; ++i) {
      _is.read(reinterpret_cast<char*>(&a[i]), _bd.typesize);
    }
  }

  // Eigen Matrix
  template<typename Mx>
  Mx loadMx(const char* fieldname) {
    Mx mx;
    loadMx(fieldname, mx);
    return mx;
  }

  //
  template<typename Mx>
  void loadMx(const char* fieldname, Mx& mx, bool stable=false) {
    using T = typename Mx::Scalar;
    if (!stable) {
      impl_find<T>(fieldname);
    }
    if (_bd.ndim < 1 || _bd.ndim > 2) throw std::runtime_error("wrong ndim. ndim = 2 or 1 expected");
    if (_bd.ndim == 1) {
      if (!(mx.rows() == 0 || mx.rows() == _bd.shape[0]) ||
        !(mx.cols() == 0 || mx.cols() == 1)) throw std::runtime_error("wrong mx.shape 1d");
      mx.resize(_bd.shape[0], 1);
    }
    else {
      if (!(mx.rows() == 0 || mx.rows() == _bd.shape[0]) ||
        !(mx.cols() == 0 || mx.cols() == _bd.shape[1])) throw std::runtime_error("wrong mx.shape 2d");
      mx.resize(_bd.shape[0], _bd.shape[1]);
    }
    assert(mx.rows() == _bd.shape[0], mx.cols() == _bd.shape[1]);
    for (uint i0 = 0, N0 = mx.rows(); i0 < N0; ++i0) {
      for (uint i1 = 0, N1 = mx.cols(); i1 < N1; ++i1) {
        _is.read(reinterpret_cast<char*>(&mx(i0, i1)), TypeInfo<T>::size);
      }
    }
  }
}; // Loader ------------------

//
inline
std::ostream& operator<<(std::ostream& os, const BlockDescriptor& dr) {
  os << "Descriptor: name: " << dr.name
    << ", type: " << dr.type
    << ", typesize: " << dr.typesize
    << ", ndim: " << dr.ndim
    << ", shape: [";
    for (uint i = 0; i < dr.ndim; ++i) os << dr.shape[i] << ", ";
    os << "]"
    << ", numel: " << dr.numel;
  return os;
}

//
inline
std::ostream& print(std::ostream& os, ReadingDescriptor& rd) {
  rd.reset();
  while (rd.next()) {
    os << rd.block() << "\n";
  }
  rd.reset();
  return os;
}
}
