/*
 * enum.h
 *
 *  Created on: Feb 8, 2013
 *      Author: radu
 */


/// \file enum.h
/// \brief Type-safe enumeration class with string-ize functions
/// \author Radu Marinescu 
///

#ifndef __IBM_MERLIN_ENUM_H
#define __IBM_MERLIN_ENUM_H

#include <cstring>
#include <iostream>

/* Type-safe enumeration class with stringize functions */

namespace merlin {

#define MER_ENUM(enum_name, v0, ...) \
  struct enum_name {                   \
    enum Type { v0, __VA_ARGS__ };    \
    Type t_;                          \
    enum_name(Type t=v0) : t_(t) {} \
    enum_name(const char* s) : t_() {  \
      /*for (int i=0;i<nType;++i) if (strcasecmp(names()[i],s)==0) { t_=Type(i); return; } */ \
      for (int i=0;names(i)[0]!=0;++i) if (strncmp(names(i),s,strlen(s))==0) { t_=Type(i); return; } \
      throw std::runtime_error("Unknown type string");                                   \
    }                                 \
    operator Type () const { return t_; }                    \
    operator char const* () const { return names(t_); }    \
    friend std::istream& operator >> (std::istream& is, enum_name& t) {    \
      std::string str; is >> str; t = enum_name(str.c_str()); return is; } \
    friend std::ostream& operator << (std::ostream& os, enum_name& t) {    \
      const char* s=(const char*)t; while (*s!=',') os<<*s++; return os; \
      /* os << (const char*)t; return os;                                */ \
    } \
private:                                                        \
    static char const* names(unsigned int i) {        \
      static char const str[] = { #v0 "," #__VA_ARGS__ ",\0" }; \
      char const* s=str; while (*s!=0 && i!=0) if (*(s++)==',') --i; \
      return s;                                              \
    }                                 \
  };                               

} // namespace

#endif // include
