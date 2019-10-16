/*
 * enum.h
 *
 *  Created on: Feb 8, 2013
 *      Author: radu
 *
 * Copyright (c) 2015, International Business Machines Corporation
 * and University of California Irvine. All rights reserved.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/// \file enum.h
/// \brief Type-safe enumeration class with string-ize functions
/// \author Radu Marinescu radu.marinescu@ie.ibm.com
///

#ifndef IBM_MERLIN_ENUM_H_
#define IBM_MERLIN_ENUM_H_

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
