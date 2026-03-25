// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "hostname.hh"

namespace Dune {
  namespace PDELab {

    //! C++ friendly wrapper around POSIX' gethostname() or GetComputerName() when compiling for Windows applications
    std::string getHostName() {

#if !defined(_WIN32)
      std::size_t bufsize = 1024;
      std::vector<char> buf(bufsize);
      while(gethostname(&buf[0], buf.size()),
            buf.back() = '\0',
            std::strlen(&buf[0]) == buf.size()-1)
      {
        buf.clear();
        buf.resize(bufsize*=2);
      }
      #else
      DWORD bufsize = MAX_COMPUTERNAME_LENGTH + 1;
      std::vector<char> buf(bufsize);
      if(!GetComputerNameA(&buf[0], &bufsize))
        return {};
      buf[bufsize] = '\0';
      #endif

      // ignore everything after the first '.', if any
      std::vector<char>::iterator end = buf.begin();
      while(*end != '\0' && *end != '.') ++end;
      std::string str(buf.begin(), end);
      return str;
    }

  } // namespace PDELab
} // namespace Dune
