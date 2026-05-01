BASE=$(shell brew --prefix)
BOOST_VERSION=
BOOST_INCLUDE = $(BASE)/include
C_PLATFORM=-pthread
GPP=/usr/bin/clang++
C_OPTIONS= -O3 -DNDEBUG -std=c++14 -fvisibility=hidden
BOOST_LIB_VERSION=
BOOST_STATIC=y

# OpenMP via Homebrew libomp (brew install libomp)
LIBOMP_PREFIX=$(shell brew --prefix libomp 2>/dev/null)
ifneq ($(LIBOMP_PREFIX),)
  C_OPTIONS    += -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
  EXTRA_LDFLAGS = -L$(LIBOMP_PREFIX)/lib
  EXTRA_LIBS    = $(LIBOMP_PREFIX)/lib/libomp.a
endif

include ../../makefile_common
