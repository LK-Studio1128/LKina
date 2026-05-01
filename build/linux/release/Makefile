BASE=/usr
BOOST_VERSION=
BOOST_INCLUDE = $(BASE)/include
C_PLATFORM=-pthread
GPP=g++
C_OPTIONS= -O3 -DNDEBUG -std=c++14
BOOST_LIB_VERSION=
BOOST_STATIC=y

include ../../makefile_common

# Debian/Ubuntu x86_64: static Boost lives in multiarch dir, not /usr/lib directly.
# Override the absolute .a paths set by makefile_common.
BOOST_LIBDIR := $(or $(wildcard /usr/lib/x86_64-linux-gnu),/usr/lib)
_MT := $(if $(wildcard $(BOOST_LIBDIR)/libboost_thread-mt.a),-mt,)
LIBS := $(wildcard $(BOOST_LIBDIR)/libboost_system$(_MT).a) \
        $(BOOST_LIBDIR)/libboost_thread$(_MT).a \
        $(BOOST_LIBDIR)/libboost_serialization$(_MT).a \
        $(BOOST_LIBDIR)/libboost_filesystem$(_MT).a \
        $(BOOST_LIBDIR)/libboost_program_options$(_MT).a
