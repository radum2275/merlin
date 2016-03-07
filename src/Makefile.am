#######################################
# The list of libraries we are building seperated by spaces.
# The 'lib_' indicates that these build products will be installed
# in the $(libdir) directory. For example /usr/lib
lib_LTLIBRARIES = libmerlin.la

#######################################
# Build information for each library

# Sources for libmerlin
libmerlin_la_SOURCES = merlin.cpp graph.cpp run.cpp

# Linker options libTestProgram
libmerlin_la_LDFLAGS = -version-info 1:0:0

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libmerlin_la_CPPFLAGS = -I$(top_srcdir)/include

