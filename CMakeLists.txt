cmake_minimum_required(VERSION 3.12)


#
# Fortran is needed so we get -lgfortran when linking to some
# blas/lapack implementations
#
if(NOT DEFINED FFPOPTVERSION)
	set(FFPOPTVERSION 0.0)
endif()

project(ffpopt VERSION "${FFPOPTVERSION}" LANGUAGES CXX)


include(FetchContent)

if(FALSE)  
  if (SKBUILD)
    #
    # This doesn't seem to install the edgembar companion library
    # How do we do that?
    #
    install(FILES src/python/bin/edgembar-amber2dats.py
      src/python/bin/edgembar-WriteGraphHtml.py
      src/python/bin/edgembar-calcamberboresch.py
      DESTINATION bin
      PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
    #
    # We can't install scripts.py to the installation prefix
    # That would be, like, installing to /usr/scripts.py
    #
    install(FILES src/python/scripts.py DESTINATION .)
  endif()
endif()

  
if(BUILD_PYTHON)
  add_subdirectory(src/python)
endif()

