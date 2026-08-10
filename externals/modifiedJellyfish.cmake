SET_PROPERTY(DIRECTORY PROPERTY "EP_BASE" ${ep_base})

SET(MODIFIED_JELLYFISH_PROJECT modified_jellyfish_project CACHE INTERNAL "modifiedJellyfish project name")
SET(MODIFIED_JELLYFISH_DIR ${CMAKE_BINARY_DIR}/externals/modified_jellyfish CACHE INTERNAL "modifiedJellyfish project directory")
SET(MODIFIED_JELLYFISH_LIB)



ExternalProject_Add(${MODIFIED_JELLYFISH_PROJECT}
	# Source lives in the src/modifiedJellyfish git submodule (repo: stefinfection/modified-jellyfish).
	# Copy it into the build tree and build there (BUILD_IN_SOURCE) so the tracked submodule working
	# copy stays pristine -- no configure/make artifacts leak back into it.
	DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory ${PROJECT_SOURCE_DIR}/src/modifiedJellyfish ${PROJECT_SOURCE_DIR}/bin/externals/modified_jellyfish/src/modified_jellyfish_project

	# git checkout / copy_directory reset file mtimes, so make would see configure as older than
	# configure.ac and try to regenerate it with autoconf (which fails). Touch the generated
	# autotools files AFTER their sources so they look up-to-date and no regeneration is attempted.
	# (The old tarball path avoided this because tar preserved the original mtime ordering.)
	PATCH_COMMAND bash -c "touch configure.ac aclocal.m4 && find . -name Makefile.am -exec touch {} + && touch configure config.h.in && find . -name Makefile.in -exec touch {} +"

        CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/bin/externals/modified_jellyfish/src/modified_jellyfish_project/configure --prefix=${PROJECT_SOURCE_DIR}/bin/externals/modified_jellyfish/src/modified_jellyfish_project/
        BUILD_IN_SOURCE 1
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        UPDATE_COMMAND ""
        PREFIX ${MODIFIED_JELLYFISH_DIR}
)

ExternalProject_Get_Property(${MODIFIED_JELLYFISH_PROJECT} INSTALL_DIR)
ExternalProject_Get_Property(${MODIFIED_JELLYFISH_PROJECT} SOURCE_DIR)
ExternalProject_Get_Property(${MODIFIED_JELLYFISH_PROJECT} BINARY_DIR)

SET(MODIFIED_JELLYFISH_INCLUDE ${SOURCE_DIR} CACHE INTERNAL "MODIFIED JELLYFISH INCLUDE")