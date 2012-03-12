if (BUILD_DOC)

FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
    MESSAGE(STATUS "Found doxygen. Documentation can be built with 'make doc' ")
    
    # this works around a bug in cmake 2.4 (Ubuntu Hardy)
    execute_process(
        COMMAND mkdir -p doc/html doc/latex
    )

    FIND_PACKAGE(LATEX)
    IF    (NOT LATEX_COMPILER)
        MESSAGE(STATUS "latex command LATEX_COMPILER not found but usually required. You will probably get warnings and user inetraction on doxy run.")
    ENDIF (NOT LATEX_COMPILER)
    IF    (NOT MAKEINDEX_COMPILER)
        MESSAGE(STATUS "makeindex command MAKEINDEX_COMPILER not found but usually required.")
    ENDIF (NOT MAKEINDEX_COMPILER)
    IF    (NOT DVIPS_CONVERTER)
        MESSAGE(STATUS "dvips command DVIPS_CONVERTER not found but usually required.")
    ENDIF (NOT DVIPS_CONVERTER)
    
    #if (EXISTS Doxyfile)
        set(DOXY_CONFIG ${CMAKE_SOURCE_DIR}/Doxyfile)
    #endif (EXISTS Doxyfile)

    add_custom_command(
        OUTPUT
            doc/latex/index.tex
            doc/html/index.html
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXY_CONFIG}
        COMMENT building LaTex & HTML docs
    )

    add_custom_target(
        doc
        DEPENDS doc/latex/index.tex
    )


    #IF (EXISTS ${PDFLATEX_COMPILER})
    #    add_custom_command(
    #        OUTPUT doc/latex/refman.pdf
    #        DEPENDS doc/latex/index.tex
    #        WORKING_DIRECTORY doc/latex
    #        COMMAND make pdf
    #        COMMENT building PDF docs
    #        COMMAND mv refman.pdf ../openvoronoi-manual.pdf
    #    )

    #    add_custom_target(
    #        doc-pdf
    #        DEPENDS doc/latex/refman.pdf
    #    )

    #    add_dependencies(doc doc-pdf)
    #ELSE (EXISTS ${PDFLATEX_COMPILER}) 
    #   message(STATUS "pdflatex compiler not found, PDF docs will not be built")
    #ENDIF (EXISTS ${PDFLATEX_COMPILER})


    #add_custom_target(
    #    doc-latex
    #    DEPENDS doc/latex/index.tex
    #)

    #install(
    #    DIRECTORY doc/latex/
    #    DESTINATION share/doc/python-opencam/pdf
    #    FILES_MATCHING PATTERN *.pdf
    #)

    #install(
    #    FILES doc/ocl-manual.pdf
    #    DESTINATION share/doc/python-opencam/pdf
    #)

    #install(
    #    DIRECTORY doc/html
    #    DESTINATION share/doc/python-opencam/
    #)

ENDIF(DOXYGEN_FOUND)
endif (BUILD_DOC)
