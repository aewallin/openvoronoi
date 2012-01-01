cmake_minimum_required(VERSION 2.6)
 
PROJECT(openvoronoi_documentation)

set( ASY_FIGURES
    pp_edge
    p_site
    l_site
    ll_edge
)

FOREACH(ASYFIG ${ASY_FIGURES})
    set(ASY_FILES ${ASY_FILES} ";" ${ASYFIG}.asy) 
    set(ASY_EPS_FILES ${ASY_EPS_FILES} ";" ${ASYFIG}.eps) 
ENDFOREACH(ASYFIG)

add_custom_target(
    figures
    ALL
    DEPENDS 
    ${ASY_EPS_FILES}
)

FOREACH(ASYFIG ${ASY_FIGURES})
    message( STATUS "figure: " ${ASYFIG} )
    message( STATUS "CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )
    message( STATUS "CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )
    add_custom_command(
        OUTPUT ${ASYFIG}.eps
        COMMAND asy ${CMAKE_CURRENT_SOURCE_DIR}/${ASYFIG}.asy 
        DEPENDS ${ASYFIG}.asy
    )
ENDFOREACH(ASYFIG)

add_custom_target(
    pdf_doc
    ALL
    DEPENDS 
    ${ASY_EPS_FILES}
    openvoronoi.pdf
)

ADD_CUSTOM_COMMAND(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/openvoronoi.pdf
   COMMAND ${CMAKE_COMMAND} -E create_symlink
       ${CMAKE_CURRENT_SOURCE_DIR}/openvoronoi.lyx
       ${CMAKE_CURRENT_BINARY_DIR}/openvoronoi.lyx
   COMMAND lyx --export pdf ${CMAKE_CURRENT_BINARY_DIR}/openvoronoi.lyx
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/openvoronoi.lyx
   ${ASY_EPS_FILES}
)

# add_custom_command(
#     OUTPUT openvoronoi.pdf
#     COMMAND lyx --export pdf ${CMAKE_CURRENT_SOURCE_DIR}/openvoronoi.lyx
#     DEPENDS openvoronoi.lyx
# )
