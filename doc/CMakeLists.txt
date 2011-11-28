cmake_minimum_required(VERSION 2.6)
 
PROJECT(openvoronoi_documentation)


add_custom_target(
    pp_edge_fig 
    ALL
    DEPENDS pp_edge.eps 
)
add_custom_command(
    OUTPUT pp_edge.eps
    COMMAND asy ${CMAKE_CURRENT_SOURCE_DIR}/pp_edge.asy 
)


add_custom_target(
    pdf_doc
    ALL
    DEPENDS 
    openvoronoi.pdf
    pp_edge.eps
)


add_custom_command(
    OUTPUT openvoronoi.pdf
    COMMAND lyx --export pdf ${CMAKE_CURRENT_SOURCE_DIR}/openvoronoi.lyx
    DEPENDS pp_edge.eps

)

# ADD_EXECUTABLE(diangle_test diangle_test.cpp)



