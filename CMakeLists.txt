cmake_minimum_required(VERSION 3.0.0)
project(merlin VERSION 0.1.0)

include(CTest)
enable_testing()
set (CMAKE_CXX_STANDARD 11)
find_package(Boost COMPONENTS program_options thread REQUIRED)

include_directories(include)

add_executable(merlin 
    src/main.cpp
    src/bte.cpp
    src/cte.cpp
    src/em.cpp
    src/gibbs.cpp
    src/graph.cpp
    src/ijgp.cpp
    src/jglp.cpp
    src/lbp.cpp
    src/merlin.cpp
    src/program_options.cpp
    src/util.cpp
    src/wmb.cpp
    include/algorithm.h
    include/base.h
    include/bte.h
    include/cte.h
    include/em.h
    include/enum.h
    include/factor_graph.h
    include/factor.h
    include/gibbs.h
    include/graph.h
    include/graphical_model.h
    include/ijgp.h
    include/index.h
    include/indexed_heap.h
    include/jglp.h
    include/lbp.h
    include/merlin.h
    include/observation.h
    include/program_options.h
    include/set.h
    include/util.h
    include/variable_set.h
    include/variable.h
    include/vector.h
    include/wmb.h    
)

set(LINK_LIBS dl ${Boost_LIBRARIES})

target_link_libraries(${PROJECT_NAME} ${LINK_LIBS} ${LINK_LIBS})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
