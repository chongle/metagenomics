project(metagenomics)

set(CMAKE_CXX_FLAGS "-Wall ")

add_executable(metagenomics ${flags} Edge.cpp HashTable.cpp main.cpp OverlapGraph.cpp OverlapGraphExt.cpp Read.cpp Dataset.cpp MatePairGraph.cpp Edge.h HashTable.h OverlapGraph.h Read.h Common.h Dataset.h MatePairGraph.h)

install(TARGETS metagenomics RUNTIME DESTINATION bin)
