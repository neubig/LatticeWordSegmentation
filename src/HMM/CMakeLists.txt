add_library(HMM
  DiscreteHMM.cpp
)

target_link_libraries(HMM
  fst
  dl
)
