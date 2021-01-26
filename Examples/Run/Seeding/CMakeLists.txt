# shared code 
add_subdirectory(Common)

add_executable(
  ActsExampleSeedingGeneric
  GenericSeedingExample.cpp)

target_link_libraries(
  ActsExampleSeedingGeneric
  PRIVATE
    ActsExamplesSeedingCommon
    ActsExamplesDetectorGeneric    
)

add_executable(
  ActsExampleSeedingTGeo
  TGeoSeedingExample.cpp)

target_link_libraries(
  ActsExampleSeedingTGeo
  PRIVATE
    ActsExamplesSeedingCommon
    ActsExamplesDetectorTGeo
)


install(
  TARGETS
    ActsExampleSeedingGeneric
    ActsExampleSeedingTGeo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
