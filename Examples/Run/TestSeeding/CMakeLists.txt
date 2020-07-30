add_executable(
  ActsExampleTestSeedAlgorithm
  TestSeedRun.cpp
  TestSeedAlgorithm.cpp)
target_link_libraries(
  ActsExampleTestSeedAlgorithm
  PRIVATE
    ActsCore
    ActsExamplesFramework 
    ActsExamplesCommon 
    ActsExamplesDetectorGeneric
    ActsExamplesMagneticField
    ActsExamplesIoCsv
    Boost::program_options)

install(
  TARGETS ActsExampleTestSeedAlgorithm
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})