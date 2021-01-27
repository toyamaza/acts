# shared code 
add_subdirectory(Common)

add_executable(
  ActsExampleTruthTracksGeneric
  GenericRecTruthTracks.cpp)

target_link_libraries(
  ActsExampleTruthTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)

add_executable(
  ActsExampleCKFTracksGeneric
  GenericRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksGeneric
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorGeneric)

add_executable(
  ActsExampleCKFTracksTGeo
  TGeoRecCKFTracks.cpp)
target_link_libraries(
  ActsExampleCKFTracksTGeo
  PRIVATE
    ActsExamplesRecTracksCommon 
    ActsExamplesDetectorTGeo)

install(
  TARGETS
    ActsExampleTruthTracksGeneric
    ActsExampleCKFTracksGeneric
    ActsExampleCKFTracksTGeo
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

add_subdirectory_if(DD4hep ACTS_BUILD_EXAMPLES_DD4HEP)
