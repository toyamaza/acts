// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilder.hpp"
#include "Acts/SpacePointFormation/SpacePointBuilderConfig.h"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/DetectorElementStub.hpp"
#include "Acts/Tests/CommonHelpers/GenerateParameters.hpp"
#include "Acts/Tests/CommonHelpers/MeasurementsCreator.hpp"
#include "Acts/Tests/CommonHelpers/TestSpacePoint.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <limits>
#include <variant>
namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

using namespace Acts::UnitLiterals;

using StraightPropagator =
    Acts::Propagator<Acts::StraightLineStepper, Acts::Navigator>;

using TestMeasurement = Acts::BoundVariantMeasurement;
using ConstantFieldStepper = Acts::EigenStepper<>;
using ConstantFieldPropagator =
    Acts::Propagator<ConstantFieldStepper, Acts::Navigator>;
// Construct initial track parameters.
CurvilinearTrackParameters makeParameters(double phi, double theta, double p,
                                          double q) {
  // create covariance matrix from reasonable standard deviations
  Acts::BoundVector stddev;
  stddev[Acts::eBoundLoc0] = 100_um;
  stddev[Acts::eBoundLoc1] = 100_um;
  stddev[Acts::eBoundTime] = 25_ns;
  stddev[Acts::eBoundPhi] = 2_degree;
  stddev[Acts::eBoundTheta] = 2_degree;
  stddev[Acts::eBoundQOverP] = 1 / 100_GeV;
  BoundSymMatrix cov = stddev.cwiseProduct(stddev).asDiagonal();
  // Let the particle starts from the origin
  Vector4 mPos4(-3_m, 0., 0., 0.);
  return CurvilinearTrackParameters(mPos4, phi, theta, p, q, cov);
}

// Create a test context
GeometryContext tgContext = GeometryContext();

const GeometryContext geoCtx;
const MagneticFieldContext magCtx;
const CalibrationContext calCtx;

// detector geometry
CubicTrackingGeometry geometryStore(geoCtx);
const auto geometry = geometryStore();

// detector resolutions
const MeasurementResolution resPixel = {MeasurementType::eLoc01,
                                        {25_um, 50_um}};
const MeasurementResolution resStrip0 = {MeasurementType::eLoc0, {100_um}};
const MeasurementResolution resStrip1 = {MeasurementType::eLoc0, {100_um}};
const MeasurementResolutionMap resolutions = {
    {GeometryIdentifier().setVolume(2), resPixel},
    {GeometryIdentifier().setVolume(3).setLayer(2), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(4), resStrip1},
    {GeometryIdentifier().setVolume(3).setLayer(6), resStrip0},
    {GeometryIdentifier().setVolume(3).setLayer(8), resStrip1},
};

// Construct a straight-line propagator.
static StraightPropagator makeStraightPropagator(
    std::shared_ptr<const Acts::TrackingGeometry> geo) {
  Acts::Navigator::Config cfg{geo};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator{cfg};
  Acts::StraightLineStepper stepper;
  return StraightPropagator(std::move(stepper), std::move(navigator));
}

// simulation propagator
const auto measPropagator = makeStraightPropagator(geometry);

std::default_random_engine rng(42);

BOOST_DATA_TEST_CASE(SpacePointBuilder_basic, bdata::xrange(1), index) {
  (void)index;

  double phi = 5._degree;
  double theta = 95._degree;
  double p = 20._GeV;
  double q = 1;

  Acts::Navigator navigator({
      geometry,
      true,  // sensitive
      true,  // material
      false  // passive
  });
  auto field =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, 2._T));
  ConstantFieldStepper stepper(std::move(field));

  ConstantFieldPropagator propagator(std::move(stepper), std::move(navigator));
  auto start = makeParameters(phi, theta, p, q);

  auto measurements =
      createMeasurements(propagator, geoCtx, magCtx, start, resolutions, rng);

  auto sourceLinks = measurements.sourceLinks;

  std::vector<const TestMeasurement*> frontMeasurements;
  std::vector<const TestMeasurement*> backMeasurements;
  std::vector<const TestMeasurement*> singleHitMeasurements;
  for (auto& sl : sourceLinks) {
    const auto geoId = sl.geometryId();

    const TestMeasurement* meas = new TestMeasurement(makeMeasurement(
        sl, sl.parameters, sl.covariance, sl.indices[0], sl.indices[1]));

    const auto volumeId = geoId.volume();

    if (volumeId == 2) {  // pixel type detector
      singleHitMeasurements.emplace_back(meas);
    } else if (volumeId == 3) {  // strip type detector

      const auto layerId = geoId.layer();

      if (layerId == 2 || layerId == 6) {
        frontMeasurements.emplace_back(meas);
      } else if (layerId == 4 || layerId == 8) {
        backMeasurements.emplace_back(meas);
      }

    }  // volume 3 (strip detector)
  }
  auto spBuilderConfig = SpacePointBuilderConfig();
  spBuilderConfig.trackingGeometry = geometry;

  auto spBuilder = Acts::SpacePointBuilder<TestSpacePoint>(spBuilderConfig);

  TestSpacePointContainer spacePoints;
  std::cout << "number of front/back measurements " << frontMeasurements.size()
            << " / " << backMeasurements.size() << std::endl;
  // pixel SP building
  spBuilder.calculateSpacePoints(tgContext, spacePoints,
                                 &singleHitMeasurements);

  // strip SP building
  spBuilder.calculateSpacePoints(tgContext, spacePoints, &frontMeasurements,
                                 &backMeasurements);

  std::cout << "Number of space points " << spacePoints.size() << std::endl;

  for (auto& sp : spacePoints) {
    std::cout << "space point (" << sp.x() << " " << sp.y() << " " << sp.z()
              << ") var: " << sp.varianceR() << " " << sp.varianceZ()
              << std::endl;
  }
  std::cout << "Space point calculated" << std::endl;
}

}  // end of namespace Test
}  // namespace Acts
