#!/usr/bin/env python3
import sys
from pathlib import Path
from typing import Optional
import argparse

from acts.examples import (
    TGeoDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
    Interval,
    Sequencer,
)


import acts

from acts import MaterialMapJsonConverter, UnitConstants as u


def runITk(
    trackingGeometry,
    decorators,
    geometrySelection: Path,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    truthSmearedSeeded=False,
    truthEstimatedSeeded=False,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    s=None,
):

    s = s or Sequencer(events=100, numThreads=-1)

    logger = acts.logging.getLogger("CKFExample")

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)

    if inputParticlePath is None:
        logger.info("Generating particles using particle gun")

        evGen = acts.examples.EventGenerator(
            level=acts.logging.INFO,
            generators=[
                acts.examples.EventGenerator.Generator(
                    multiplicity=acts.examples.FixedMultiplicityGenerator(n=2),
                    vertex=acts.examples.GaussianVertexGenerator(
                        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
                    ),
                    particles=acts.examples.ParametricParticleGenerator(
                        p=(1 * u.GeV, 10 * u.GeV),
                        eta=(-2, 2),
                        phi=(0, 360 * u.degree),
                        randomizeCharge=True,
                        numParticles=4,
                    ),
                )
            ],
            outputParticles="particles_input",
            randomNumbers=rnd,
        )
        s.addReader(evGen)
        inputParticles = evGen.config.outputParticles
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        inputParticles = "particles_read"
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection=inputParticles,
                orderedEvents=False,
            )
        )

    # Selector
    selector = acts.examples.ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        outputParticles="particles_selected",
    )
    s.addAlgorithm(selector)

    # Simulation
    simAlg = acts.examples.FatrasSimulation(
        level=acts.logging.INFO,
        inputParticles=selector.config.outputParticles,
        outputParticlesInitial="particles_initial",
        outputParticlesFinal="particles_final",
        outputSimHits="simhits",
        randomNumbers=rnd,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        generateHitsOnSensitive=True,
    )
    s.addAlgorithm(simAlg)

    # Run the sim hits smearing
    digiCfg = acts.examples.DigitizationConfig(
        acts.examples.readDigiConfigFromJson(str(digiConfigFile)),
        trackingGeometry=trackingGeometry,
        randomNumbers=rnd,
        inputSimHits=simAlg.config.outputSimHits,
    )
    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, acts.logging.INFO)
    s.addAlgorithm(digiAlg)

    # Run the particle selection
    # The pre-selection will select truth particles satisfying provided criteria
    # from all particles read in by particle reader for further processing. It
    # has no impact on the truth hits themselves
    selAlg = acts.examples.TruthSeedSelector(
        level=acts.logging.INFO,
        ptMin=500 * u.MeV,
        nHitsMin=9,
        inputParticles=simAlg.config.outputParticlesInitial,
        inputMeasurementParticlesMap=digiCfg.outputMeasurementParticlesMap,
        outputParticles="particles_seed_selected",
    )
    s.addAlgorithm(selAlg)

    inputParticles = selAlg.config.outputParticles

    # Create starting parameters from either particle smearing or combined seed
    # finding and track parameters estimation
    if truthSmearedSeeded:
        logger.info("Using smeared truth particles for seeding")
        # Run particle smearing
        ptclSmear = acts.examples.ParticleSmearing(
            level=acts.logging.INFO,
            inputParticles=inputParticles,
            outputTrackParameters="smearedparameters",
            randomNumbers=rnd,
            # gaussian sigmas to smear particle parameters
            sigmaD0=20 * u.um,
            sigmaD0PtA=30 * u.um,
            sigmaD0PtB=0.3 / 1 * u.GeV,
            sigmaZ0=20 * u.um,
            sigmaZ0PtA=30 * u.um,
            sigmaZ0PtB=0.3 / 1 * u.GeV,
            sigmaPhi=1 * u.degree,
            sigmaTheta=1 * u.degree,
            sigmaPRel=0.01,
            sigmaT0=1 * u.ns,
            initialVarInflation=[1, 1, 1, 1, 1, 1],
        )
        outputTrackParameters = ptclSmear.config.outputTrackParameters
        s.addAlgorithm(ptclSmear)
    else:
        # Create space points
        spAlg = acts.examples.SpacePointMaker(
            level=acts.logging.INFO,
            inputSourceLinks=digiAlg.config.outputSourceLinks,
            inputMeasurements=digiAlg.config.outputMeasurements,
            outputSpacePoints="spacepoints",
            trackingGeometry=trackingGeometry,
            geometrySelection=acts.examples.readJsonGeometryList(
                str(geometrySelection)
            ),
        )
        s.addAlgorithm(spAlg)

        # Run either: truth track finding or seeding
        if truthEstimatedSeeded:
            logger.info("Using truth track finding from space points for seeding")
            # Use truth tracking
            truthTrackFinder = acts.examples.TruthTrackFinder(
                level=acts.logging.INFO,
                inputParticles=inputParticles,
                inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
                outputProtoTracks="prototracks",
            )
            s.addAlgorithm(truthTrackFinder)
            inputProtoTracks = truthTrackFinder.config.outputProtoTracks
            inputSeeds = ""
        else:
            logger.info("Using seeding")
            # Use seeding
            gridConfig = acts.SpacePointGridConfig(
                bFieldInZ=1.99724 * u.T,
                minPt=500 * u.MeV,
                rMax=200 * u.mm,
                zMax=2000 * u.mm,
                zMin=-2000 * u.mm,
                deltaRMax=60 * u.mm,
                cotThetaMax=7.40627,  # 2.7 eta
            )

            seedFilterConfig = acts.SeedFilterConfig(
                maxSeedsPerSpM=1, deltaRMin=1 * u.mm
            )

            seedFinderConfig = acts.SeedfinderConfig(
                rMax=gridConfig.rMax,
                deltaRMin=seedFilterConfig.deltaRMin,
                deltaRMax=gridConfig.deltaRMax,
                deltaRMinTopSP=seedFilterConfig.deltaRMin,
                deltaRMinBottomSP=seedFilterConfig.deltaRMin,
                deltaRMaxTopSP=gridConfig.deltaRMax,
                deltaRMaxBottomSP=gridConfig.deltaRMax,
                collisionRegionMin=-250 * u.mm,
                collisionRegionMax=250 * u.mm,
                zMin=gridConfig.zMin,
                zMax=gridConfig.zMax,
                maxSeedsPerSpM=seedFilterConfig.maxSeedsPerSpM,
                cotThetaMax=gridConfig.cotThetaMax,
                sigmaScattering=50,
                radLengthPerSeed=0.1,
                minPt=gridConfig.minPt,
                bFieldInZ=gridConfig.bFieldInZ,
                beamPos=acts.Vector2(0 * u.mm, 0 * u.mm),
                impactMax=3 * u.mm,
            )
            seeding = acts.examples.SeedingAlgorithm(
                level=acts.logging.INFO,
                inputSpacePoints=[spAlg.config.outputSpacePoints],
                outputSeeds="seeds",
                outputProtoTracks="prototracks",
                gridConfig=gridConfig,
                seedFilterConfig=seedFilterConfig,
                seedFinderConfig=seedFinderConfig,
            )
            s.addAlgorithm(seeding)
            inputProtoTracks = seeding.config.outputProtoTracks
            inputSeeds = seeding.config.outputSeeds

        # Write truth track finding / seeding performance
        trackFinderPerformanceWriter = acts.examples.TrackFinderPerformanceWriter(
            level=acts.logging.INFO,
            inputProtoTracks=inputProtoTracks,
            inputParticles=inputParticles,  # the original selected particles after digitization
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            filePath=str(outputDir / "performance_seeding_trees.root"),
        )
        s.addWriter(trackFinderPerformanceWriter)

        # Estimate track parameters from seeds
        paramEstimation = acts.examples.TrackParamsEstimationAlgorithm(
            level=acts.logging.INFO,
            inputSeeds=inputSeeds,
            inputProtoTracks=inputProtoTracks,
            inputSpacePoints=[spAlg.config.outputSpacePoints],
            inputSourceLinks=digiCfg.outputSourceLinks,
            outputTrackParameters="estimatedparameters",
            outputProtoTracks="prototracks_estimated",
            trackingGeometry=trackingGeometry,
            magneticField=field,
            bFieldMin=0.1 * u.T,
            deltaRMax=100.0 * u.mm,
            deltaRMin=10.0 * u.mm,
            sigmaLoc0=25.0 * u.um,
            sigmaLoc1=100.0 * u.um,
            sigmaPhi=0.02 * u.degree,
            sigmaTheta=0.02 * u.degree,
            sigmaQOverP=0.1 / 1.0 * u.GeV,
            sigmaT0=1400.0 * u.s,
            initialVarInflation=[1, 1, 1, 1, 1, 1],
        )
        s.addAlgorithm(paramEstimation)
        outputTrackParameters = paramEstimation.config.outputTrackParameters

    # Setup the track finding algorithm with CKF
    # It takes all the source links created from truth hit smearing, seeds from
    # truth particle smearing and source link selection config
    trackFinder = acts.examples.TrackFindingAlgorithm(
        level=acts.logging.INFO,
        measurementSelectorCfg=acts.MeasurementSelector.Config(
            [(acts.GeometryIdentifier(), ([], [15.0], [10]))]
        ),
        inputMeasurements=digiAlg.config.outputMeasurements,
        inputSourceLinks=digiAlg.config.outputSourceLinks,
        inputInitialTrackParameters=outputTrackParameters,
        outputTrajectories="trajectories",
        findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
            trackingGeometry, field
        ),
    )
    s.addAlgorithm(trackFinder)

    # write track states from CKF
    trackStatesWriter = acts.examples.RootTrajectoryStatesWriter(
        level=acts.logging.INFO,
        inputTrajectories=trackFinder.config.outputTrajectories,
        # @note The full particles collection is used here to avoid lots of warnings
        # since the unselected CKF track might have a majority particle not in the
        # filtered particle collection. This could be avoided when a seperate track
        # selection algorithm is used.
        inputParticles=selector.config.outputParticles,
        inputSimHits=simAlg.config.outputSimHits,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
        filePath=str(outputDir / "trackstates_ckf.root"),
        treeName="trackstates",
    )
    s.addWriter(trackStatesWriter)

    # write track summary from CKF
    trackSummaryWriter = acts.examples.RootTrajectorySummaryWriter(
        level=acts.logging.INFO,
        inputTrajectories=trackFinder.config.outputTrajectories,
        # @note The full particles collection is used here to avoid lots of warnings
        # since the unselected CKF track might have a majority particle not in the
        # filtered particle collection. This could be avoided when a seperate track
        # selection algorithm is used.
        inputParticles=selector.config.outputParticles,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        filePath=str(outputDir / "tracksummary_ckf.root"),
        treeName="tracksummary",
    )
    s.addWriter(trackSummaryWriter)

    # Write CKF performance data
    ckfPerfWriter = acts.examples.CKFPerformanceWriter(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        inputTrajectories=trackFinder.config.outputTrajectories,
        inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
        # The bottom seed could be the first, second or third hits on the truth track
        nMeasurementsMin=selAlg.config.nHitsMin - 3,
        ptMin=0.4 * u.GeV,
        filePath=str(outputDir / "performance_ckf.root"),
    )
    s.addWriter(ckfPerfWriter)

    if outputCsv:
        csv_dir = outputDir / "csv"
        csv_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Writing CSV files")
        csvMTJWriter = acts.examples.CsvMultiTrajectoryWriter(
            level=acts.logging.INFO,
            inputTrajectories=trackFinder.config.outputTrajectories,
            inputMeasurementParticlesMap=digiAlg.config.outputMeasurementParticlesMap,
            outputDir=str(csv_dir),
        )
        s.addWriter(csvMTJWriter)
        
    return s


def buildITkGeometry(geo_dir: Path, material: bool = True):
    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet
    equidistant = TGeoDetector.Config.BinningType.equidistant
    arbitrary = TGeoDetector.Config.BinningType.arbitrary

    logger = acts.logging.getLogger("buildITkGeometry")

    matDeco = None
    if material:
        file = geo_dir / "atlas/itk-hgtd/material-maps-ITk-HGTD.json"
        logger.info("Adding material from %s", file.absolute())
        matDeco = acts.IMaterialDecorator.fromFile(
            file,
            level=acts.logging.INFO,
        )

    return TGeoDetector.create(
        fileName=str(geo_dir / "atlas/itk-hgtd/ATLAS-ITk-HGTD.tgeo.root"),
        mdecorator=matDeco,
        buildBeamPipe=True,
        unitScalor=1.0,  # explicit units
        beamPipeRadius=23.934 * u.mm,
        beamPipeHalflengthZ=3000.0 * u.mm,
        beamPipeLayerThickness=0.8 * u.mm,
        surfaceLogLevel=acts.logging.WARNING,
        layerLogLevel=acts.logging.WARNING,
        volumeLogLevel=acts.logging.WARNING,
        volumes=[
            Volume(
                name="InnerPixels",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("Pixel::Pixel"),
                sensitiveNames=LayerTriplet(["Pixel::siLog"]),
                sensitiveAxes=LayerTriplet("YZX"),
                rRange=LayerTriplet((0 * u.mm, 135 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -250 * u.mm),
                    central=(-250 * u.mm, 250 * u.mm),
                    positive=(250 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(negative=-1.0, central=5 * u.mm, positive=-1.0),
                splitTolZ=LayerTriplet(
                    negative=10 * u.mm, central=-1.0, positive=10 * u.mm
                ),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(6, equidistant), (10, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(12, equidistant), (6, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),
            Volume(
                name="OuterPixels",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("Pixel::Pixel"),
                sensitiveNames=LayerTriplet(["Pixel::siLog"]),
                sensitiveAxes=LayerTriplet("YZX"),
                rRange=LayerTriplet((135 * u.mm, 350 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -377 * u.mm),
                    central=(-377 * u.mm, 377 * u.mm),
                    positive=(377 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=15 * u.mm, central=5 * u.mm, positive=15 * u.mm
                ),
                splitTolZ=LayerTriplet(
                    negative=20 * u.mm, central=-1.0, positive=20 * u.mm
                ),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),
            Volume(
                name="Strips",
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet(
                    negative="*",
                    central="SCT::SCT_Barrel",
                    positive="*",
                ),
                sensitiveNames=LayerTriplet(
                    negative=["SCT::ECSensor*"],
                    central=["SCT::BRLSensor*"],
                    positive=["SCT::ECSensor*"],
                ),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    negative=(-1.0, 1050 * u.mm),
                    central=(380 * u.mm, 1050 * u.mm),
                    positive=(-1.0, 1050 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -1400 * u.mm),
                    central=(-1400 * u.mm, 1400 * u.mm),
                    positive=(1400 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=-1.0,
                    central=35 * u.mm,
                    positive=-1.0,
                ),
                splitTolZ=LayerTriplet(
                    negative=35 * u.mm, central=-1.0, positive=35 * u.mm
                ),
                binning0=LayerTriplet(
                    negative=[(-1, arbitrary)],
                    central=[(0, equidistant)],
                    positive=[(-1, arbitrary)],
                ),
                binning1=LayerTriplet(
                    negative=[(-1, arbitrary)],
                    central=[(28, equidistant)] * 4,
                    positive=[(-1, arbitrary)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=True,
                barrelMap={"MS": 2, "SS": 4},
                discMap={
                    "EC0": [
                        [384.5, 403.481],
                        [403.481, 427.462],
                        [427.462, 456.442],
                        [456.442, 488.423],
                    ],
                    "EC1": [
                        [489.823, 507.916],
                        [507.916, 535.009],
                        [535.009, 559.101],
                        [559.101, 574.194],
                    ],
                    "EC2": [[575.594, 606.402], [606.402, 637.209]],
                    "EC3": [
                        [638.609, 670.832],
                        [670.832, 697.055],
                        [697.055, 723.278],
                        [723.278, 755.501],
                    ],
                    "EC4": [[756.901, 811.482], [811.482, 866.062]],
                    "EC5": [[867.462, 907.623], [907.623, 967.785]],
                },
            ),
            Volume(
                name="HGTD",
                binToleranceR=(15 * u.mm, 15 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.25 * u.mm, 0.25 * u.mm),
                layers=LayerTriplet(positive=True, central=False, negative=True),
                subVolumeName=LayerTriplet("HGTD::HGTD"),
                sensitiveNames=LayerTriplet(["HGTD::HGTDSiSensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    negative=(0 * u.mm, 1050 * u.mm),
                    positive=(0 * u.mm, 1050 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-4000 * u.mm, -3000 * u.mm),
                    positive=(3000 * u.mm, 4000 * u.mm),
                ),
                splitTolR=LayerTriplet(-1.0),
                splitTolZ=LayerTriplet(negative=10 * u.mm, positive=10 * u.mm),
                binning0=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                binning1=LayerTriplet(
                    negative=[(0, equidistant)],
                    central=[(0, equidistant)],
                    positive=[(0, equidistant)],
                ),
                cylinderDiscSplit=False,
                cylinderNZSegments=0,
                cylinderNPhiSegments=0,
                discNRSegments=0,
                discNPhiSegments=0,
                itkModuleSplit=False,
                barrelMap={},
                discMap={},
            ),
        ],
    )


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to construct the ITk geometry and write it out to CSV and OBJ formats"
    )
    p.add_argument(
        "geo_dir",
        help="Input directory containing the ITk standalone geometry. Get in touch if you don't have this.",
    )
    p.add_argument(
        "--output-dir",
        default=Path.cwd(),
        type=Path,
        help="Directory to write outputs to",
    )
    p.add_argument(
        "--output-csv", action="store_true", help="Write geometry in CSV format."
    )
    p.add_argument(
        "--output-obj", action="store_true", help="Write geometry in OBJ format."
    )
    p.add_argument(
        "--output-json",
        action="store_true",
        help="Write geometry and material in JSON format.",
    )
    p.add_argument(
        "--no-material", action="store_true", help="Decorate material to the geometry"
    )

    args = p.parse_args()
    args.output_dir.mkdir(exist_ok=True, parents=True)

    geo_example_dir = Path(args.geo_dir)
    assert geo_example_dir.exists(), "Detector example input directory missing"

    detector, trackingGeometry, decorators = buildITkGeometry(
        geo_example_dir,
        material=not args.no_material,
    )
    
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    runITk(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        field = field,
        geometrySelection= geo_example_dir / "atlas/itk-hgtd/geoSelection-ITk.json",
        digiConfigFile= geo_example_dir / "atlas/itk-hgtd/itk-smearing-config.json",
        outputCsv=True,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        inputParticlePath=inputParticlePath,
        outputDir=Path.cwd(),
    ).run()
