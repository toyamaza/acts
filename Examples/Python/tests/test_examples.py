from pathlib import Path
import os
import json
import functools
import subprocess

import pytest

from helpers import (
    geant4Enabled,
    rootEnabled,
    dd4hepEnabled,
    hepmc3Enabled,
    AssertCollectionExistsAlg,
    isCI,
    doHashChecks,
)

pytestmark = pytest.mark.skipif(not rootEnabled, reason="ROOT not set up")


import acts
from acts.examples import (
    Sequencer,
    GenericDetector,
    AlignedDetector,
    RootParticleWriter,
)

from common import getOpenDataDetector

u = acts.UnitConstants


@pytest.fixture
def field():
    return acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))


@pytest.fixture
def seq():
    return Sequencer(events=10, numThreads=1)


def assert_csv_output(csv_path, stem):
    __tracebackhide__ = True
    # print(list(csv_path.iterdir()))
    assert len([f for f in csv_path.iterdir() if f.name.endswith(stem + ".csv")]) > 0
    assert all([f.stat().st_size > 100 for f in csv_path.iterdir()])


def assert_entries(root_file, tree_name, exp):
    __tracebackhide__ = True
    import ROOT

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(True)

    rf = ROOT.TFile.Open(str(root_file))
    keys = [k.GetName() for k in rf.GetListOfKeys()]
    assert tree_name in keys
    assert rf.Get(tree_name).GetEntries() == exp, f"{root_file}:{tree_name}"


@pytest.mark.slow
def test_pythia8(tmp_path, seq):
    from pythia8 import runPythia8

    (tmp_path / "csv").mkdir()

    assert not (tmp_path / "pythia8_particles.root").exists()
    assert len(list((tmp_path / "csv").iterdir())) == 0

    events = seq.config.events

    runPythia8(str(tmp_path), s=seq).run()

    del seq

    fp = tmp_path / "pythia8_particles.root"
    assert fp.exists()
    assert fp.stat().st_size > 2**10 * 50
    assert_entries(fp, "particles", events)

    assert len(list((tmp_path / "csv").iterdir())) > 0
    assert_csv_output(tmp_path / "csv", "particles")


def test_fatras(trk_geo, tmp_path, field, assert_root_hash):
    from fatras import runFatras

    csv = tmp_path / "csv"
    csv.mkdir()

    nevents = 10

    root_files = [
        (
            "fatras_particles_final.root",
            "particles",
            nevents,
        ),
        (
            "fatras_particles_initial.root",
            "particles",
            nevents,
        ),
        (
            "hits.root",
            "hits",
            115,
        ),
    ]

    assert len(list(csv.iterdir())) == 0
    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    seq = Sequencer(events=nevents)
    runFatras(trk_geo, field, str(tmp_path), s=seq).run()

    del seq

    assert_csv_output(csv, "particles_final")
    assert_csv_output(csv, "particles_initial")
    assert_csv_output(csv, "hits")
    for f, tn, exp_entries in root_files:
        rfp = tmp_path / f
        assert rfp.exists()
        assert rfp.stat().st_size > 2**10 * 10

        assert_entries(rfp, tn, exp_entries)
        assert_root_hash(f, rfp)


def test_seeding(tmp_path, trk_geo, field, assert_root_hash):
    from seeding import runSeeding

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    csv = tmp_path / "csv"
    csv.mkdir()

    seq = Sequencer(events=10, numThreads=1)

    root_files = [
        (
            "estimatedparams.root",
            "estimatedparams",
            371,
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_tracks",
            371,
        ),
        (
            "performance_seeding_hists.root",
            None,
            0,
        ),
        (
            "evgen_particles.root",
            "particles",
            seq.config.events,
        ),
        (
            "fatras_particles_final.root",
            "particles",
            seq.config.events,
        ),
        (
            "fatras_particles_initial.root",
            "particles",
            seq.config.events,
        ),
    ]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(csv.iterdir())) == 0

    runSeeding(trk_geo, field, outputDir=str(tmp_path), s=seq).run()

    del seq

    for fn, tn, exp_entries in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 100

        if tn is not None:
            assert_entries(fp, tn, exp_entries)
            assert_root_hash(fn, fp)

    assert_csv_output(csv, "evgen_particles")
    assert_csv_output(csv, "evgen_particles")
    assert_csv_output(csv, "fatras_particles_final")
    assert_csv_output(csv, "fatras_particles_initial")


def test_propagation(tmp_path, trk_geo, field, seq, assert_root_hash):
    from propagation import runPropagation

    obj = tmp_path / "obj"
    obj.mkdir()

    root_files = [
        (
            "propagation_steps.root",
            "propagation_steps",
            10000,
        )
    ]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    assert len(list(obj.iterdir())) == 0

    runPropagation(trk_geo, field, str(tmp_path), s=seq).run()

    for fn, tn, ee in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 2**10 * 50
        assert_entries(fp, tn, ee)
        assert_root_hash(fn, fp)

    assert len(list(obj.iterdir())) > 0


@pytest.mark.slow
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_recording(tmp_path, material_recording, assert_root_hash):

    root_files = [
        (
            "geant4_material_tracks.root",
            "material-tracks",
            200,
        )
    ]

    for fn, tn, ee in root_files:
        fp = material_recording / fn
        assert fp.exists()
        assert fp.stat().st_size > 2**10 * 50
        assert_entries(fp, tn, ee)
        assert_root_hash(fn, fp)


@pytest.mark.slow
@pytest.mark.skipif(not hepmc3Enabled, reason="HepMC3 plugin not available")
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_event_recording(tmp_path):

    script = (
        Path(__file__).parent.parent.parent.parent
        / "Examples"
        / "Scripts"
        / "Python"
        / "event_recording.py"
    )
    assert script.exists()

    env = os.environ.copy()
    env["NEVENTS"] = "1"
    subprocess.check_call([str(script)], cwd=tmp_path, env=env)

    from acts.examples.hepmc3 import HepMC3AsciiReader

    out_path = tmp_path / "hepmc3"
    # out_path.mkdir()

    assert len([f for f in out_path.iterdir() if f.name.endswith("events.hepmc3")]) > 0
    assert all([f.stat().st_size > 100 for f in out_path.iterdir()])

    s = Sequencer(numThreads=1)

    s.addReader(
        HepMC3AsciiReader(
            level=acts.logging.INFO,
            inputDir=str(out_path),
            inputStem="events",
            outputEvents="hepmc-events",
        )
    )

    alg = AssertCollectionExistsAlg(
        "hepmc-events", name="check_alg", level=acts.logging.INFO
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 1


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.parametrize("revFiltMomThresh", [0 * u.GeV, 1 * u.TeV])
def test_truth_tracking(tmp_path, assert_root_hash, revFiltMomThresh):
    from truth_tracking import runTruthTracking

    detector, trackingGeometry, _ = getOpenDataDetector()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    seq = Sequencer(events=10, numThreads=1)

    root_files = [
        ("trackstates_fitter.root", "trackstates", 19),
        ("tracksummary_fitter.root", "tracksummary", 10),
        ("performance_track_finder.root", "track_finder_tracks", 19),
        ("performance_track_fitter.root", None, -1),
    ]

    for fn, _, _ in root_files:
        fp = tmp_path / fn
        assert not fp.exists()

    runTruthTracking(
        trackingGeometry,
        field,
        digiConfigFile=Path(
            "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json",
        ),
        outputDir=tmp_path,
        reverseFilteringMomThreshold=revFiltMomThresh,
        s=seq,
    )

    seq.run()

    del seq

    for fn, tn, ee in root_files:
        fp = tmp_path / fn
        assert fp.exists()
        assert fp.stat().st_size > 1024
        if tn is not None:
            assert_entries(fp, tn, ee)
            assert_root_hash(fn, fp)


def test_particle_gun(tmp_path, assert_root_hash):
    from particle_gun import runParticleGun

    s = Sequencer(events=20, numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "particles.root"

    assert not csv_dir.exists()
    assert not root_file.exists()

    runParticleGun(str(tmp_path), s=s).run()

    assert csv_dir.exists()
    assert root_file.exists()

    assert len([f for f in csv_dir.iterdir() if f.name.endswith("particles.csv")]) > 0
    assert all([f.stat().st_size > 100 for f in csv_dir.iterdir()])

    assert root_file.stat().st_size > 200
    assert_entries(root_file, "particles", 20)
    assert_root_hash(root_file.name, root_file)


@pytest.mark.slow
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_material_mapping(material_recording, tmp_path, assert_root_hash):
    map_file = tmp_path / "material-map_tracks.root"
    assert not map_file.exists()

    s = Sequencer(numThreads=1)

    detector, trackingGeometry, decorators = getOpenDataDetector()

    from material_mapping import runMaterialMapping

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir=str(tmp_path),
        inputDir=material_recording,
        s=s,
    )

    s.run()

    # MaterialMapping alg only writes on destruct.
    # See https://github.com/acts-project/acts/issues/881
    del s

    mat_file = tmp_path / "material-map.json"

    assert mat_file.exists()
    assert mat_file.stat().st_size > 10

    with mat_file.open() as fh:
        assert json.load(fh)

    assert map_file.exists()
    assert_entries(map_file, "material-tracks", 200)
    assert_root_hash(map_file.name, map_file)

    val_file = tmp_path / "propagation-material.root"
    assert not val_file.exists()

    # test the validation as well

    # we need to destroy the ODD to reload with material
    # del trackingGeometry
    # del detector

    detector, trackingGeometry, decorators = getOpenDataDetector(
        mdecorator=acts.IMaterialDecorator.fromFile(mat_file)
    )

    from material_validation import runMaterialValidation

    s = Sequencer(events=10, numThreads=1)

    field = acts.NullBField()

    runMaterialValidation(
        trackingGeometry, decorators, field, outputDir=str(tmp_path), s=s
    )

    s.run()

    assert val_file.exists()
    assert_entries(val_file, "material-tracks", 10000)
    assert_root_hash(val_file.name, val_file)


@pytest.mark.slow
@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
def test_volume_material_mapping(material_recording, tmp_path, assert_root_hash):
    map_file = tmp_path / "material-map-volume_tracks.root"
    assert not map_file.exists()

    s = Sequencer(numThreads=1)

    geo_map = Path(__file__).parent / "geometry-volume-map.json"
    assert geo_map.exists()
    assert geo_map.stat().st_size > 10
    with geo_map.open() as fh:
        assert json.load(fh)

    detector, trackingGeometry, decorators = getOpenDataDetector(
        mdecorator=acts.IMaterialDecorator.fromFile(geo_map)
    )

    from material_mapping import runMaterialMapping

    runMaterialMapping(
        trackingGeometry,
        decorators,
        mapName="material-map-volume",
        outputDir=str(tmp_path),
        inputDir=material_recording,
        s=s,
    )

    s.run()

    # MaterialMapping alg only writes on destruct.
    # See https://github.com/acts-project/acts/issues/881
    del s

    mat_file = tmp_path / "material-map-volume.json"

    assert mat_file.exists()
    assert mat_file.stat().st_size > 10

    with mat_file.open() as fh:
        assert json.load(fh)

    assert map_file.exists()
    assert_entries(map_file, "material-tracks", 200)
    assert_root_hash(map_file.name, map_file)

    val_file = tmp_path / "propagation-volume-material.root"
    assert not val_file.exists()

    # test the validation as well

    # we need to destroy the ODD to reload with material
    # del trackingGeometry
    # del detector

    detector, trackingGeometry, decorators = getOpenDataDetector(
        mdecorator=acts.IMaterialDecorator.fromFile(mat_file)
    )

    from material_validation import runMaterialValidation

    s = Sequencer(events=10, numThreads=1)

    field = acts.NullBField()

    runMaterialValidation(
        trackingGeometry,
        decorators,
        field,
        outputDir=str(tmp_path),
        outputName="propagation-volume-material",
        s=s,
    )

    s.run()

    assert val_file.exists()
    assert_root_hash(val_file.name, val_file)


@pytest.mark.parametrize(
    "geoFactory,nobj",
    [
        (GenericDetector.create, 450),
        pytest.param(
            getOpenDataDetector,
            540,
            marks=pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up"),
        ),
        (functools.partial(AlignedDetector.create, iovSize=1), 450),
    ],
)
def test_geometry_example(geoFactory, nobj, tmp_path):
    detector, trackingGeometry, decorators = geoFactory()

    from geometry import runGeometry

    json_dir = tmp_path / "json"
    csv_dir = tmp_path / "csv"
    obj_dir = tmp_path / "obj"

    for d in (json_dir, csv_dir, obj_dir):
        d.mkdir()

    events = 5

    kwargs = dict(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        events=events,
        outputDir=str(tmp_path),
    )

    runGeometry(outputJson=True, **kwargs)
    runGeometry(outputJson=False, **kwargs)

    assert len(list(obj_dir.iterdir())) == nobj
    assert all(f.stat().st_size > 200 for f in obj_dir.iterdir())

    assert len(list(csv_dir.iterdir())) == 3 * events
    assert all(f.stat().st_size > 200 for f in csv_dir.iterdir())

    detector_files = [csv_dir / f"event{i:>09}-detectors.csv" for i in range(events)]
    for detector_file in detector_files:
        assert detector_file.exists()
        assert detector_file.stat().st_size > 200

    contents = [f.read_text() for f in detector_files]
    ref = contents[0]
    for c in contents[1:]:
        if isinstance(detector, AlignedDetector):
            assert c != ref, "Detector writeout is expected to be different"
        else:
            assert c == ref, "Detector writeout is expected to be identical"

    if not isinstance(detector, AlignedDetector):
        for f in [json_dir / f"event{i:>09}-detector.json" for i in range(events)]:
            assert detector_file.exists()
            with f.open() as fh:
                data = json.load(fh)
                assert data
        material_file = tmp_path / "geometry-map.json"
        assert material_file.exists()
        assert material_file.stat().st_size > 200


def test_digitization_example(trk_geo, tmp_path, assert_root_hash):
    from digitization import configureDigitization

    s = Sequencer(events=10, numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "measurements.root"

    assert not root_file.exists()
    assert not csv_dir.exists()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    configureDigitization(trk_geo, field, outputDir=tmp_path, s=s)

    s.run()

    assert root_file.exists()
    assert csv_dir.exists()

    assert len(list(csv_dir.iterdir())) == 3 * s.config.events
    assert all(f.stat().st_size > 50 for f in csv_dir.iterdir())
    for tn, nev in (
        (8, 407),
        (9, 0),
        (12, 11),
        (13, 375),
        (14, 2),
        (16, 25),
        (17, 146),
        (18, 9),
    ):
        assert_entries(root_file, f"vol{tn}", nev)

    assert_root_hash(root_file.name, root_file)


def test_digitization_example_input(trk_geo, tmp_path, assert_root_hash):
    from particle_gun import runParticleGun
    from digitization import configureDigitization

    ptcl_dir = tmp_path / "ptcl"
    ptcl_dir.mkdir()
    pgs = Sequencer(events=20, numThreads=-1)
    runParticleGun(str(ptcl_dir), s=pgs)
    pgs.run()

    s = Sequencer(numThreads=-1)

    csv_dir = tmp_path / "csv"
    root_file = tmp_path / "measurements.root"

    assert not root_file.exists()
    assert not csv_dir.exists()

    assert_root_hash(
        "particles.root",
        ptcl_dir / "particles.root",
    )

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    configureDigitization(
        trk_geo,
        field,
        outputDir=tmp_path,
        particlesInput=ptcl_dir / "particles.root",
        s=s,
    )

    s.run()

    assert root_file.exists()
    assert csv_dir.exists()

    assert len(list(csv_dir.iterdir())) == 3 * pgs.config.events
    assert all(f.stat().st_size > 50 for f in csv_dir.iterdir())
    for tn, nev in (
        (7, 0),
        (8, 193),
        (9, 0),
        (12, 1),
        (13, 183),
        (14, 6),
        (16, 3),
        (17, 76),
        (18, 10),
    ):
        assert_entries(root_file, f"vol{tn}", nev)
    assert_root_hash(root_file.name, root_file)


def test_digitization_config_example(trk_geo, tmp_path):
    from digitization_config import runDigitizationConfig

    out_file = tmp_path / "output.json"
    assert not out_file.exists()

    input = (
        Path(__file__).parent
        / "../../../Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
    )
    assert input.exists(), input.resolve()

    runDigitizationConfig(trk_geo, input=input, output=out_file)

    assert out_file.exists()

    with out_file.open() as fh:
        data = json.load(fh)
    assert len(data.keys()) == 2
    assert data["acts-geometry-hierarchy-map"]["format-version"] == 0
    assert (
        data["acts-geometry-hierarchy-map"]["value-identifier"]
        == "digitization-configuration"
    )
    assert len(data["entries"]) == 27


def test_ckf_tracks_example_full_seeding(tmp_path, assert_root_hash):
    csv = tmp_path / "csv"

    assert not csv.exists()

    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        (
            "performance_ckf.root",
            None,
            None,
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_tracks",
            368,
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_particles",
            80,
        ),
        (
            "trackstates_ckf.root",
            "trackstates",
            368,
        ),
        (
            "tracksummary_ckf.root",
            "tracksummary",
            10,
        ),
    ]

    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            Path(__file__).parent.parent.parent.parent
            / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            Path(__file__).parent.parent.parent.parent
            / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)
            assert_root_hash(rf, rp)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 300 for f in csv.iterdir()])


def test_ckf_tracks_example_truth_estimate(tmp_path, assert_root_hash):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        ("performance_ckf.root", None, None),
        (
            "performance_seeding_trees.root",
            "track_finder_tracks",
            80,
        ),
        (
            "performance_seeding_trees.root",
            "track_finder_particles",
            80,
        ),
        (
            "trackstates_ckf.root",
            "trackstates",
            80,
        ),
        (
            "tracksummary_ckf.root",
            "tracksummary",
            10,
        ),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            Path(__file__).parent.parent.parent.parent
            / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            Path(__file__).parent.parent.parent.parent
            / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=True,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)
            assert_root_hash(rf, rp)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 100 for f in csv.iterdir()])


def test_ckf_tracks_example_truth_smeared(tmp_path, assert_root_hash):
    # the example as written is only compatible with the generic detector
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))
    events = 10
    s = Sequencer(events=events, numThreads=1)  # Digitization is not thread-safe

    root_files = [
        ("performance_ckf.root", None, None),
        (
            "trackstates_ckf.root",
            "trackstates",
            80,
        ),
        (
            "tracksummary_ckf.root",
            "tracksummary",
            10,
        ),
    ]

    csv = tmp_path / "csv"

    assert not csv.exists()
    for rf, _, _ in root_files:
        assert not (tmp_path / rf).exists()

    from ckf_tracks import runCKFTracks

    runCKFTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=Path(
            Path(__file__).parent.parent.parent.parent
            / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        ),
        digiConfigFile=Path(
            Path(__file__).parent.parent.parent.parent
            / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        ),
        outputCsv=True,
        outputDir=tmp_path,
        truthSmearedSeeded=True,
        truthEstimatedSeeded=False,
        s=s,
    )
    s.run()

    del s  # files are closed in destructors, not great

    assert csv.exists()
    for rf, tn, nume in root_files:
        rp = tmp_path / rf
        assert rp.exists()
        if tn is not None and nume is not None:
            assert_entries(rp, tn, nume)
            assert_root_hash(rf, rp)

    assert len([f for f in csv.iterdir() if f.name.endswith("CKFtracks.csv")]) == events
    assert all([f.stat().st_size > 300 for f in csv.iterdir()])


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep not set up")
@pytest.mark.slow
# @pytest.mark.filterwarnings("ignore::UserWarning")
def test_vertex_fitting(tmp_path):
    detector, trackingGeometry, decorators = getOpenDataDetector()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    from vertex_fitting import runVertexFitting, VertexFinder

    s = Sequencer(events=100)

    runVertexFitting(
        field,
        vertexFinder=VertexFinder.Truth,
        outputDir=Path.cwd(),
        s=s,
    )

    alg = AssertCollectionExistsAlg(["fittedVertices"], name="check_alg")
    s.addAlgorithm(alg)

    s.run()
    assert alg.events_seen == s.config.events


import itertools


@pytest.mark.parametrize(
    "finder,inputTracks,entries",
    [
        ("Truth", False, 100),
        # ("Truth", True, 0), # this combination seems to be not working
        ("Iterative", False, 100),
        ("Iterative", True, 100),
        ("AMVF", False, 100),
        ("AMVF", True, 100),
    ],
)
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_vertex_fitting_reading(
    tmp_path, ptcl_gun, rng, finder, inputTracks, entries, assert_root_hash
):

    ptcl_file = tmp_path / "particles.root"

    detector, trackingGeometry, decorators = GenericDetector.create()
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    from vertex_fitting import runVertexFitting, VertexFinder

    inputTrackSummary = None
    if inputTracks:
        from truth_tracking import runTruthTracking

        s2 = Sequencer(numThreads=1, events=100)
        runTruthTracking(
            trackingGeometry,
            field,
            digiConfigFile=Path(
                "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
            ),
            outputDir=tmp_path,
            s=s2,
        )
        s2.run()
        del s2
        inputTrackSummary = tmp_path / "tracksummary_fitter.root"
        assert inputTrackSummary.exists()
        assert ptcl_file.exists()
    else:
        s0 = Sequencer(events=100, numThreads=1)
        evGen = ptcl_gun(s0)
        s0.addWriter(
            RootParticleWriter(
                level=acts.logging.INFO,
                inputParticles=evGen.config.outputParticles,
                filePath=str(ptcl_file),
            )
        )
        s0.run()
        del s0

        assert ptcl_file.exists()

    finder = VertexFinder[finder]

    s3 = Sequencer(numThreads=1)

    runVertexFitting(
        field,
        inputParticlePath=ptcl_file,
        inputTrackSummary=inputTrackSummary,
        outputDir=tmp_path,
        vertexFinder=finder,
        s=s3,
    )

    alg = AssertCollectionExistsAlg(["fittedVertices"], name="check_alg")
    s3.addAlgorithm(alg)

    s3.run()

    vertexing_file = tmp_path / "performance_vertexing.root"
    assert vertexing_file.exists()

    assert_entries(vertexing_file, "vertexing", entries)
    assert_root_hash(vertexing_file.name, vertexing_file)
