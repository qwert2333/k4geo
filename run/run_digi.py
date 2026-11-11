import os

from Gaudi.Configuration import *

# Loading the input SIM file, defining output file
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = "ALFA_o1_v00-onlyCALO_gun_3to5_100evt.root"
io_svc.Output = "Digi_ALFA_o1_v00-onlyCALO_gun_3to5_100evt.root"


################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# ALFA calorimeter-only geometry: Grainita ECAL barrel + endcaps, solenoid, HCal barrel + endcap, magnetic field.
# The truth tracking below needs the ECAL LayeredCalorimeterData extensions and the field from this geometry.
compact_file = "../FCCee/ALFA/compact/ALFA_o1_v00/ALFA_o1_v00-onlyCALO.xml"
geoservice.detectors = [compact_file]
geoservice.OutputLevel = INFO


################ Calorimeter digitization
# SiPM emulation: one FastSiPMDigiWithEdep instance per calorimeter collection.
# Output collections: <prefix>DigiHit, <prefix>DigiHit_SimHit_link, <prefix>DigiHit_MCParticle_link
from Configurables import FastSiPMDigiWithEdep


def make_sipm_digi(name, input_collection, readout, output_prefix, params):
    digi = FastSiPMDigiWithEdep(name)
    digi.InputSimCaloHitCollection = input_collection
    digi.OutputCaloHitCollection = output_prefix + "DigiHit"
    digi.OutputCaloSimLinkCollection = output_prefix + "DigiHit_SimHit_link"
    digi.OutputCaloMCPLinkCollection = output_prefix + "DigiHit_MCParticle_link"
    digi.ReadOutName = readout
    for key, value in params.items():
        setattr(digi, key, value)
    digi.OutputLevel = INFO
    return digi


# Super simplified digitization model: only keep scintillation p.e. (light yield) stochastics and 10 p.e. threshold. 
lightyield_only_params = dict(
    UseDigi=True,
    SiPMGainMean=1.0,
    SiPMGainSigma=0.0,
    ADCGainRatio=1.0,
    Pedestal=0.0,
    PedestalResolution=0.0,
    ADC=1000000000,
    ADCSwitch=1000000000,
    ReadoutThreshold=10,  # p.e.
)

# Grainita ECAL: 10 p.e./MeV
ECAL_LIGHT_YIELD = 10.0  # p.e./MeV
HCAL_MIP_ENERGY = 0.47  # MeV
HCAL_LIGHT_YIELD = 128  # p.e./MeV = 60 p.e./MIP

grainita_params = dict(lightyield_only_params, LightYield=ECAL_LIGHT_YIELD)
hcal_params = dict(lightyield_only_params, LightYield=HCAL_LIGHT_YIELD)

ecal_barrel_digi = make_sipm_digi(
    "EcalBarrelSiPMDigi",
    "GrainitaCalorimeterHitsRaw",
    "GrainitaEcalBarrelRO",
    "GrainitaEcalBarrel",
    dict(grainita_params, OuterRadius=2645),  # unit mm, fiber readout radius (as outerRadius in the SIM steering)
)
# ECAL endcap: the transport time in FastSiPMDigiWithEdep uses |OuterRadius - R| (barrel geometry);
# OuterRadius is left at its default here.
ecal_endcap_digi = make_sipm_digi(
    "EcalEndcapSiPMDigi",
    "GrainitaEndcapCalorimeterHitsRaw",
    "GrainitaEcalEndcapRO",
    "GrainitaEcalEndcap",
    grainita_params,
)
hcal_barrel_digi = make_sipm_digi(
    "HCalBarrelSiPMDigi",
    "HCalBarrelCollection",
    "HCalBarrelCollection",
    "HCalBarrel",
    hcal_params,
)
hcal_endcap_digi = make_sipm_digi(
    "HCalEndcapSiPMDigi",
    "HCalEndcapCollection",
    "HCalEndcapCollection",
    "HCalEndcap",
    hcal_params,
)
digi_algs = [ecal_barrel_digi, ecal_endcap_digi, hcal_barrel_digi, hcal_endcap_digi]


################ Truth tracks
# Truth tracks from the generated particles (TracksFromGenParticles, k4RecTracker) with the track state at the
# ECAL inner face, from the LayeredCalorimeterData extensions of the Grainita ECAL barrel and endcaps.
# Calorimeter-only geometry: no SimTrackerHits, the helix starts at the production vertex.
# Needs a k4RecTracker build with the UseSimTrackerHits property (e.g. k4RecTracker_main).
from Configurables import TracksFromGenParticles
tracksFromGenParticles = TracksFromGenParticles(
    "TracksFromGenParticles",
    InputGenParticles=["MCParticles"],
    InputSimTrackerHits=[],
    UseSimTrackerHits=False,
    OutputTracks=["TracksFromGenParticles"],
    OutputMCRecoTrackParticleAssociation=["TracksFromGenParticlesAssociation"],
    ExtrapolateToECal=True,
    KeepOnlyBestExtrapolation=False,
    OutputLevel=INFO,
)

# RNG for sipm emulation (TODO harmonize RNG with other modules)
from Configurables import HepRndm__Engine_CLHEP__RanluxEngine_ as RndmEngine
rndmEngine = RndmEngine('RndmGenSvc.Engine',
  SetSingleton = True,
  Seeds = [ int(os.environ.get("DIGI_SEED", "1234567")) ] # default seed is 1234567
)


from Configurables import RndmGenSvc
rndmGenSvc = RndmGenSvc("RndmGenSvc",
  Engine = rndmEngine.name()
)

################ Output
io_svc.outputCommands = [
  "drop *",
  "keep EventHeader",
  "keep MCParticles",
  "keep GrainitaEcalBarrelDigiHit",
  "keep GrainitaEcalEndcapDigiHit",
  "keep HCalBarrelDigiHit",
  "keep HCalEndcapDigiHit",
  "keep TracksFromGenParticles",
  "keep TracksFromGenParticlesAssociation",
  "keep GrainitaEcalBarrelDigiHit_MCParticle_link",
  "keep GrainitaEcalEndcapDigiHit_MCParticle_link",
  "keep HCalBarrelDigiHit_MCParticle_link",
  "keep HCalEndcapDigiHit_MCParticle_link",
]

# Profiling
from Configurables import AuditorSvc, ChronoAuditor, UniqueIDGenSvc
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

from k4FWCore import ApplicationMgr
application_mgr = ApplicationMgr(
    TopAlg = digi_algs + [tracksFromGenParticles],
    EvtSel = 'NONE',
    EvtMax = -1,
    ExtSvc = [
        EventDataSvc("EventDataSvc"),
        geoservice,
        audsvc,
        UniqueIDGenSvc("uidSvc"),
        rndmEngine,
        rndmGenSvc
    ],
    StopOnSignal = True,
)

for algo in application_mgr.TopAlg:
    algo.AuditExecute = True
