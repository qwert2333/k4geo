import os

from Gaudi.Configuration import *

# Loading the input SIM file, defining output file
from k4FWCore import IOSvc
from Configurables import EventDataSvc
io_svc = IOSvc("IOSvc")
io_svc.Input = "AlfaCalSimulation.root"
io_svc.Output = "ALFA_CaloDigi_pi-_10GeV.root"

################## Simulation setup
# Detector geometry
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
compact_file = "FCCee/ALFA/compact/ALFA_o1_v00/ALFA_o1_v00.xml"
geoservice.detectors = [compact_file]
geoservice.OutputLevel = INFO


from Configurables import TracksFromGenParticles
tracksFromGenParticles = TracksFromGenParticles("CreateTracksFromGenParticles",
                            InputGenParticles=["MCParticles"],
                            InputSimTrackerHits=["OTBarCollection"],
                            OutputTracks=["TracksFromGenParticles"],
                            OutputMCRecoTrackParticleAssociation=["TracksFromGenParticlesAssociation"],
                            ExtrapolateToECal=False,
                            KeepOnlyBestExtrapolation=False,
                            TrackerIDs=[67],
                            OutputLevel=INFO)

################ Dual-readout calorimeter
# SiPM emulation
from Configurables import SimpleSiPMDigiAlg
simpledigi = SimpleSiPMDigiAlg("simpledigi_ecalbarrel")
simpledigi.InputSimCaloHitCollection = "GrainitaCalorimeterHits"
simpledigi.OutputCaloHitCollection = "GrainitaEcalBarrelDigiHit"
simpledigi.OutputCaloSimLinkCollection = "GrainitaEcalBarrelDigiHit_SimHit_link"
simpledigi.OutputCaloMCPLinkCollection = "GrainitaEcalBarrelDigiHit_MCParticle_link"
simpledigi.ReadOutName = "GrainitaEcalBarrelRO"
simpledigi.UseDigi = True
simpledigi.TimeWindow = 25000  # unit ns
simpledigi.SiPMCT = 0.01
simpledigi.SiPMDCR = 90e-6  # unit GHz
simpledigi.SiPMPixel = 10000
simpledigi.LightYield = 10  # ph / MeV
simpledigi.WriteNtuple = False
simpledigi.OutFileName = "DigiTuple_EcalBarrel_pi-_30GeV.root"
simpledigi.OutputLevel = INFO


from Configurables import SimpleCalibAlg
simplecali = SimpleCalibAlg("SimpleCalibAlg")
simplecali.InputCaloHitCollection = "GrainitaEcalBarrelDigiHit"
simplecali.OutputCaloHitCollection = "GrainitaEcalBarrelCalibHit"
simplecali.ReadOutName = "GrainitaEcalBarrelRO"
simplecali.InputFormat = "Energy"
simplecali.CalibrationConstant = 1.
simplecali.ApplyRhoPitchCorrection = True
#simplecali.AttLength = 0.34   # unit in cm
simplecali.OutputLevel = INFO


simpledigi_hcalbarrel = SimpleSiPMDigiAlg("simpledigi_hcalbarrel")
simpledigi_hcalbarrel.InputSimCaloHitCollection = "HCalBarrelCollection"
simpledigi_hcalbarrel.OutputCaloHitCollection = "HcalBarrelDigiHit"
simpledigi_hcalbarrel.OutputCaloSimLinkCollection = "HcalBarrelDigiHit_SimHit_link"
simpledigi_hcalbarrel.OutputCaloMCPLinkCollection = "HcalBarrelDigiHit_MCParticle_link"
simpledigi_hcalbarrel.ReadOutName = "HCalBarrelCollection"
simpledigi_hcalbarrel.UseDigi = True
simpledigi_hcalbarrel.TimeWindow = 25000  # unit ns
simpledigi_hcalbarrel.SiPMCT = 0.01
simpledigi_hcalbarrel.SiPMDCR = 90e-6  # unit GHz
simpledigi_hcalbarrel.SiPMPixel = 10000
simpledigi_hcalbarrel.LightYield = 100  # ph / MeV
simpledigi_hcalbarrel.WriteNtuple = False
simpledigi_hcalbarrel.OutFileName = "DigiTuple_HcalBarrel_pi-_30GeV.root"
simpledigi_hcalbarrel.OutputLevel = INFO


simpledigi_hcalendcap = SimpleSiPMDigiAlg("simpledigi_hcalendcap")
simpledigi_hcalendcap.InputSimCaloHitCollection = "HCalEndcapCollection"
simpledigi_hcalendcap.OutputCaloHitCollection = "HcalEndcapDigiHit"
simpledigi_hcalendcap.OutputCaloSimLinkCollection = "HcalEndcapDigiHit_SimHit_link"
simpledigi_hcalendcap.OutputCaloMCPLinkCollection = "HcalEndcapDigiHit_MCParticle_link"
simpledigi_hcalendcap.ReadOutName = "HCalEndcapCollection"
simpledigi_hcalendcap.UseDigi = True
simpledigi_hcalendcap.TimeWindow = 25000  # unit ns
simpledigi_hcalendcap.SiPMCT = 0.01
simpledigi_hcalendcap.SiPMDCR = 90e-6  # unit GHz
simpledigi_hcalendcap.SiPMPixel = 10000
simpledigi_hcalendcap.LightYield = 100  # ph / MeV
simpledigi_hcalendcap.WriteNtuple = False
simpledigi_hcalendcap.OutFileName = "DigiTuple_HcalEndcap_pi-_30GeV.root"
simpledigi_hcalendcap.OutputLevel = INFO



# RNG for sipm emulation (TODO harmonize RNG with other modules)
from Configurables import HepRndm__Engine_CLHEP__RanluxEngine_ as RndmEngine
rndmEngine = RndmEngine('RndmGenSvc.Engine',
  SetSingleton = True,
  Seeds = [ int(os.environ.get("DIGI_SEED", "1234567")) ] # default seed is 1234567
)

#from Configurables import CreateTruthLinks
#createTruthLinks = CreateTruthLinks("CreateTruthLinks",
#    cell_hit_links=["GrainitaCaloSiPMreadoutDigiHit_link"],
#    clusters=["TopoClusterAll"],
#    mcparticles="MCParticles",
#    cell_mcparticle_links="CaloHitMCParticleLinks",
#    cluster_mcparticle_links="ClusterMCParticleLinks",
#    OutputLevel=INFO
#)

from Configurables import RndmGenSvc
rndmGenSvc = RndmGenSvc("RndmGenSvc",
  Engine = rndmEngine.name()
)

################ Output
io_svc.outputCommands = [
  "keep *",
]

# Profiling
from Configurables import AuditorSvc, ChronoAuditor, UniqueIDGenSvc
chra = ChronoAuditor()
audsvc = AuditorSvc()
audsvc.Auditors = [chra]

from k4FWCore import ApplicationMgr
application_mgr = ApplicationMgr(
    TopAlg = [
        tracksFromGenParticles,
        simpledigi, simplecali, simpledigi_hcalbarrel, simpledigi_hcalendcap
    ],
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
