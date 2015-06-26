import os,sys

onlyHasT=True
removeClones=True
#To do loops on a variable for a gaudijob
#TolTySlope = int(os.environ['TOLTYSLOPE'])
#TolTyOffset = int(os.environ['TOLTYOFFSET'])

#Number of Events to process
nEvts = 1000

Eta25Cut = True
# Activate the profiler for time improvements
CallGrind= False

#Produce all the Efficiency plots 
Plotta = True

#Configure HitManager to fix the Size for instance, or...add some smearing
#ConfigureManager = False
#fracPos = 0.125
#managerPlots = False
#FixedSize = 3


#To tune the Manager (Only Local Brunel v47r4 & Brunel v47r5)


#Seeding Versions to run
New      = True
New_PrSeeding = False

Options = "Debug_1000"

#######################################################################
# General imports
#######################################################################
from GaudiKernel.ProcessJobOptions import importOptions
from Gaudi.Configuration import *
from Configurables import Brunel , LHCbApp , DDDBConf
from Configurables import CondDB
CondDB().Upgrade     = True

#######################################################################
# Brunel settings
#######################################################################

# Brunel Settings

Brunel().EvtMax     = nEvts
Brunel().PrintFreq  = nEvts/10
Brunel().DataType   = "Upgrade"
Brunel().WithMC     = True 
Brunel().Simulation = True
Brunel().OutputType = "NONE"
#Bsphiphi+NOW-1000ev-Extended_NewFTDet.digi
EventSelector().Input = ["DATAFILE='Bsphiphi+NOW__Digi_0to1_ADCThreshold_3_ClusterMinWidth_1_ClusterMaxWidth_4_-1000ev-Extended.digi'"]
from Configurables import TrackSys
TrackSys().TrackTypes= [ "Seeding"]#,"Forward"]#,"Downstream"]?
from Configurables import RecSysConf
Brunel().RecoSequence = ["Decoding","Tr"]
# More options: "Vertex","RICH","CALO","MUON","PROTO","SUMMARY" # commented out while waiting for a fix
Brunel().Detectors = [ 'FT', 'Magnet' ]
CondDB().AllLocalTagsByDataType = ["VP_UVP+RICH_2019+UT_UUT", "FT_StereoAngle5"]
CondDB().LoadCALIBDB = 'HLT1' 


#, 'UT' ] Load Detectors
if (Plotta):
   HistogramPersistencySvc().OutputFile = "SciFiPlots_"+Options+".root"
   
   
#######################################################################
   # File / tag settings
#######################################################################
   
from Configurables import CondDB
CondDB().Upgrade     = True

LHCbApp.DDDBtag     = "dddb-20131025"
LHCbApp().CondDBtag = "sim-20130830-vc-md100"
CondDB().AllLocalTagsByDataType = ["VP_UVP+RICH_2019+UT_UUT", "FT_StereoAngle5"]
Brunel().InputType  = "DIGI" 


#######################################################################
# Misc settings
#######################################################################

from Configurables import RecMoniConf
RecMoniConf().MoniSequence = [ ]

from Configurables import L0Conf
L0Conf().EnsureKnownTCK=False

#######################################################################
# Define sequences and run
#######################################################################

def doIt():
   #------------------------------
   #Configure PrChecker
   #------------------------------
   
   from Configurables import GaudiSequencer
   from Configurables import PrChecker
   from Configurables import IdealStateCreator
   from Configurables import PrLHCbID2MCParticle
   # Produce tuple for hit and cluster study (Extra package PrClusterResidual need)                                    
   from Configurables import PrCheatedSciFiTracking
   cheated= PrCheatedSciFiTracking()
   cheated.MinXHits = 5
   cheated.MinStereoHits = 5
   cheated.MinTotHits = 9
   cheated.SpecialCase = True
   cheated.SpecialCase2 = False #6+3
   cheated.NoDoubleCounting = True
   cheated.DecodeData = True
   cheated.CheckReco= True
   GaudiSequencer("MCLinksUnpackSeq").Members = []
   GaudiSequencer("TrSeedingSeq").Members=[]
   GaudiSequencer("TrSeedingSeq").Members += [ cheated]

   from Configurables import PrClustersResidual
   cluster = PrClustersResidual()
   cluster.DoClusterResidual =True
   cluster.DoTrackStudy = True
   cluster.DecodeData = False
   #cluster.OutputLevel = 
   cluster.OnlyHasT = onlyHasT
   cluster.RemoveClones = removeClones
   from Configurables import PrChecker
   prChecker = PrChecker()
   GaudiSequencer("CheckPatSeq").Members= [  ]
   GaudiSequencer("CheckPatSeq").Members += [ cluster ]
   from Configurables import IdealStateCreator
   if (Plotta):
      prChecker.WriteTTrackHistos = True
      prChecker.Eta25Cut = Eta25Cut
      prChecker.UseElectrons = False
      prChecker.TriggerNumbers = True
   GaudiSequencer("CheckPatSeq").Members +=[prChecker]
   from Configurables import MCParticle2MCHitAlg, IdealStateCreator, PrPlotFTHits
   # Define the algorithms
   FTAssoc = MCParticle2MCHitAlg( "MCP2FTMCHitAlg", 
                                  MCHitPath = "MC/FT/Hits", 
                                  OutputData = "/Event/MC/Particles2MCFTHits" )
   
   # tell the Data On Demand Service about them
   DataOnDemandSvc().AlgMap[ "/Event/Link/MC/Particles2MCFTHits" ]    = FTAssoc
   DataOnDemandSvc().NodeMap[ "/Event/Link" ]    = "DataObject"
   DataOnDemandSvc().NodeMap[ "/Event/Link/MC" ] = "DataObject"

   
   
   

   
appendPostConfigAction(doIt)
NTupleSvc().Output = ["FILE1 DATAFILE='SciFi-Tuple-"+Options+".root' TYP='ROOT' OPT='NEW'"]



#---------------------------------
# Callgrind configuration
#---------------------------------

from Configurables import CallgrindProfile
def addProfile():
      p = CallgrindProfile()
      p.StartFromEventN = 20
      p.StopAtEventN = 80
      p.DumpAtEventN= 80
      p.DumpName = "PRTEST"
      GaudiSequencer("InitBrunelSeq").Members.insert(0,p )



# Now add the profiling algorithm to the sequence
if (CallGrind):
  appendPostConfigAction(addProfile)

