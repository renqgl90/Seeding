import os,sys

#TolTySlope = int(os.environ['TOLTYSLOPE'])
#TolTyOffset = int(os.environ['TOLTYOFFSET'])

#Number of Events to process
nEvts = 1000
Eta25Cut = True

# Activate the profiler
CallGrind= False

Plotta = True #Produce all the Efficiency plots
ConfigureManager = True
fracPos = 0.125
#managerPlots = False
FixedSize = 3

#To tune the Manager (Only Local Brunel v47r4 & Brunel v47r5)
#Seeding Versions to run

New      = True
Options = "HybridSeeding_1000evts_Test_RecoverTheEff_Case0"

#######################################################################
# General imports
#######################################################################
from GaudiKernel.ProcessJobOptions import importOptions
from Gaudi.Configuration import *
from Configurables import Brunel , LHCbApp , DDDBConf


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
EventSelector().Input = ["DATAFILE='./Bsphiphi+NOW__Digi_0to1_ADCThreshold_3_ClusterMinWidth_1_ClusterMaxWidth_4_-1000ev-Extended.digi'"]#Bsphiphi+NOW__Digi_0to1_ADCThreshold_3_ClusterMinWidth_1_ClusterMaxWidth_4_-1000ev-Extended.digi'"]
from Configurables import TrackSys
#TrackSys().TrackTypes= [ "Velo","Forward","Seeding","Matching","Downstream"]
TrackSys().TrackTypes= [ "Seeding"]

from Configurables import RecSysConf
Brunel().RecoSequence = ["Decoding","Tr"]
# More options: "Vertex","RICH","CALO","MUON","PROTO","SUMMARY" # commented out while waiting for a fix
Brunel().Detectors = ['FT', 'Magnet']

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
   GaudiSequencer("CheckPatSeq").Members += []
   prChecker = PrChecker()
   from Configurables import IdealStateCreator

   if (Plotta):
      prChecker.WriteTTrackHistos = 2
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
   
   #---------------------------------
   #Configure the HitManager
   #---------------------------------
   #if (ConfigureManager):
#      from Configurables import PrFTHitManager
#      manager = PrFTHitManager("PrFTHitManager")
      
#      manager.fracPosOffset = fracPos
#      manager.doTuple = manager
#      manager.HackSize1 = False
#      manager.FixError = False
#      manager.SizeFix = FixedSize

   
   #---------------------------------
   #Configure the Seeding Tracking
   #---------------------------------
   seedingSeq = GaudiSequencer("TrSeedingSeq")
   GaudiSequencer("TrBestSeq").Members = []
   GaudiSequencer("TrSeedingSeq").Members = []
   #if you do truthmatching in the pat reco
   GaudiSequencer("MCLinksUnpackSeq").Members =[]
   GaudiSequencer("RecoTrSeq").Members += [ seedingSeq ]
   from Configurables import PrHybridSeeding
   seeding = PrHybridSeeding()
   #seeding.OutputLevel = DEBUG #Uncomment this line for debug      
   seeding.InputName ="" #Standalone seeding Put "Forward to get forward imput"
   seeding.MaxNHits = 12 #Force algorithm to find 12 Hits track when > 12
   seeding.DecodeData = True # Switch it off if Runned after Forward
   seeding.XOnly = False
   #N Cases
   seeding.NCases = 3
   seeding.MinXPlanes = 4
   #Clones Kill
   seeding.RemoveClonesX = True
   seeding.RemoveClones = True
   seeding.minNCommonUV = 7 #>=7
   seeding.minCommonX = [2,2,2] # N Common X (Remove Clones [Case0,Case1,Case2])
   seeding.RemoveClonesUpDown = False # should Speed Up
   seeding.minNCommonUVUpDown = 2 #>=2
   
   #Flag Hits
   seeding.FlagHits = True
   seeding.SizeToFlag = 12 # >=size
   seeding.RemoveFlagged = True # Case 1 and 2 will not use flagged Hits
   #If Flag Size = 11
   seeding.Flag_MaxChi2 = 0.3 # track.chi2(hit) = 0.3
   seeding.Flag_MaxX0  = 200  # <200
   #dRatio Business
   seeding.UseCubicCorrection = True # dRatio correction in the fit
   seeding.dRatio = -0.000262
   seeding.UseCorrPosition = True # dRatio( x,y)
   seeding.UseCorrSlopes = False #  dRatio( bx, by)
   seeding.CConst = 2.458e8 # BackwardProjection ( To be used somewhere in the algo )
   #Recover Track ( to be fully implemented properly)
   seeding.RecoverTrack = False
   seeding.ChiDoFRecover = -1.0
   #Case 0,1,2 parameters 1st Last search

   #XZ-Search
   # 1st - Last Layer
   seeding.L0_AlphaCorr = [120.64, 510.64, 730.64 ]
   seeding.L0_tolHp     = [280.0 , 540.0 , 1080.0 ]
   # ParabolaSeedHits
   seeding.x0Corr        =     [0.002152, 0.001534, 0.001534]
   
   seeding.X0SlopeChange =     [500.    ,    500. ,    500. ]
   seeding.x0Cut        =      [4000.   ,    4000.,   4000. ]
   seeding.TolAtX0Cut   =      [12.0    ,     8.0 ,    8.0  ]
   seeding.ToleranceX0Up =     [ 0.75   ,     0.75,    0.75 ]

   seeding.X0SlopeChangeDown = [ 1500.  ,    2000.,    1500. ]
   seeding.TolAtX0CutOpp     = [3.0     ,      2.0,     2.0  ]
   seeding.ToleranceX0Down   = [ 0.75   ,      0.75,    0.75 ]
   
   seeding.TolXRemaining     = [1.0     ,      1.0 ,     1.0]
   seeding.maxParabolaSeedHits = 12 # Hits are sorted by the distance from
   #  x0 = first-last projection to z=0 ; then hits sorted by xParabola - (x0+tx_1stLast*z_PlaneParabola + x0Corr*x0)  and we keep the first maxParabolaSeedHits list
   seeding.maxChi2HitsX     = [5.5      , 5.5     , 5.5   ]
   seeding.maxChi2DoFX      = [4.0      , 5.0     , 6.0   ]
   
   #Add Stereo Part
   seeding.DoAsymm = True
   seeding.TriangleFix = True
   seeding.TriangleFix2ndOrder = True
   seeding.yMin = -1.0
   seeding.yMin_TrFix = -2.0
   seeding.yMax =    2700
   seeding.yMax_TrFix = 30.0
   seeding.RemoveHole = True
   #Hough Cluster Size settings Hits are sorted by y/z where y is computed from the XZ-segment processed
   seeding.TolTyOffset = [ 0.002  , 0.002 , 0.0035  ]
   seeding.TolTySlope  = [ 0.0    , 0.0   , 0.015   ]

   #Once you find The UV hits you check the Line Y ?
   seeding.UseLineY = True
   
   #9 and 10 hits
   seeding.Chi2LowLine =           [ 5.0 , 6.0  , 7.0  ] #Chi2PerDoF when 10 Hits or 9
   seeding.maxChi2Hits_less11Hit = [ 2.5 , 2.5  , 2.5  ]
   seeding.maxYatZeroLow =         [ 50., 50.   , 50.  ]
   seeding.maxYatzRefLow =         [500., 500.  , 500. ]
   #11 and 12 Hits
   seeding.Chi2HighLine =          [30.0, 50.0 , 80.0] #Chi2PerDoF when 11 or 12 hits Using LineChi2DoF + XZChi2DoF
   seeding.maxChi2Hits_11and12Hit =[ 5.5, 5.5, 5.5    ] 
   
   #Make the Full Fit
   seeding.maxChi2PerDoF          = [4.0, 6.0, 7.0]
      
   seeding.RecoverTrack = False #Get only tracks with 6 UV or 6 X
   seeding.ChiDoFRecover = -1.0 #Change it only if RecoverTrack = True
   
   #TruthMatching Settings (Comment all of them if you want to run the normal seeding 
   if(ConfigureManager):
      from Configurables import PrFTHitManager
      #seeding.addTool(manager)
      #seeding.PrFTHitManager.SizeFix = FixedSize
      #seeding.PrFTHitManager.FixError = True
   seedingSeq.Members = [seeding]

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
# to run it gaudirun.py -T --profilerName=valgrindcachegrind --profilerExtraOptions="__instr-atstart=no -v __smc-check=all-non-file __dump-instr=yes __trace-jump=yes" BrunelStep_DevelopSeeding.py`
# kcachegrind valgrindcallgrind.output.log`
# 

if (CallGrind):
  appendPostConfigAction(addProfile)

