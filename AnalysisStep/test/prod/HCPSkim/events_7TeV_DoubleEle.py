LEPTON_SETUP = 2011

IsMC = False
PD = "DoubleEle"
MCFILTER = ""
ELECORRTYPE = "Jan16ReReco"
APPLYELEREGRESSION = True
APPLYELECALIB = True
APPLYMUCORR = True

import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZAnalysis/AnalysisStep/test/"
execfile(PyFilePath + "analyzer.py")        
execfile(PyFilePath + "prod/json_2011.py")      

process.patElectronsWithRegression.rhoCollection = cms.InputTag('kt6PFJetsForIso:rho')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source.fileNames = cms.untracked.vstring(
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_1655.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_1810.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_1853.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_200.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_290.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_299.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_397.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_441.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_530.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_540.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_550.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_749.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_871.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_116.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_92.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_928.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_1262.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_134.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011A-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_1392.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_104.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_912.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_916.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_94.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_941.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_183.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_186.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_187.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_232.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_245.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_264.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_298.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_12.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_121.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_344.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_366.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_370.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_392.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_401.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_126.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_445.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_447.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_477.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_53.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_556.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_627.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_634.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_667.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_692.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_745.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_751.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_763.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_767.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_84.root',
                                'root://eoscms//eos/cms/store/cmst3/user/cmgtools/CMG/DoubleElectron/Run2011B-16Jan2012-v1/AOD/V5/PAT_CMG_V5_9_0/cmgTuple_881.root',
                                )


process.source.eventsToProcess = cms.untracked.VEventRange(
    '163334:221:129514273',
    '163334:499:286336207',
    '163659:463:344708580',
    '163758:151:113529826',
    '163795:34:30998576',
    '163817:174:155679852',
    '165633:303:394010457',
    '165970:236:275108397',
    '166408:724:917379387',
    '166438:79:78213037',
    '166438:768:862270386',
    '166512:281:337493970',
    '166554:333:395098004',
    '166922:85:112725318',
    '166950:1373:1491724484',
    '167281:386:480301165',
    '167282:44:44166176',
    '167807:750:966824024',
    '171050:459:591031316',
    '171106:127:141954801',
    '171369:150:160966858',
    '172163:128:191231387',
    '172208:75:66033190',
    '172401:7:3729470',
    '172620:242:218903169',
    '172799:11:10347106',
    '172802:125:107360878',
    '172819:220:298086610',
    '172868:689:933807102',
    '172868:1822:2299160459',
    '172949:840:1188043146',
    '172952:466:559839432',
    '172992:836:1153485608',
    '173243:12:16706390',
    '173657:85:65557571',
    '173659:270:389185367',
    '173692:2066:2722114329',
    '175906:190:227517585',
    '175921:220:297753357',
    '175921:349:495614354',
    '175974:9:7526662',
    '176201:182:261184429',
    '176201:354:562295642',
    '176207:206:256888239',
    '176304:300:418052877',
    '176309:224:257489763',
    '176309:950:1340034258',
    '176468:128:215855118',
    '176548:231:403771114',
    '176799:24:35688265',
    '176886:260:427567024',
    '176886:631:1057019814',
    '177074:384:588602439',
    '177074:582:931848091',
    '177139:183:290826062',
    '177222:227:339499459',
    '177318:169:270676815',
    '177449:68:58273256',
    '177782:99:72158025',
    '177790:168:222240677',
    '177790:527:657843813',
    '177875:133:148667118',
    '177875:289:419370375',
    '178100:236:326364918',
    '178116:428:695859609',
    '178116:437:709511403',
    '178162:10:10608364',
    '178421:86:87514902',
    '178421:973:1450980155',
    '178421:1087:1610336854',
    '178424:585:666626491',
    '178479:210:298608854',
    '178479:369:589085976',
    '178479:470:757111474',
    '178703:137:191352626',
    '178708:354:573962528',
    '178731:192:248562036',
    '178786:197:277942410',
    '178866:82:140063742',
    '178970:103:122998167',
    '178970:398:658883361',
    '179434:52:86225612',
    '179452:1056:1459855927',
    '179476:30:30532070',
    '179563:871:1409064222',
    '180076:46:79350642',
    '180076:271:456795917',
    '180250:28:45096064',
    '180250:326:591651181',
    '180250:496:905064541',
)
