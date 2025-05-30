#include "MassInv.h"

// Gaudi
#include "GaudiKernel/ITHistSvc.h"

// Key4HEP
#include "k4Interface/IGeoSvc.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDSegmentation/MultiSegmentation.h"

// EDM4HEP
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/Cluster.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/MCParticleCollection.h"

// ROOT
#include "TFile.h"
#include "TFitResult.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TSystem.h"

DECLARE_COMPONENT(MassInv)

MassInv::MassInv(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc), m_histSvc("THistSvc", "MassInv"), m_geoSvc("GeoSvc", "MassInv"),
      m_hEnergyPreAnyCorrections(nullptr), m_hEnergyPostAllCorrections(nullptr), m_hPileupEnergy(nullptr),
      m_hUpstreamEnergy(nullptr) {
  declareProperty("clusters", m_inClusters, "Input clusters (input)");
  declareProperty("correctedClusters", m_correctedClusters, "Corrected clusters (output)");
  declareProperty("particle", m_particle, "Generated single-particle event (input)");
  declareProperty("towerTool", m_towerTool, "Handle for the tower building tool");
}

StatusCode MassInv::initialize() {
  {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure())
      return sc;
  }

  int energyStart = 0;
  int energyEnd = 0;
  if (m_energy == 0) {
    energyStart = 0;
    energyEnd = 1000;
  } else {
    energyStart = 0.2 * m_energy;
    energyEnd = 1.4 * m_energy;
  }

  // create control histograms
  m_hEnergyPreAnyCorrections =
      new TH1F("energyPreAnyCorrections", "Energy of cluster before any correction", 3000, energyStart, energyEnd);
  if (m_histSvc->regHist("/rec/energyPreAnyCorrections", m_hEnergyPreAnyCorrections).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hEnergyPostAllCorrections =
      new TH1F("energyPostAllCorrections", "Energy of cluster after all corrections", 3000, energyStart, energyEnd);
  if (m_histSvc->regHist("/rec/energyPostAllCorrections", m_hEnergyPostAllCorrections).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hEnergyPostAllCorrectionsAndScaling =
      new TH1F("energyPostAllCorrectionsAndScaling", "Energy of cluster after all corrections and scaling", 3000,
               energyStart, energyEnd);
  if (m_histSvc->regHist("/rec/energyPostAllCorrectionsAndScaling", m_hEnergyPostAllCorrectionsAndScaling)
          .isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hPileupEnergy = new TH1F("pileupCorrectionEnergy", "Energy added to a cluster as a correction for correlated noise",
                             1000, -10, 10);
  if (m_histSvc->regHist("/rec/pileupCorrectionEnergy", m_hPileupEnergy).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hUpstreamEnergy = new TH1F("upstreamCorrectionEnergy",
                               "Energy added to a cluster as a correction for upstream material", 1000, -10, 10);
  if (m_histSvc->regHist("/rec/upstreamCorrectionEnergy", m_hUpstreamEnergy).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hDiffEta =
      new TH1F("diffEta", "#eta resolution", 10 * ceil(2 * m_etaMax / m_dEta), -m_etaMax / 10., m_etaMax / 10.);
  if (m_histSvc->regHist("/rec/diffEta", m_hDiffEta).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  for (uint i = 0; i < m_numLayers; i++) {
    m_hDiffEtaLayer.push_back(new TH1F(("diffEtaLayer" + std::to_string(i)).c_str(),
                                       ("#eta resolution for layer " + std::to_string(i)).c_str(),
                                       10 * ceil(2 * m_etaMax / m_dEta), -m_etaMax / 10., m_etaMax / 10.));
    if (m_histSvc->regHist("/rec/diffEta_layer" + std::to_string(i), m_hDiffEtaLayer.back()).isFailure()) {
      error() << "Couldn't register histogram" << endmsg;
      return StatusCode::FAILURE;
    }
    m_hDiffEtaHitLayer.push_back(new TH1F(("diffEtaHitLayer" + std::to_string(i)).c_str(),
                                          ("#eta hot distribution for layer " + std::to_string(i)).c_str(),
                                          10 * ceil(2 * m_etaMax / m_dEta), -m_etaMax / 10., m_etaMax / 10.));
    if (m_histSvc->regHist("/rec/diffEtaHit_layer" + std::to_string(i), m_hDiffEtaHitLayer.back()).isFailure()) {
      error() << "Couldn't register histogram" << endmsg;
      return StatusCode::FAILURE;
    }
  }
  m_hEta = new TH1F("eta", "#eta", 1000, -m_etaMax, m_etaMax);
  if (m_histSvc->regHist("/rec/eta", m_hEta).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hDiffPhi =
      new TH1F("diffPhi", "#varphi resolution", 10 * ceil(2 * m_phiMax / m_dPhi), -m_phiMax / 10., m_phiMax / 10.);
  if (m_histSvc->regHist("/rec/diffPhi", m_hDiffPhi).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hPhi = new TH1F("phi", "#varphi", 1000, -m_phiMax, m_phiMax);
  if (m_histSvc->regHist("/rec/phi", m_hPhi).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  double thetaMin = 2. * atan(exp(-static_cast<double>(m_etaMax)));
  double thetaMax = 2. * atan(exp(static_cast<double>(m_etaMax)));
  m_hDiffTheta = new TH1F("diffTheta", "#theta resolution", 10 * ceil((thetaMax - thetaMin) / 0.01), -0.25, 0.25);
  if (m_histSvc->regHist("/rec/diffTheta", m_hDiffTheta).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hDiffTheta2point =
      new TH1F("diffTheta2point", "#theta resolution", 10 * ceil((thetaMax - thetaMin) / 0.01), -0.25, 0.25);
  if (m_histSvc->regHist("/rec/diffTheta2point", m_hDiffTheta2point).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hDiffZ = new TH1F("diffZ", "z resolution", 1e4, -10, 10);
  if (m_histSvc->regHist("/rec/diffZ", m_hDiffZ).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hNumCells = new TH1F("numCells", "number of cells", 2000, -0.5, 1999.5);
  if (m_histSvc->regHist("/rec/numCells", m_hNumCells).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hDiPT = new TH1F("diPT", "transverse momentum of diparticles", 5000, 0, 500);
  if (m_histSvc->regHist("/rec/diPT", m_hDiPT).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hDiPTScaled = new TH1F("diPTScaled",
                           ("transverse momentum of diparticles with cluster energy scaled to " +
                            std::to_string(round(1. / m_response * 100) / 100))
                               .c_str(),
                           5000, 0, 500);
  if (m_histSvc->regHist("/rec/diPTScaled", m_hDiPTScaled).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInv = new TH1F("massInv", "invariant mass", 5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInv", m_hMassInv).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaled = new TH1F(
      "massInvScaled",
      ("invariant mass with cluster energy scaled to " + std::to_string(round(1. / m_response * 100) / 100)).c_str(),
      5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaled", m_hMassInvScaled).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaled100 = new TH1F("massInvScaled100",
                                 ("invariant mass for pT>100GeV with cluster energy scaled to " +
                                  std::to_string(round(1. / m_response * 100) / 100))
                                     .c_str(),
                                 5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaled100", m_hMassInvScaled100).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaled200 = new TH1F("massInvScaled200",
                                 ("invariant mass for pT>200GeV with cluster energy scaled to " +
                                  std::to_string(round(1. / m_response * 100) / 100))
                                     .c_str(),
                                 5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaled200", m_hMassInvScaled200).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaled300 = new TH1F("massInvScaled300",
                                 ("invariant mass for pT>300GeV with cluster energy scaled to " +
                                  std::to_string(round(1. / m_response * 100) / 100))
                                     .c_str(),
                                 5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaled300", m_hMassInvScaled300).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated = new TH1F(
      "massInvScaledIsolated",
      ("invariant mass with cluster energy scaled to " + std::to_string(round(1. / m_response * 100) / 100)).c_str(),
      5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated", m_hMassInvScaledIsolated).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated100 = new TH1F("massInvScaledIsolated100",
                                         ("invariant mass for pT>100GeV with cluster energy scaled to " +
                                          std::to_string(round(1. / m_response * 100) / 100))
                                             .c_str(),
                                         5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated100", m_hMassInvScaledIsolated100).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated200 = new TH1F("massInvScaledIsolated200",
                                         ("invariant mass for pT>200GeV with cluster energy scaled to " +
                                          std::to_string(round(1. / m_response * 100) / 100))
                                             .c_str(),
                                         5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated200", m_hMassInvScaledIsolated200).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated300 = new TH1F("massInvScaledIsolated300",
                                         ("invariant mass for pT>300GeV with cluster energy scaled to " +
                                          std::to_string(round(1. / m_response * 100) / 100))
                                             .c_str(),
                                         5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated300", m_hMassInvScaledIsolated300).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated2 = new TH1F(
      "massInvScaledIsolated2",
      ("invariant mass with cluster energy scaled to " + std::to_string(round(1. / m_response * 100) / 100)).c_str(),
      5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated2", m_hMassInvScaledIsolated2).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated2100 = new TH1F("massInvScaledIsolated2100",
                                          ("invariant mass for pT>100GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated2100", m_hMassInvScaledIsolated2100).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated2200 = new TH1F("massInvScaledIsolated2200",
                                          ("invariant mass for pT>200GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated2200", m_hMassInvScaledIsolated2200).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated2300 = new TH1F("massInvScaledIsolated2300",
                                          ("invariant mass for pT>300GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated2300", m_hMassInvScaledIsolated2300).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated3 = new TH1F(
      "massInvScaledIsolated3",
      ("invariant mass with cluster energy scaled to " + std::to_string(round(1. / m_response * 100) / 100)).c_str(),
      5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated3", m_hMassInvScaledIsolated3).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated3100 = new TH1F("massInvScaledIsolated3100",
                                          ("invariant mass for pT>100GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated3100", m_hMassInvScaledIsolated3100).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated3200 = new TH1F("massInvScaledIsolated3200",
                                          ("invariant mass for pT>200GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated3200", m_hMassInvScaledIsolated3200).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated3300 = new TH1F("massInvScaledIsolated3300",
                                          ("invariant mass for pT>300GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated3300", m_hMassInvScaledIsolated3300).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated4 = new TH1F(
      "massInvScaledIsolated4",
      ("invariant mass with cluster energy scaled to " + std::to_string(round(1. / m_response * 100) / 100)).c_str(),
      5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated4", m_hMassInvScaledIsolated4).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated4100 = new TH1F("massInvScaledIsolated4100",
                                          ("invariant mass for pT>100GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated4100", m_hMassInvScaledIsolated4100).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated4200 = new TH1F("massInvScaledIsolated4200",
                                          ("invariant mass for pT>200GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated4200", m_hMassInvScaledIsolated4200).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated4300 = new TH1F("massInvScaledIsolated4300",
                                          ("invariant mass for pT>300GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated4300", m_hMassInvScaledIsolated4300).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated5 = new TH1F(
      "massInvScaledIsolated5",
      ("invariant mass with cluster energy scaled to " + std::to_string(round(1. / m_response * 100) / 100)).c_str(),
      5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated5", m_hMassInvScaledIsolated5).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated5100 = new TH1F("massInvScaledIsolated5100",
                                          ("invariant mass for pT>100GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated5100", m_hMassInvScaledIsolated5100).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated5200 = new TH1F("massInvScaledIsolated5200",
                                          ("invariant mass for pT>200GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated5200", m_hMassInvScaledIsolated5200).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledIsolated5300 = new TH1F("massInvScaledIsolated5300",
                                          ("invariant mass for pT>300GeV with cluster energy scaled to " +
                                           std::to_string(round(1. / m_response * 100) / 100))
                                              .c_str(),
                                          5000, 0, 500);
  if (m_histSvc->regHist("/rec/massInvScaledIsolated5300", m_hMassInvScaledIsolated5300).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hMassInvScaledPt = new TH2F(
      "massInvPtScaled",
      ("invariant mass vs p_T with cluster energy scaled to " + std::to_string(round(1. / m_response * 100) / 100))
          .c_str(),
      5000, 0, 500, 5000, 0, 1000);
  if (m_histSvc->regHist("/rec/massInPtScaled", m_hMassInvScaledPt).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hHCalEnergy = new TH1F("HCalenergy", "Energy deposited in HCal behind EM clusters", 10000, 0, 100);
  if (m_histSvc->regHist("/rec/energyHCal", m_hHCalEnergy).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  m_hHCalTotalEnergy = new TH1F("HCalenergyTotal", "Total deposited energy in HCal", 10000, 0, 1000);
  if (m_histSvc->regHist("/rec/energyTotalHCal", m_hHCalTotalEnergy).isFailure()) {
    error() << "Couldn't register histogram" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_etaRecalcLayerWeights.size() < m_numLayers) {
    error() << "m_etaRecalcLayerWeights size is smaller than numLayers." << endmsg;
    return StatusCode::FAILURE;
  }
  for (uint iSys = 0; iSys < m_systemId.size(); iSys++) {
    // check if readouts exist
    if (m_geoSvc->getDetector()->readouts().find(m_readoutName[iSys]) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutName[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    // retrieve PhiEta segmentation
    m_segmentationPhiEta[m_systemId[iSys]] = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
        m_geoSvc->getDetector()->readout(m_readoutName[iSys]).segmentation().segmentation());
    m_segmentationMulti[m_systemId[iSys]] = dynamic_cast<dd4hep::DDSegmentation::MultiSegmentation*>(
        m_geoSvc->getDetector()->readout(m_readoutName[iSys]).segmentation().segmentation());
    if (m_segmentationPhiEta[m_systemId[iSys]] == nullptr && m_segmentationMulti[m_systemId[iSys]] == nullptr) {
      error() << "There is no phi-eta or multi- segmentation." << endmsg;
      return StatusCode::FAILURE;
    }
    m_decoder.insert(
        std::make_pair(m_systemId[iSys], m_geoSvc->getDetector()->readout(m_readoutName[iSys]).idSpec().decoder()));
  }
  // Initialize random service
  m_randSvc = service("RndmGenSvc", false);
  if (!m_randSvc) {
    error() << "Couldn't get RndmGenSvc!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  {
    StatusCode sc = m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.));
    if (sc.isFailure()) {
      error() << "Failed to initialize Gaussian random number generator!" << endmsg;
    }
  }

  // open and check file, read the histograms with noise constants
  if (initNoiseFromFile().isFailure()) {
    error() << "Couldn't open file with noise constants!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // calculate borders of eta bins:
  if (m_etaValues.size() != m_presamplerShiftP0.size() && m_etaValues.size() != m_presamplerShiftP1.size() &&
      m_etaValues.size() != m_presamplerScaleP0.size() && m_etaValues.size() != m_presamplerScaleP1.size()) {
    error() << "Sizes of parameter vectors for upstream energy correction should be the same" << endmsg;
    return StatusCode::FAILURE;
  }
  // if only one eta, leave border vector empty
  for (uint iEta = 1; iEta < m_etaValues.size(); iEta++) {
    m_etaBorders.push_back(m_etaValues[iEta - 1] + 0.5 * (m_etaValues[iEta] - m_etaValues[iEta - 1]));
  }
  // push values for the last eta bin, width as the last one
  if (m_etaValues.size() > 1) {
    m_etaBorders.push_back(m_etaValues[m_etaValues.size() - 1] +
                           0.5 * (m_etaValues[m_etaValues.size() - 1] - m_etaValues[m_etaValues.size() - 2]));
  } else {
    // high eta to ensure all particles fall below
    m_etaBorders.push_back(100);
  }
  // OPTIMISATION OF CLUSTER SIZE
  // sanity check
  if (!(m_nEtaFinal.size() == m_numLayers && m_nPhiFinal.size() == m_numLayers)) {
    error() << "Size of optimised window should be equal to number of layers: " << endmsg;
    error() << "Size of windows in eta:  " << m_nEtaFinal.size() << "\tsize of windows in phi:  " << m_nPhiFinal.size()
            << "number of layers:  " << m_numLayers << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_nEtaFinal.size() == m_numLayers) {
    for (uint iLayer = 0; iLayer < m_numLayers; iLayer++) {
      m_halfPhiFin.push_back(floor(m_nPhiFinal[iLayer] / 2));
      m_halfEtaFin.push_back(floor(m_nEtaFinal[iLayer] / 2));
    }
  }

  // INITIALIZE

  if (!m_towerTool.retrieve()) {
    error() << "Unable to retrieve the tower building tool." << endmsg;
    return StatusCode::FAILURE;
  }
  // Get number of calorimeter towers
  m_towerTool->towersNumber(m_nEtaTower, m_nPhiTower);
  debug() << "Number of calorimeter towers (eta x phi) : " << m_nEtaTower << " x " << m_nPhiTower << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode MassInv::execute(const EventContext&) const {
  // Get the input collection with clusters
  const edm4hep::ClusterCollection* inClusters = m_inClusters.get();
  edm4hep::ClusterCollection* correctedClusters = m_correctedClusters.createAndPut();

  // for single particle events compare with truth particles
  TVector3 momentum;
  double phiVertex = 0;
  double etaVertex = 0;
  double thetaVertex = 0;
  const auto particle = m_particle.get();
  if (particle->size() == 1) {
    for (const auto& part : *particle) {
      momentum = TVector3(part.getMomentum().x, part.getMomentum().y, part.getMomentum().z);
      etaVertex = momentum.Eta();
      phiVertex = momentum.Phi();
      thetaVertex = 2 * atan(exp(-etaVertex));
      verbose() << " vertex eta " << etaVertex << "   phi = " << phiVertex << " theta = " << thetaVertex << endmsg;
    }
  }

  // TODO change that so all systems can be used
  uint systemId = m_systemId[0];
  const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo* segmentation = nullptr;
  if (m_segmentationPhiEta[systemId] != nullptr) {
    segmentation = m_segmentationPhiEta[systemId];
  }

  std::vector<TLorentzVector> clustersMassInv;
  std::vector<TLorentzVector> clustersMassInvScaled;
  std::vector<TLorentzVector> clustersMassInvScaled2;
  std::vector<TLorentzVector> clustersMassInvScaled3;
  std::vector<TLorentzVector> clustersMassInvScaled4;
  std::vector<TLorentzVector> clustersMassInvScaled5;
  for (const auto& cluster : *inClusters) {
    double oldEnergy = 0;
    TVector3 pos(cluster.getPosition().x, cluster.getPosition().y, cluster.getPosition().z);
    double oldEta = pos.Eta();
    double oldPhi = pos.Phi();
    for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); cell++) {
      oldEnergy += cell->getEnergy();
    }
    verbose() << " OLD ENERGY = " << oldEnergy << " from " << cluster.hits_size() << " cells" << endmsg;
    verbose() << " OLD CLUSTER ENERGY = " << cluster.getEnergy() << endmsg;

    // Do everything only using the first defined calorimeter (default: Ecal barrel)
    double oldEtaId = -1;
    double oldPhiId = -1;
    if (m_segmentationPhiEta[systemId] != nullptr) {
      oldEtaId = int(floor((oldEta + 0.5 * segmentation->gridSizeEta() - segmentation->offsetEta()) /
                           segmentation->gridSizeEta()));
      oldPhiId = int(floor((oldPhi + 0.5 * segmentation->gridSizePhi() - segmentation->offsetPhi()) /
                           segmentation->gridSizePhi()));
    }

    // 0. Create new cluster, copy information from input
    auto newCluster = correctedClusters->create();
    double energy = 0;
    newCluster.setPosition(cluster.getPosition());
    for (auto cell = cluster.hits_begin(); cell != cluster.hits_end(); cell++) {
      if (m_segmentationMulti[systemId] != nullptr) {
        segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
            &m_segmentationMulti[systemId]->subsegmentation(cell->getCellID()));
        oldEtaId = int(floor((oldEta + 0.5 * segmentation->gridSizeEta() - segmentation->offsetEta()) /
                             segmentation->gridSizeEta()));
        oldPhiId = int(floor((oldPhi + 0.5 * segmentation->gridSizePhi() - segmentation->offsetPhi()) /
                             segmentation->gridSizePhi()));
      }
      if (m_decoder[systemId]->get(cell->getCellID(), "system") == systemId) {
        uint layerId = m_decoder[systemId]->get(cell->getCellID(), "layer");
        uint etaId = m_decoder[systemId]->get(cell->getCellID(), "eta");
        uint phiId = m_decoder[systemId]->get(cell->getCellID(), "phi");
        if (etaId >= (oldEtaId - m_halfEtaFin[layerId]) && etaId <= (oldEtaId + m_halfEtaFin[layerId]) &&
            phiId >= phiNeighbour((oldPhiId - m_halfPhiFin[layerId]), segmentation->phiBins()) &&
            phiId <= phiNeighbour((oldPhiId + m_halfPhiFin[layerId]), segmentation->phiBins())) {
          if (m_ellipseFinalCluster) {
            if (pow((etaId - oldEtaId) / (m_nEtaFinal[layerId] / 2.), 2) +
                    pow((phiId - oldPhiId) / (m_nPhiFinal[layerId] / 2.), 2) <
                1) {
              newCluster.addToHits(*cell);
              energy += cell->getEnergy();
            }
          } else {
            newCluster.addToHits(*cell);
            energy += cell->getEnergy();
          }
        }
      }
    }
    newCluster.setEnergy(energy);

    // 1. Correct eta position with log-weighting
    double newEta = 0;
    std::vector<double> sumEtaLayer;
    if (fabs(oldEta) > 1.5) {
      newEta = oldEta;
    } else {
      double sumEnFirstLayer = 0;
      // get current pseudorapidity
      std::vector<double> sumEnLayer;
      std::vector<double> sumWeightLayer;
      sumEnLayer.assign(m_numLayers, 0);
      sumEtaLayer.assign(m_numLayers, 0);
      sumWeightLayer.assign(m_numLayers, 0);
      // first check the energy deposited in each layer
      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++) {
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();
        uint layer = m_decoder[systemId]->get(cID, m_layerFieldName) + m_firstLayerId;
        sumEnLayer[layer] += cell->getEnergy();
      }
      sumEnFirstLayer = sumEnLayer[0];
      // repeat but calculating eta barycentre in each layer
      for (auto cell = newCluster.hits_begin(); cell != newCluster.hits_end(); cell++) {
        if (m_segmentationMulti[systemId] != nullptr) {
          segmentation = dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
              &m_segmentationMulti[systemId]->subsegmentation(cell->getCellID()));
        }
        dd4hep::DDSegmentation::CellID cID = cell->getCellID();
        uint layer = m_decoder[systemId]->get(cID, m_layerFieldName) + m_firstLayerId;
        double weightLog = std::max(0., m_etaRecalcLayerWeights[layer] + log(cell->getEnergy() / sumEnLayer[layer]));
        double eta = segmentation->eta(cell->getCellID());
        sumEtaLayer[layer] += (weightLog * eta);
        sumWeightLayer[layer] += weightLog;
        m_hDiffEtaHitLayer[layer]->Fill(eta - etaVertex);
      }
      // calculate eta position weighting with energy deposited in layer
      // this energy is a good estimator of 1/sigma^2 of (eta_barycentre-eta_MC) distribution
      for (uint iLayer = 0; iLayer < m_numLayers; iLayer++) {
        if (sumWeightLayer[iLayer] > 1e-10) {
          sumEtaLayer[iLayer] /= sumWeightLayer[iLayer];
          newEta += sumEtaLayer[iLayer] * sumEnLayer[iLayer];
          m_hDiffEtaLayer[iLayer]->Fill(sumEtaLayer[iLayer] - etaVertex);
        }
      }
      newEta /= energy;
      // alter Cartesian position of a cluster using new eta position
      double radius = pos.Perp();
      edm4hep::Vector3f newClusterPosition =
          edm4hep::Vector3f(radius * cos(oldPhi), radius * sin(oldPhi), radius * sinh(newEta));
      newCluster.setPosition(newClusterPosition);

      // 3. Correct for energy upstream
      // correct for presampler based on energy in the first layer layer:
      // check eta of the cluster and get correction parameters:
      double P00 = 0, P01 = 0, P10 = 0, P11 = 0;
      for (uint iEta = 0; iEta < m_etaBorders.size(); iEta++) {
        if (fabs(newEta) < m_etaBorders[iEta]) {
          P00 = m_presamplerShiftP0[iEta];
          P01 = m_presamplerShiftP1[iEta];
          P10 = m_presamplerScaleP0[iEta];
          P11 = m_presamplerScaleP1[iEta];
          break;
        }

        // if eta is larger than the last available eta values, take the last known parameters
        if (fabs(newEta) > m_etaBorders.back()) {
          error() << "cluster eta = " << newEta << " is larger than last defined eta values." << endmsg;
          return StatusCode::FAILURE;
        }
        double presamplerShift = P00 + P01 * cluster.getEnergy();
        double presamplerScale = P10 + P11 * sqrt(cluster.getEnergy());
        double energyFront = presamplerShift + presamplerScale * sumEnFirstLayer * m_samplingFraction[0];
        m_hUpstreamEnergy->Fill(energyFront);
        newCluster.setEnergy(newCluster.getEnergy() + energyFront);
      }
    }

    // 2. Correct energy for pileup noise
    uint numCells = newCluster.hits_size();
    double noise = 0;
    if (m_constPileupNoise == 0) {
      noise = getNoiseRMSPerCluster(newEta, numCells) * m_gauss.shoot() * std::sqrt(static_cast<int>(m_mu));
      verbose() << " NUM CELLS = " << numCells << "   cluster noise RMS = " << getNoiseRMSPerCluster(newEta, numCells)
                << " scaled to PU " << m_mu << "  = " << noise << endmsg;
    } else {
      noise = m_constPileupNoise * m_gauss.shoot() * std::sqrt(static_cast<int>(m_mu));
    }
    newCluster.setEnergy(newCluster.getEnergy() + noise);
    m_hPileupEnergy->Fill(noise);

    // Fill histograms
    m_hEnergyPreAnyCorrections->Fill(oldEnergy);
    m_hEnergyPostAllCorrections->Fill(newCluster.getEnergy());
    m_hEnergyPostAllCorrectionsAndScaling->Fill(newCluster.getEnergy() / m_response);

    // Position resolution
    m_hEta->Fill(newEta);
    m_hPhi->Fill(oldPhi);
    verbose() << " energy " << energy << "   numCells = " << numCells << " old energy = " << oldEnergy << " newEta "
              << newEta << "   phi = " << oldPhi << " theta = " << 2 * atan(exp(-newEta)) << endmsg;
    m_hNumCells->Fill(numCells);
    // // Calculate pointing resolution
    // TGraphErrors gZR = TGraphErrors();
    // for (uint iLayer = 0; iLayer < m_numLayers; iLayer++) {
    //   double theta = 2. * atan( exp(-sumEtaLayer[iLayer]) );
    //   gZR.SetPoint(iLayer, m_layerR[iLayer], m_layerR[iLayer] / tan(theta));
    //   gZR.SetPointError(iLayer, m_layerWidth[iLayer], 1/sumEnLayer[iLayer]);
    // }
    // auto result = gZR.Fit("pol1","S");
    // m_hDiffZ->Fill(result.Get()->Parameter(1) - zVertex);

    // Fill histograms for single particle events
    if (particle->size() == 1) {
      m_hDiffEta->Fill(newEta - etaVertex);
      m_hDiffPhi->Fill(oldPhi - phiVertex);
      m_hDiffTheta->Fill(2 * atan(exp(-newEta)) - thetaVertex);
    }

    // For invariant mass calculation
    if (m_energyAsThreshold) {
      if (energy / m_response > m_massInvThreshold) {
        TLorentzVector vec;
        vec.SetPtEtaPhiE(newCluster.getEnergy() / cosh(newEta), newEta, oldPhi, newCluster.getEnergy());
        clustersMassInv.push_back(vec);
        TLorentzVector vecScaled;
        vecScaled.SetPtEtaPhiE(newCluster.getEnergy() / m_response / cosh(newEta), newEta, oldPhi,
                               newCluster.getEnergy() / m_response);
        clustersMassInvScaled.push_back(vecScaled);
        clustersMassInvScaled2.push_back(vecScaled);
        clustersMassInvScaled3.push_back(vecScaled);
        clustersMassInvScaled4.push_back(vecScaled);
        clustersMassInvScaled5.push_back(vecScaled);
      }
    } else {
      if (energy / m_response / cosh(newEta) > m_massInvThreshold) {
        TLorentzVector vec;
        vec.SetPtEtaPhiE(newCluster.getEnergy() / cosh(newEta), newEta, oldPhi, newCluster.getEnergy());
        clustersMassInv.push_back(vec);
        TLorentzVector vecScaled;
        vecScaled.SetPtEtaPhiE(newCluster.getEnergy() / m_response / cosh(newEta), newEta, oldPhi,
                               newCluster.getEnergy() / m_response);
        clustersMassInvScaled.push_back(vecScaled);
        clustersMassInvScaled2.push_back(vecScaled);
        clustersMassInvScaled3.push_back(vecScaled);
        clustersMassInvScaled4.push_back(vecScaled);
        clustersMassInvScaled5.push_back(vecScaled);
      }
    }
    debug() << "pt of candidate: E " << energy << endmsg;
    debug() << "pt of candidate:resp  " << m_response << endmsg;
    debug() << "pt of candidate: newEta " << newEta << endmsg;
    debug() << "pt of candidate:en clu " << newCluster.getEnergy() << endmsg;
  }
  debug() << "Number of ALL candidates: " << inClusters->size() << endmsg;

  for (const auto& candidate1 : clustersMassInv) {
    for (const auto& candidate2 : clustersMassInv) {
      if (candidate1 != candidate2) {
        m_hMassInv->Fill((candidate1 + candidate2).Mag() * m_massInvCorrection);
        m_hDiPT->Fill((candidate1 + candidate2).Pt());
      }
    }
  }
  debug() << "Number of photon candidates: " << clustersMassInvScaled.size() << endmsg;
  // Sort (descending)
  if (clustersMassInvScaled.size() > 1) {
    std::sort(clustersMassInvScaled.begin(), clustersMassInvScaled.end(),
              [](TLorentzVector photon1, TLorentzVector photon2) { return photon1.Pt() > photon2.Pt(); });
    double diPhotonMass = (clustersMassInvScaled[0] + clustersMassInvScaled[1]).Mag() * m_massInvCorrection;
    double diPhotonPt = (clustersMassInvScaled[0] + clustersMassInvScaled[1]).Pt();
    m_hDiPTScaled->Fill(diPhotonPt);
    m_hMassInvScaled->Fill(diPhotonMass);
    m_hMassInvScaledPt->Fill(diPhotonMass, diPhotonPt);
    if (diPhotonPt > 100) {
      m_hMassInvScaled100->Fill(diPhotonMass);
    }
    if (diPhotonPt > 200) {
      m_hMassInvScaled200->Fill(diPhotonMass);
    }
    if (diPhotonPt > 300) {
      m_hMassInvScaled300->Fill(diPhotonMass);
    }

    // create towers
    m_towers.assign(m_nEtaTower, std::vector<float>(m_nPhiTower, 0));
    m_towerTool->buildTowers(m_towers);
    // check all isolation windows around photons
    debug() << "Number of photon candidates: " << clustersMassInvScaled.size() << endmsg;
    for (uint iCluster = 0; iCluster < m_etaSizes.size(); iCluster++) {
      debug() << "Size of the reconstruction window (eta,phi) " << m_etaSizes[iCluster] << ", " << m_phiSizes[iCluster]
              << endmsg;
      int halfEtaWin = floor(m_etaSizes[iCluster] / 2.);
      int halfPhiWin = floor(m_phiSizes[iCluster] / 2.);
      debug() << "Half-size of the reconstruction window (eta,phi) " << halfEtaWin << ", " << halfPhiWin << endmsg;
      for (auto photonCandidate = clustersMassInvScaled.begin(); photonCandidate != clustersMassInvScaled.end();
           photonCandidate++) {
        uint photonIdEta = m_towerTool->idEta(photonCandidate->Eta());
        uint photonIdPhi = m_towerTool->idPhi(photonCandidate->Phi());
        // LOOK AROUND
        double sumWindow = 0;
        for (size_t iEtaWindow = photonIdEta - halfEtaWin; iEtaWindow <= photonIdEta + halfEtaWin; iEtaWindow++) {
          for (size_t iPhiWindow = photonIdPhi - halfPhiWin; iPhiWindow <= photonIdPhi + halfPhiWin; iPhiWindow++) {
            sumWindow += m_towers[iEtaWindow][phiNeighbour(iPhiWindow, m_nPhiTower)];
          }
        }
        m_hHCalEnergy->Fill(sumWindow);
        if (sumWindow > m_hcalEnergyThreshold) {
          clustersMassInvScaled.erase(photonCandidate);
          photonCandidate--;
        }
      }
      for (auto photonCandidate = clustersMassInvScaled2.begin(); photonCandidate != clustersMassInvScaled2.end();
           photonCandidate++) {
        uint photonIdEta = m_towerTool->idEta(photonCandidate->Eta());
        uint photonIdPhi = m_towerTool->idPhi(photonCandidate->Phi());
        // LOOK AROUND
        double sumWindow = 0;
        for (size_t iEtaWindow = photonIdEta - halfEtaWin; iEtaWindow <= photonIdEta + halfEtaWin; iEtaWindow++) {
          for (size_t iPhiWindow = photonIdPhi - halfPhiWin; iPhiWindow <= photonIdPhi + halfPhiWin; iPhiWindow++) {
            sumWindow += m_towers[iEtaWindow][phiNeighbour(iPhiWindow, m_nPhiTower)];
          }
        }
        m_hHCalEnergy->Fill(sumWindow);
        if (sumWindow > m_hcalEnergyThreshold * 0.1) {
          clustersMassInvScaled2.erase(photonCandidate);
          photonCandidate--;
        }
      }
      for (auto photonCandidate = clustersMassInvScaled3.begin(); photonCandidate != clustersMassInvScaled3.end();
           photonCandidate++) {
        uint photonIdEta = m_towerTool->idEta(photonCandidate->Eta());
        uint photonIdPhi = m_towerTool->idPhi(photonCandidate->Phi());
        // LOOK AROUND
        double sumWindow = 0;
        for (size_t iEtaWindow = photonIdEta - halfEtaWin; iEtaWindow <= photonIdEta + halfEtaWin; iEtaWindow++) {
          for (size_t iPhiWindow = photonIdPhi - halfPhiWin; iPhiWindow <= photonIdPhi + halfPhiWin; iPhiWindow++) {
            sumWindow += m_towers[iEtaWindow][phiNeighbour(iPhiWindow, m_nPhiTower)];
          }
        }
        m_hHCalEnergy->Fill(sumWindow);
        if (sumWindow > m_hcalEnergyThreshold * 0.2) {
          clustersMassInvScaled3.erase(photonCandidate);
          photonCandidate--;
        }
      }
      for (auto photonCandidate = clustersMassInvScaled4.begin(); photonCandidate != clustersMassInvScaled4.end();
           photonCandidate++) {
        uint photonIdEta = m_towerTool->idEta(photonCandidate->Eta());
        uint photonIdPhi = m_towerTool->idPhi(photonCandidate->Phi());
        // LOOK AROUND
        double sumWindow = 0;
        for (size_t iEtaWindow = photonIdEta - halfEtaWin; iEtaWindow <= photonIdEta + halfEtaWin; iEtaWindow++) {
          for (size_t iPhiWindow = photonIdPhi - halfPhiWin; iPhiWindow <= photonIdPhi + halfPhiWin; iPhiWindow++) {
            sumWindow += m_towers[iEtaWindow][phiNeighbour(iPhiWindow, m_nPhiTower)];
          }
        }
        m_hHCalEnergy->Fill(sumWindow);
        if (sumWindow > m_hcalEnergyThreshold * 0.3) {
          clustersMassInvScaled4.erase(photonCandidate);
          photonCandidate--;
        }
      }
      for (auto photonCandidate = clustersMassInvScaled5.begin(); photonCandidate != clustersMassInvScaled5.end();
           photonCandidate++) {
        uint photonIdEta = m_towerTool->idEta(photonCandidate->Eta());
        uint photonIdPhi = m_towerTool->idPhi(photonCandidate->Phi());
        // LOOK AROUND
        double sumWindow = 0;
        for (size_t iEtaWindow = photonIdEta - halfEtaWin; iEtaWindow <= photonIdEta + halfEtaWin; iEtaWindow++) {
          for (size_t iPhiWindow = photonIdPhi - halfPhiWin; iPhiWindow <= photonIdPhi + halfPhiWin; iPhiWindow++) {
            sumWindow += m_towers[iEtaWindow][phiNeighbour(iPhiWindow, m_nPhiTower)];
          }
        }
        m_hHCalEnergy->Fill(sumWindow);
        if (sumWindow > m_hcalEnergyThreshold * 0.4) {
          clustersMassInvScaled5.erase(photonCandidate);
          photonCandidate--;
        }
      }
    }
    double diPhotonMassIsolated = (clustersMassInvScaled[0] + clustersMassInvScaled[1]).Mag() * m_massInvCorrection;
    double diPhotonPtIsolated = (clustersMassInvScaled[0] + clustersMassInvScaled[1]).Pt();
    m_hMassInvScaledIsolated->Fill(diPhotonMassIsolated);
    if (diPhotonPtIsolated > 100) {
      m_hMassInvScaledIsolated100->Fill(diPhotonMassIsolated);
    }
    if (diPhotonPtIsolated > 200) {
      m_hMassInvScaledIsolated200->Fill(diPhotonMassIsolated);
    }
    if (diPhotonPtIsolated > 300) {
      m_hMassInvScaledIsolated300->Fill(diPhotonMassIsolated);
    }
    double diPhotonMassIsolated2 = (clustersMassInvScaled2[0] + clustersMassInvScaled2[1]).Mag() * m_massInvCorrection;
    double diPhotonPtIsolated2 = (clustersMassInvScaled2[0] + clustersMassInvScaled2[1]).Pt();
    m_hMassInvScaledIsolated2->Fill(diPhotonMassIsolated2);
    if (diPhotonPtIsolated2 > 100) {
      m_hMassInvScaledIsolated2100->Fill(diPhotonMassIsolated2);
    }
    if (diPhotonPtIsolated2 > 200) {
      m_hMassInvScaledIsolated2200->Fill(diPhotonMassIsolated2);
    }
    if (diPhotonPtIsolated2 > 300) {
      m_hMassInvScaledIsolated2300->Fill(diPhotonMassIsolated2);
    }
    double diPhotonMassIsolated3 = (clustersMassInvScaled3[0] + clustersMassInvScaled3[1]).Mag() * m_massInvCorrection;
    double diPhotonPtIsolated3 = (clustersMassInvScaled3[0] + clustersMassInvScaled3[1]).Pt();
    m_hMassInvScaledIsolated3->Fill(diPhotonMassIsolated3);
    if (diPhotonPtIsolated3 > 100) {
      m_hMassInvScaledIsolated3100->Fill(diPhotonMassIsolated3);
    }
    if (diPhotonPtIsolated3 > 200) {
      m_hMassInvScaledIsolated3200->Fill(diPhotonMassIsolated3);
    }
    if (diPhotonPtIsolated3 > 300) {
      m_hMassInvScaledIsolated3300->Fill(diPhotonMassIsolated3);
    }
    double diPhotonMassIsolated4 = (clustersMassInvScaled4[0] + clustersMassInvScaled4[1]).Mag() * m_massInvCorrection;
    double diPhotonPtIsolated4 = (clustersMassInvScaled4[0] + clustersMassInvScaled4[1]).Pt();
    m_hMassInvScaledIsolated4->Fill(diPhotonMassIsolated4);
    if (diPhotonPtIsolated4 > 100) {
      m_hMassInvScaledIsolated4100->Fill(diPhotonMassIsolated4);
    }
    if (diPhotonPtIsolated4 > 200) {
      m_hMassInvScaledIsolated4200->Fill(diPhotonMassIsolated4);
    }
    if (diPhotonPtIsolated4 > 300) {
      m_hMassInvScaledIsolated4300->Fill(diPhotonMassIsolated4);
    }
    double diPhotonMassIsolated5 = (clustersMassInvScaled5[0] + clustersMassInvScaled5[1]).Mag() * m_massInvCorrection;
    double diPhotonPtIsolated5 = (clustersMassInvScaled5[0] + clustersMassInvScaled5[1]).Pt();
    m_hMassInvScaledIsolated5->Fill(diPhotonMassIsolated5);
    if (diPhotonPtIsolated5 > 100) {
      m_hMassInvScaledIsolated5100->Fill(diPhotonMassIsolated5);
    }
    if (diPhotonPtIsolated5 > 200) {
      m_hMassInvScaledIsolated5200->Fill(diPhotonMassIsolated5);
    }
    if (diPhotonPtIsolated5 > 300) {
      m_hMassInvScaledIsolated5300->Fill(diPhotonMassIsolated5);
    }
    debug() << "Number of photon candidates: " << clustersMassInvScaled.size() << endmsg;
    debug() << "Number of photon candidates: " << clustersMassInvScaled2.size() << endmsg;
    debug() << "Number of photon candidates: " << clustersMassInvScaled3.size() << endmsg;
    debug() << "Number of photon candidates: " << clustersMassInvScaled4.size() << endmsg;
    debug() << "Number of photon candidates: " << clustersMassInvScaled5.size() << endmsg;
  }
  return StatusCode::SUCCESS;
}

StatusCode MassInv::finalize() { return Gaudi::Algorithm::finalize(); }

StatusCode MassInv::initNoiseFromFile() {
  // Check if file exists
  if (m_noiseFileName.empty()) {
    error() << "Name of the file with the noise values not provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (gSystem->AccessPathName(m_noiseFileName.value().c_str())) {
    error() << "Provided file with the noise values not found!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  }
  std::unique_ptr<TFile> inFile(TFile::Open(m_noiseFileName.value().c_str(), "READ"));
  if (inFile->IsZombie()) {
    error() << "Unable to open the file with the noise values!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Using the following file with the noise constants: " << m_noiseFileName.value() << endmsg;
  }

  std::string pileupParamHistoName;
  // Read the histograms with parameters for the pileup noise from the file
  for (unsigned i = 0; i < 2; i++) {
    pileupParamHistoName = m_pileupHistoName + std::to_string(i);
    debug() << "Getting histogram with a name " << pileupParamHistoName << endmsg;
    m_histoPileupConst.push_back(*dynamic_cast<TH1F*>(inFile->Get(pileupParamHistoName.c_str())));
    if (m_histoPileupConst.at(i).GetNbinsX() < 1) {
      error() << "Histogram  " << pileupParamHistoName
              << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // Check if we have same number of histograms (all layers) for pileup and electronics noise
  if (m_histoPileupConst.size() == 0) {
    error() << "No histograms with noise found!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

double MassInv::getNoiseRMSPerCluster(double aEta, uint aNumCells) const {
  double param0 = 0.;
  double param1 = 0.;

  // All histograms have same binning, all bins with same size
  // Using the histogram of the first parameter to get the bin size
  unsigned index = 0;
  if (m_histoPileupConst.size() != 0) {
    int Nbins = m_histoPileupConst.at(index).GetNbinsX();
    double deltaEtaBin =
        (m_histoPileupConst.at(index).GetBinLowEdge(Nbins) + m_histoPileupConst.at(index).GetBinWidth(Nbins) -
         m_histoPileupConst.at(index).GetBinLowEdge(1)) /
        Nbins;
    double etaFirtsBin = m_histoPileupConst.at(index).GetBinLowEdge(1);
    // find the eta bin for the cell
    int ibin = floor((fabs(aEta) - etaFirtsBin) / deltaEtaBin) + 1;
    verbose() << "Current eta = " << aEta << " bin = " << ibin << endmsg;
    if (ibin > Nbins) {
      debug() << "eta outside range of the histograms! Cell eta: " << aEta << " Nbins in histogram: " << Nbins
              << endmsg;
      ibin = Nbins;
    }
    param0 = m_histoPileupConst.at(0).GetBinContent(ibin);
    param1 = m_histoPileupConst.at(1).GetBinContent(ibin);
    verbose() << "p0 = " << param0 << " param1 = " << param1 << endmsg;
  } else {
    debug() << "No histograms with noise constants!!!!! " << endmsg;
  }
  double pileupNoise = param0 * pow(aNumCells * (m_dEta / 0.01), param1);

  return pileupNoise;
}

unsigned int MassInv::phiNeighbour(int aIPhi, int aMaxPhi) const {
  if (aIPhi < 0) {
    return aMaxPhi + aIPhi;
  } else if (aIPhi >= aMaxPhi) {
    return aIPhi % aMaxPhi;
  }
  return aIPhi;
}
