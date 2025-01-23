#ifndef RECCALORIMETER_MUONCALOHITDIGI_H
#define RECCALORIMETER_MUONCALOHITDIGI_H

// GAUDI
#include "Gaudi/Property.h"
#include "Gaudi/Algorithm.h"

// K4FWCORE
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4Interface/IGeoSvc.h"

// edm4hep
#include "DD4hep/Detector.h"  // for dd4hep::VolumeManager
#include "DD4hep/Readout.h"
#include "DD4hep/Volumes.h"
#include "DDSegmentation/Segmentation.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include <DDRec/DetectorData.h>

/** @class MuonCaloHitDigi
 *
 *  Algorithm for creating digitized/reconstructed Muon tagger hits from Geant4 hits (edm4hep::SimTrackerHit).
 *
 *  @author Archil Durglishvili
 *  @date   2025-01
 *
 */

class MuonCaloHitDigi : public Gaudi::Algorithm {
public:
  explicit MuonCaloHitDigi(const std::string&, ISvcLocator*);
  virtual ~MuonCaloHitDigi();
  /**  Initialize.
   *   @return status code
   */
  virtual StatusCode initialize() final;
  /**  Execute.
   *   @return status code
   */
  virtual StatusCode execute(const EventContext&) const final;
  /**  Finalize.
   *   @return status code
   */
  virtual StatusCode finalize() final;

private:
  // Input SimTrackerHit collection name
  mutable DataHandle<edm4hep::SimTrackerHitCollection> m_input_sim_hits{"inputSimHits", Gaudi::DataHandle::Reader, this};
  // Output Calo hit collection name
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_output_digi_hits{"outputDigiHits", Gaudi::DataHandle::Writer, this};
  MetaDataHandle<std::string> m_cellsCellIDEncoding{m_output_digi_hits, edm4hep::labels::CellIDEncoding, Gaudi::DataHandle::Writer};

  // Name of the Muon Barrel det.
  Gaudi::Property<std::string> m_barrelDetectorName{this, "barrelDetectorName", "MuonTaggerBarrel"};
  // Name of the Muon Endcap det.
  Gaudi::Property<std::string> m_endcapDetectorName{this, "endcapDetectorName", "MuonTaggerEndcap"};
  // Detector readout names
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "MuonTaggerPhiTheta", "Name of the detector readout"};
  // Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  // Decoder for the cellID
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  // Theta-phi segmentation
  dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo* m_segmentation = nullptr;
  // Barrel systemId
  Gaudi::Property<uint> m_barrelSysId{this, "barrelSysId", 12};
  // vector to store layer positions for Barrel and Endcap
  std::vector<std::vector<double> > m_layerPositions;
  // layers retrieved from the geometry 
  const std::vector<dd4hep::rec::LayeredCalorimeterStruct::Layer>* m_layersRetrieved = nullptr;
  // Maps of cell IDs (corresponding to DD4hep IDs) on final energies
  mutable std::unordered_map<uint64_t, double> m_cellsMap;
};

#endif /* RECCALORIMETER_MUONCALOHITDIGI_H */
