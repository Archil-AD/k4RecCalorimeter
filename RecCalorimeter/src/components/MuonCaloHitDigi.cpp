#include "MuonCaloHitDigi.h"
#include "DDSegmentation/Segmentation.h"
#include "detectorSegmentations/FCCSWGridPhiTheta_k4geo.h"
#include <DDRec/DetectorData.h>
#include "DD4hep/Detector.h"

using dd4hep::DetElement;

DECLARE_COMPONENT(MuonCaloHitDigi)

MuonCaloHitDigi::MuonCaloHitDigi(const std::string& aName, ISvcLocator* aSvcLoc)
    : Gaudi::Algorithm(aName, aSvcLoc), m_geoSvc("GeoSvc", "MuonCaloHitDigi") {
  declareProperty("inputSimHits", m_input_sim_hits, "Input sim tracker hit collection name");
  declareProperty("outputDigiHits", m_output_digi_hits, "Output calo hit collection name");
}

MuonCaloHitDigi::~MuonCaloHitDigi() {}

StatusCode MuonCaloHitDigi::initialize() {
  // check if readout exists
  if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
    error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }

  // get PhiTheta segmentation
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiTheta_k4geo *>(
        m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());

  if (m_segmentation == nullptr)
  {
    error() << "There is no phi-theta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  // set the cellID decoder
  m_decoder = m_geoSvc->getDetector()->readout(m_readoutName).idSpec().decoder(); // Can be used to access e.g. layer index: m_decoder->get(cellID, "layer"),
  m_cellsCellIDEncoding.put( m_decoder->fieldDescription() );

  //--------------------------------------
  // retrieve layer positions from the LayeredCalorimeterData extension
  dd4hep::Detector* detector = m_geoSvc->getDetector();
  if (!detector)
  {
    error() << "Unable to retrieve the detector." << endmsg;
    return StatusCode::FAILURE;
  }

  DetElement caloDetElemBarrel = detector->detector(m_barrelDetectorName);
  if (!caloDetElemBarrel.isValid())
  {
    error() << "Unable to retrieve the detector element: " << m_barrelDetectorName << endmsg;
    return StatusCode::FAILURE;
  }

  DetElement caloDetElemEndcap = detector->detector(m_endcapDetectorName);
  if (!caloDetElemEndcap.isValid())
  {
    error() << "Unable to retrieve the detector element: " << m_endcapDetectorName << endmsg;
    return StatusCode::FAILURE;
  }

  std::vector<dd4hep::rec::LayeredCalorimeterData* > theExtension;
  theExtension.push_back( caloDetElemBarrel.extension<dd4hep::rec::LayeredCalorimeterData>() );
  theExtension.push_back( caloDetElemEndcap.extension<dd4hep::rec::LayeredCalorimeterData>() );
  if (!theExtension[0] || !theExtension[1])
  {
    error() << "The detector element does not have the required LayeredCalorimeterData extension." << endmsg;
    return StatusCode::FAILURE;
  }

  // prepare layer positions vectors for barrel end endcap
  m_layerPositions.resize(2);

  // fill layer positions
  for(int iDet = 0; iDet < 2; iDet++)
  {
    m_layersRetrieved = &(theExtension[iDet]->layers) ;
    for (unsigned int idxLayer = 0; idxLayer < m_layersRetrieved->size(); ++idxLayer) {
      const dd4hep::rec::LayeredCalorimeterStruct::Layer & theLayer = m_layersRetrieved->at(idxLayer);
      // layer inner position (for barrel it's inner radius while for endcap - smallest distance in z)
      double layerInnerPosition = theLayer.distance;
      // layer thickness
      double layerThickness = theLayer.sensitive_thickness;

      // Calculate the position of the center of current layer
      double position = layerInnerPosition + layerThickness / 2.;

      m_layerPositions[iDet].push_back(position);
    }
  }
  //--------------------------------------

  return StatusCode::SUCCESS;
}

StatusCode MuonCaloHitDigi::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
  const edm4hep::SimTrackerHitCollection* input_sim_hits = m_input_sim_hits.get();
  verbose() << "Input Sim Hit collection size: " << input_sim_hits->size() << endmsg;

  // Digitize the sim hits
  edm4hep::CalorimeterHitCollection* output_digi_hits = m_output_digi_hits.createAndPut();

  for (const auto& input_sim_hit : *input_sim_hits) {
    // retrieve the cellID
    dd4hep::DDSegmentation::CellID cellID = input_sim_hit.getCellID();
    debug() << "Digitisation of " << m_readoutName << ", cellID: " << cellID << endmsg;
    m_cellsMap[cellID] += input_sim_hit.getEDep();
  }

  // create cells
  for (const auto& cell : m_cellsMap) {
    dd4hep::DDSegmentation::CellID cellID = cell.first;
    // get systemId
    int systemId = m_decoder->get(cellID, "system");
    // get layer
    int layer = m_decoder->get(cellID, "layer");

    double radius = 0;

    // get radius for barrel
    if(systemId==m_barrelSysId)
    {
      radius = m_layerPositions[0][layer];
    }
    // get radius for endcap
    else
    {
      auto theta = m_segmentation->theta(cellID);
      radius = std::sin(theta)/std::cos(theta) * m_layerPositions[1][layer];
    }

    // get local position (for r=1)
    auto localPos = m_segmentation->position(cellID);

    // create the cell position
    // NOTE! AD: z-position of some cells at the adges might be out of range. This is a bad feature of FCCSWGridPhiTheta_k4geo segmentation.
    edm4hep::Vector3f hitPosition(localPos.x() * radius, localPos.y() * radius, localPos.z() * radius);

    // create output digitized hit
    auto output_digi_hit = output_digi_hits->create();

    output_digi_hit.setCellID(cellID);
    output_digi_hit.setEnergy(cell.second);
    output_digi_hit.setPosition(hitPosition);

    debug() << "Position of digi hit (mm) : \t" << output_digi_hit.getPosition().x/dd4hep::mm
                                      << "\t" << output_digi_hit.getPosition().y/dd4hep::mm
                                      << "\t" << output_digi_hit.getPosition().z/dd4hep::mm << endmsg;
  }
  debug() << "Output Cell collection size: " << output_digi_hits->size() << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode MuonCaloHitDigi::finalize() { return StatusCode::SUCCESS; }
