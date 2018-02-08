#include "RewriteHCalBarrelBitfield.h"

// FCCSW
#include "DetInterface/IGeoSvc.h"

// datamodel
#include "datamodel/CaloHitCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DDSegmentation/Segmentation.h"

#include <math.h>

DECLARE_ALGORITHM_FACTORY(RewriteHCalBarrelBitfield)

RewriteHCalBarrelBitfield::RewriteHCalBarrelBitfield(const std::string& aName, ISvcLocator* aSvcLoc) : GaudiAlgorithm(aName, aSvcLoc) {
  declareProperty("inhits", m_inHits, "Hit collection with old segmentation (input)");
  declareProperty("outhits", m_outHits, "Hit collection with modified segmentation (output)");
}

RewriteHCalBarrelBitfield::~RewriteHCalBarrelBitfield() {}

StatusCode RewriteHCalBarrelBitfield::initialize() {
  if (GaudiAlgorithm::initialize().isFailure()) return StatusCode::FAILURE;
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  // check if readouts exist
  if (m_geoSvc->lcdd()->readouts().find(m_oldReadoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_oldReadoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_geoSvc->lcdd()->readouts().find(m_newReadoutName) == m_geoSvc->lcdd()->readouts().end()) {
    error() << "Readout <<" << m_newReadoutName << ">> does not exist." << endmsg;
    return StatusCode::FAILURE;
  }
  // Take readout, bitfield from GeoSvc
  m_oldDecoder = m_geoSvc->lcdd()->readout(m_oldReadoutName).idSpec().decoder();
  // segmentation identifiers to be overwritten
  if (m_oldIdentifiers.size() == 0) {
    // it is not an error, maybe no segmentation was used previously
    info() << "No identifiers to remove. Only rewritting the readout" << endmsg;
  }
  // create detector identifiers (= all bitfield ids - segmentation ids)
  for (uint itField = 0; itField < m_oldDecoder->size(); itField++) {
    std::string field = (*m_oldDecoder)[itField].name();
    auto iter = std::find(m_oldIdentifiers.begin(), m_oldIdentifiers.end(), field);
    if (iter == m_oldIdentifiers.end()) {
      m_detectorIdentifiers.push_back(field);
    }
  }
  std::vector<std::string> newFields;
  m_newDecoder =  m_geoSvc->lcdd()->readout(m_newReadoutName).idSpec().decoder();
  for (uint itField = 0; itField < m_newDecoder->size(); itField++) {
    newFields.push_back((*m_newDecoder)[itField].name());
  }
  for (const auto& detectorField : m_detectorIdentifiers) {
    auto iter = std::find(newFields.begin(), newFields.end(), detectorField);
    if (iter == newFields.end()) {
      error() << "New readout does not contain field <<" << detectorField << ">> that describes the detector ID."
	      << endmsg;
      return StatusCode::FAILURE;
    }
  }
  info() << "Rewritting the readout bitfield." << endmsg;
  info() << "Old bitfield:\t" << m_oldDecoder->fieldDescription() << endmsg;
  info() << "New bitfield:\t" << m_newDecoder->fieldDescription() << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode RewriteHCalBarrelBitfield::execute() {
  const auto inHits = m_inHits.get();
  auto outHits = m_outHits.createAndPut();
  // loop over positioned hits to get the energy deposits: position and cellID
  // cellID contains the volumeID that needs to be copied to the new id
  uint64_t oldid = 0;
  uint debugIter = 0;
  uint offsetPhiIds = 3*256/4.+1; 
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridPhiEta*>(
      m_geoSvc->lcdd()->readout(m_oldReadoutName).segmentation().segmentation());

  for (const auto& hit : *inHits) {
    fcc::CaloHit newHit = outHits->create();
    newHit.energy(hit.energy());
    newHit.time(hit.time());
    m_oldDecoder->setValue(hit.cellId());
    auto inSegPhi = m_segmentation->phi(hit.cellId());
    if (debugIter < m_debugPrint) {
      debug() << "OLD: " << m_oldDecoder->valueString() << endmsg;
    }
    // now rewrite all fields except for those to be removed
    for (const auto& detectorField : m_detectorIdentifiers) {
      oldid = (*m_oldDecoder)[detectorField];
      (*m_newDecoder)[detectorField] = oldid;
      if (detectorField == "module"){
	// get phi id
	auto phiId = (*m_oldDecoder)["phi"];
	uint newModuleId = abs(phiId - offsetPhiIds);
	(*m_newDecoder)[detectorField] = newModuleId;
	if( (inSegPhi<0.025 && inSegPhi>0.025) || (inSegPhi > (M_PI/2.-0.025)&&inSegPhi<(M_PI/2.+0.025)) || (inSegPhi > (M_PI-0.025)&&inSegPhi<(M_PI+0.025)) || (inSegPhi > (-M_PI/2.-0.025)&&inSegPhi<(-M_PI/2.+0.025))){
	  info() << "Phi of Cell   : " << inSegPhi << endmsg;
	  info() << "old module id : " << oldid << endmsg;
	  info() << "old phi    id : " << phiId << endmsg;
	  info() << "new module id : " << newModuleId << endmsg;
	}
      }else{ 
	(*m_newDecoder)[detectorField] = oldid;
      } 
    }
    newHit.cellId(m_newDecoder->getValue());
    if (debugIter < m_debugPrint) {
      debug() << "NEW: " << m_newDecoder->valueString() << endmsg;
      debugIter++;
    }
  }
  return StatusCode::SUCCESS;
}

StatusCode RewriteHCalBarrelBitfield::finalize() { return GaudiAlgorithm::finalize(); }
