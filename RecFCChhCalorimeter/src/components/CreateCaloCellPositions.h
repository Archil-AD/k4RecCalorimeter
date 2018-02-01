#ifndef DETCOMPONENTS_CREATECELLPOSITIONS_H
#define DETCOMPONENTS_CREATECELLPOSITIONS_H

// FCCSW
#include "FWCore/DataHandle.h"
#include "RecInterface/ICellPositionsTool.h"

// Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"

#include "datamodel/CaloHit.h"
#include "datamodel/CaloHitCollection.h"
#include "datamodel/PositionedCaloHit.h"
#include "datamodel/PositionedCaloHitCollection.h"

class IGeoSvc;

/** @class CreateCaloCellPositions Reconstruction/RecCalorimeter/src/components/CreateCaloCellPositions.h CreateCaloCellPositions.h
 *
 *  Retrieve positions of the cells from cell ID.
 *  This algorithm saves the centre position of the volume. Defined for all Calo-Subsystems within tools.
 *
 *  @author Coralie Neubueser
 *
 */

class CreateCaloCellPositions : public GaudiAlgorithm {

public:
  CreateCaloCellPositions(const std::string& name, ISvcLocator* svcLoc);
  /**  Initialize.
   *   @return status code
   */
  StatusCode initialize();
  /**  Execute.
   *   @return status code
   */
  StatusCode execute();
  /**  Finalize.
   *   @return status code
   */
  StatusCode finalize();

private:
  /// Handle for tool to get positions in ECal Barrel
  ToolHandle<ICellPositionsTool> m_cellPositionsECalBarrelTool{"CellPositionsECalBarrelTool", this};
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalBarrelTool{"CellPositionsHCalBarrelNoSegTool", this};
  /// Handle for tool to get positions in HCal Barrel and Ext Barrel, no Segmentation
  ToolHandle<ICellPositionsTool> m_cellPositionsHCalExtBarrelTool{"CellPositionsHCalBarrelNoSegTool", this};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsEMECTool{"CellPositionsCaloDiscsTool", this};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsHECTool{"CellPositionsCaloDiscsTool", this};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsEMFwdTool{"CellPositionsCaloDiscsTool", this};
  /// Handle for tool to get positions in Calo Discs
  ToolHandle<ICellPositionsTool> m_cellPositionsHFwdTool{"CellPositionsCaloDiscsTool", this};
  /// Decoder for system ID
  DD4hep::DDSegmentation::BitField64* m_decoder = new DD4hep::DDSegmentation::BitField64("system:4");
  /// Input collection
  DataHandle<fcc::CaloHitCollection> m_hits{"hits/hits", Gaudi::DataHandle::Reader, this};
  /// Output collection
  DataHandle<fcc::PositionedCaloHitCollection> m_positionedHits{"hits/positionedHits", Gaudi::DataHandle::Writer, this};
};

#endif /* DETCOMPONENTS_CREATECELLPOSITIONS_H */
