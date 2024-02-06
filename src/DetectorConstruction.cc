//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include <cuda_runtime.h>
#include "SEventConfig.hh"
#include "G4CXOpticks.hh"
#include "U4Physics.hh"
#include <cuda_runtime.h>
#include <globals.hh>
#include "U4Scint.h"
#include "U4Material.hh"
#include "MaterialsList.h"
#include "OpticalMaterialProperties.h"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "U4SensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include "PhotonSD.hh"
#include "G4SDManager.hh"
namespace B1 {

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    using namespace CLHEP;

    G4VPhysicalVolume *DetectorConstruction::Construct() {




        // Option to switch on/off checking of volumes overlaps
        //
        G4bool checkOverlaps = true;

        //
        // World
        //
        G4Material *GXenon = materials::GXe(1 * bar, 293 * kelvin);
        GXenon->SetMaterialPropertiesTable(opticalprops::GXe(1 * bar, 293 * kelvin, 30000));

        G4Material *MgF2 = materials::MgF2();
        MgF2->SetMaterialPropertiesTable(opticalprops::MgF2());


        G4Material *Steel = materials::Steel();
        Steel->SetMaterialPropertiesTable(opticalprops::STEEL());

        G4double env_sizeXY = 20 * cm, env_sizeZ = 30 * cm;
        G4double world_sizeXY = 1.2 * env_sizeXY;
        G4double world_sizeZ = 1.2 * env_sizeZ;
        auto solidWorld = new G4Box("World",                           // its name
                                    0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size


        auto CubeDetector_Solid = new G4Box("DetectorSolid", 1 * cm, 1 * cm, 1/2 * mm);  // its size
        auto CubeDetectorLogic = new G4LogicalVolume(CubeDetector_Solid, MgF2, "Detector_Logic");
        auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
                                              GXenon,                                       // its material
                                              "World");                                        // its name



        auto physWorld = new G4PVPlacement(nullptr,  // no rotation
                                           G4ThreeVector(),                           // at (0,0,0)
                                           logicWorld,                                // its logical volume
                                           "World",                                   // its name
                                           nullptr,                                   // its mother  volume
                                           false,                                     // no boolean operation
                                           0,                                         // copy number
                                           checkOverlaps);                            // overlaps checking

        //
        // Envelope
        //
        auto solidEnv = new G4Box("GasXe_Solid",                    // its name
                                  0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

        auto logicEnv = new G4LogicalVolume(solidEnv, GXenon, "GasXe_Logic");


        auto Steel_Cover_solid = new G4Box("SteelCover_Solid",                    // its name
                                           0.5 * env_sizeXY + 5 * CLHEP::mm, 0.5 * env_sizeXY + 10 * CLHEP::mm,
                                           0.5 * env_sizeZ + 5 * CLHEP::mm);  // its size
        auto Steel_Cover_logic = new G4LogicalVolume(Steel_Cover_solid,  // its solid
                                                     Steel,                                     // its material
                                                     "SteelCover_logic");                                 // its name
        //
        //always return the physical World
        //
        auto SteelPlace = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "SteelCover", Steel_Cover_logic, physWorld, 0, 0,
                                            0);
        auto GasXePlace = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "GasXe", logicEnv, SteelPlace, 0, 0, 0);
        auto DetectorPlace = new G4PVPlacement(0, G4ThreeVector(0, 0, -10 * cm), "Detector", CubeDetectorLogic, GasXePlace, 0, 0,0);

        logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
        Steel_Cover_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
        G4VisAttributes gasVis = new G4VisAttributes();
        gasVis.SetColor(0, 1, 1);
        gasVis.SetForceCloud(true);
        logicEnv->SetVisAttributes(gasVis);
        G4VisAttributes detVis = new G4VisAttributes();
        detVis.SetColor(0, 0, 1);
        detVis.SetForceCloud(true);
        CubeDetectorLogic->SetVisAttributes(detVis);
        //G4OpticalSurface *OpSteelSurf = new G4OpticalSurface("SteelSurface", unified, polished, dielectric_metal);
        //OpSteelSurf->SetMaterialPropertiesTable(opticalprops::STEEL());


        G4OpticalSurface *opXenon_Glass = new G4OpticalSurface("DetectorSurface");
        opXenon_Glass->SetMaterialPropertiesTable(opticalprops::MgF2());
        opXenon_Glass->SetModel(glisur);                  // SetModel
        opXenon_Glass->SetType(dielectric_dielectric);   // SetType
        opXenon_Glass->SetFinish(ground);                 // SetFinish
        opXenon_Glass->SetPolish(0);

        new G4LogicalBorderSurface("DetectorSurface", GasXePlace, DetectorPlace, opXenon_Glass);

        Detector=CubeDetectorLogic;
        G4SDManager *Mang = G4SDManager::GetSDMpointer();
        PhotonSD *SD = new PhotonSD("PMT1");
        Mang->AddNewDetector(SD);
        Detector->SetSensitiveDetector(SD);

        cudaDeviceReset();
        G4CXOpticks::Get()->SetGeometry(physWorld);



        return physWorld;
    }

    void DetectorConstruction::ConstructSDandField() {


    }
}