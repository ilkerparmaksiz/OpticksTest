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
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
using namespace CLHEP;
G4VPhysicalVolume* DetectorConstruction::Construct()
{


  // OpticalProperties
    constexpr G4double optPhotMinE_ =  0.2  * eV;
    constexpr G4double optPhotMaxE_ = 11.5  * eV;
    constexpr G4double noAbsLength_ = 1.e8  * m;

    // Constant that allows to convert nm to eV:
    // nm_to_eV_ / wavelength (nm) = energy (eV)
    constexpr G4double nm_to_eV_ = h_Planck * c_light * 1.e6;
    G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();

    // REFRACTIVE INDEX
    const G4int ri_entries = 200;
    G4double eWidth = (optPhotMaxE_ - optPhotMinE_) / ri_entries;

    std::vector<G4double> ri_energy;
    for (int i=0; i<ri_entries; i++) {
        ri_energy.push_back(optPhotMinE_ + i * eWidth);
    }

    std::vector<G4double> rIndex;
    for (int i=0; i<ri_entries; i++) {
        G4double wl = h_Planck * c_light / ri_energy[i] * 1000; // in micron
        // From refractiveindex.info
        rIndex.push_back(1 + 0.012055*(0.2075*pow(wl,2)/(91.012*pow(wl,2)-1) +
                                       0.0415*pow(wl,2)/(87.892*pow(wl,2)-1) +
                                       4.3330*pow(wl,2)/(214.02*pow(wl,2)-1)));
        //G4cout << "* GAr rIndex:  " << std::setw(5) << ri_energy[i]/eV
        //       << " eV -> " << rIndex[i] << G4endl;
    }
    mpt->AddProperty("RINDEX", ri_energy, rIndex);

    // ABSORPTION LENGTH
    std::vector<G4double> abs_energy = {optPhotMinE_, optPhotMaxE_};
    std::vector<G4double> absLength  = {noAbsLength_, noAbsLength_};
    mpt->AddProperty("ABSLENGTH", abs_energy, absLength);

    // EMISSION SPECTRUM
    G4double Wavelength_peak  = 128.000 * nm;
    G4double Wavelength_sigma =   2.929 * nm;
    G4double Energy_peak  = (h_Planck*c_light / Wavelength_peak);
    G4double Energy_sigma = (h_Planck*c_light * Wavelength_sigma / pow(Wavelength_peak,2));
    //G4cout << "*** GAr Energy_peak: " << Energy_peak/eV << " eV   Energy_sigma: "
    //       << Energy_sigma/eV << " eV" << G4endl;

    // Sampling from ~110 nm to 150 nm <----> from ~11.236 eV to 8.240 eV
    const G4int sc_entries = 380;
    std::vector<G4double> sc_energy;
    std::vector<G4double> intensity;
    for (int i=0; i<sc_entries; i++){
        sc_energy.push_back(8.240*eV + 0.008*i*eV);
        intensity.push_back(exp(-pow(Energy_peak/eV-sc_energy[i]/eV,2) /
                                (2*pow(Energy_sigma/eV, 2)))/(Energy_sigma/eV*sqrt(pi*2.)));
        //G4cout << "* GAr energy: " << std::setw(6) << sc_energy[i]/eV << " eV  ->  "
        //       << std::setw(6) << intensity[i] << G4endl;
    }
    mpt->AddProperty("SCINTILLATIONCOMPONENT1", sc_energy, intensity);
    mpt->AddProperty("SCINTILLATIONCOMPONENT2", sc_energy, intensity);
    mpt->AddProperty("ELSPECTRUM"             , sc_energy, intensity, 1);

    // CONST PROPERTIES
    mpt->AddConstProperty("SCINTILLATIONYIELD", 10/MeV);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1",   6.*ns);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT2",   3480.*ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", .136);
    mpt->AddConstProperty("SCINTILLATIONYIELD2", .864);

    mpt->AddProperty("FASTCOMPONENT",   sc_energy,intensity,1);
    mpt->AddProperty("SLOWCOMPONENT",  sc_energy,intensity,1);
    mpt->AddProperty("REEMISSIONPROB",  sc_energy,intensity,1);
    mpt->AddConstProperty("RESOLUTIONSCALE",    1.0);

    // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
    G4String name = "GasAr";

    G4Material* env_mat = G4Material::GetMaterial(name, false);

    if (env_mat == 0) {
        G4NistManager* nist = G4NistManager::Instance();

        G4double density = 1.60279*kg/m3;
        double pressure=10*bar;
        if (pressure/bar > 0.9 && pressure/bar < 1.1)
            density = 1.60279*kg/m3;
        else if (pressure/bar > 1.9 && pressure/bar < 2.1)
            density = 3.20719*kg/m3;
        else if (pressure/bar > 4.9 && pressure/bar < 5.1)
            density = 8.032*kg/m3;
        else if (pressure/bar > 9.9 && pressure/bar < 10.1)
            density = 16.1118*kg/m3;
        else if (pressure/bar > 14.9 && pressure/bar < 15.1)
            density = 24.2369 *kg/m3;
        else if (pressure/bar > 19.9 && pressure/bar < 20.1)
            density = 32.4066*kg/m3;
        else if (pressure/bar > 29.9 && pressure/bar < 30.1)
            density = 48.8708*kg/m3;
        else if (pressure/bar > 39.9 && pressure/bar < 40.1)
            density = 65.494*kg/m3;
        else
            G4Exception("[ArgonProperties]", "ArgonDensity()", FatalException,
                        "Unknown argon density for this pressure!");
        env_mat = new G4Material(name, density, 1,
                             kStateGas, 293*keV, pressure);

        G4Element* Ar = nist->FindOrBuildElement("Ar");

        env_mat->AddElement(Ar,1);
    }
     // Set  Gas Argon Materialies table
     env_mat->SetMaterialPropertiesTable(mpt);



    // Composition ranges correspond to stainless steel grade 304L
    name = "STEEL";
    G4Material* steelmat = G4Material::GetMaterial(name, false);

    if (steelmat == 0) {

        steelmat = new G4Material(name, 8000*kg/m3, 4);

        G4NistManager* nist = G4NistManager::Instance();

        G4Element* Fe = nist->FindOrBuildElement("Fe");
        steelmat->AddElement(Fe, 0.66);

        G4Element* Cr = nist->FindOrBuildElement("Cr");
        steelmat->AddElement(Cr, 0.20);

        G4Element* Mn = nist->FindOrBuildElement("Mn");
        steelmat->AddElement(Mn, 0.02);

        G4Element* Ni = nist->FindOrBuildElement("Ni");
        steelmat->AddElement(Ni, 0.12);
    }

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
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
  auto solidEnv = new G4Box("GasAr_Solid",                    // its name
    0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ);  // its size

  auto logicEnv = new G4LogicalVolume(solidEnv,env_mat,"GasAr_Logic");


    auto Steel_Cover_solid = new G4Box("SteelCover_Solid",                    // its name
                              0.5 * env_sizeXY+5*CLHEP::mm, 0.5 * env_sizeXY+10*CLHEP::mm, 0.5 * env_sizeZ+5*CLHEP::mm);  // its size
    auto Steel_Cover_logic = new G4LogicalVolume(Steel_Cover_solid,  // its solid
                                                 steelmat,                                     // its material
                                        "SteelCover_logic");                                 // its name
  //
  //always return the physical World
  //
    auto SteelPlace=new G4PVPlacement(0,G4ThreeVector (),"SteelCover",Steel_Cover_logic,physWorld,0,0,0);
    auto GasArPlace=new G4PVPlacement(0,G4ThreeVector (),"GasAr",logicEnv,SteelPlace,0,0,0);

  std::cout <<"Setting our detector geometry with opticks" <<std::endl;
  //std::cout << U4Physics::LEVEL <<std::endl;
  //;

  G4CXOpticks::Get()->SetGeometry(physWorld);


    std::cout << SEventConfig::Desc() <<std::endl;
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
