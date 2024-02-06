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
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
// #include "Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4PVPlacement.hh"
#include "U4Recorder.hh"
#include "G4CXOpticks.hh"
#include "SEvt.hh"

#include "G4AnalysisManager.hh"
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{

  // Analysis Manager
  auto ana = G4AnalysisManager::Instance();
    // Create Directories
    ana->SetVerboseLevel(1);
    ana->SetNtupleMerging(true);
    // Create NTuple
    ana->CreateNtuple("G4Photons","GEANT4 Photon Info");
    ana->CreateNtupleIColumn("Event");
    ana->CreateNtupleDColumn("fx");
    ana->CreateNtupleDColumn("fy");
    ana->CreateNtupleDColumn("fz");
    ana->CreateNtupleDColumn("ft");
    ana->CreateNtupleSColumn("Volume");
    ana->CreateNtupleSColumn("Process");
    ana->FinishNtuple();

    ana->CreateNtuple("particles","Particle Names");
    ana->CreateNtupleIColumn("Event");
    ana->CreateNtupleIColumn("ID");
    ana->CreateNtupleSColumn("name");
    ana->CreateNtupleDColumn("x");
    ana->CreateNtupleDColumn("y");
    ana->CreateNtupleDColumn("z");
    ana->CreateNtupleDColumn("t");
    ana->CreateNtupleSColumn("Volume");

    ana->CreateNtuple("Opticks","Opticks Photon Hits");
    ana->CreateNtupleIColumn("Event");
    ana->CreateNtupleIColumn("id");
    ana->CreateNtupleFColumn("x");
    ana->CreateNtupleFColumn("y");
    ana->CreateNtupleFColumn("z");
    ana->CreateNtupleFColumn("t");
    ana->CreateNtupleFColumn("mx");
    ana->CreateNtupleFColumn("my");
    ana->CreateNtupleFColumn("mz");
    ana->CreateNtupleFColumn("px");
    ana->CreateNtupleFColumn("py");
    ana->CreateNtupleFColumn("pz");
    ana->CreateNtupleFColumn("wavelength");





    ana->FinishNtuple();
    ana->SetNtupleActivation(1);





}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
  // inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(false);


 auto ana=G4AnalysisManager::Instance();

 G4String FileName ="Opticks_CPU.root";
 ana->OpenFile(FileName);
 G4cout << "Will be Writing to file " << FileName << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  // Initialize
  auto ana = G4AnalysisManager::Instance();

  // Save into File
  //ana->SetCompressionLevel(3);
  ana->Write();
  ana->CloseFile();
  ana->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
