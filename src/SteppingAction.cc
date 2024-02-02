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
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalPhoton.hh"
#include "U4.hh"
#include "U4VPrimaryGenerator.h"
#include "SEvt.hh"
#include "G4Ions.hh"
#include "G4Alpha.hh"
#include "G4CXOpticks.hh"
namespace {G4Mutex opticks_mt =G4MUTEX_INITIALIZER;}
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detConstruction->GetScoringVolume();
  }

  // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
  auto atrack=step->GetTrack();

    if(atrack->GetParticleDefinition()==G4Alpha::Definition()){

        // Opticks
      U4::CollectGenstep_DsG4Scintillation_r4695(atrack,step,3000,0,4*ns);
       //SEvt::AddTorchGenstep();
        // Opticks


        /*
        G4cout<<" Opticks End of Event Action" <<G4endl;
        G4AutoLock lock(&opticks_mt);
        G4CXOpticks * g4cx=G4CXOpticks::Get();
        G4RunManager* rm     = G4RunManager::GetRunManager();
        const G4Event* event = rm->GetCurrentEvent();
        G4int eventID        = event->GetEventID();
        G4int ngenstep=SEvt::GetNumGenstepFromGenstep(eventID);
        G4int nphotons=SEvt::GetNumPhotonCollected(eventID);

        G4cout << "Number of Steps Generated " <<ngenstep << G4endl;
        G4cout << "Number of Photons Generated " <<nphotons << G4endl;

        // Simulate the photons
        if(nphotons>0 and ngenstep>0){
            std::cout<<g4cx->desc()<<std::endl;
            std::cout<<"--- G4Optickx ---" << g4cx->descSimulate() <<std::endl;
            g4cx->simulate(eventID,1); // For Simulation
            cudaDeviceSynchronize();
            //g4cx->render();  // For Rendering
        }
        G4int nHits=SEvt::GetNumHit(0);
        G4cout << "Number of Hits  " <<nHits << G4endl;
        */
      //SEvt::AddTorchGenstep();

      //U4::CollectGenstep_G4Cerenkov_modified(atrack,step,2000,1,0,3.14,1,0,10,2000);


    }



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
