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
#include "G4AnalysisManager.hh"
#include "SEventConfig.hh"
#include "G4Scintillation.hh"
#include "G4EventManager.hh"

namespace {G4Mutex opticks_mt =G4MUTEX_INITIALIZER;}
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{

    Counter=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

    auto ana=G4AnalysisManager::Instance();
    auto run=G4RunManager::GetRunManager();

  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
  G4String VolumeName="None";
  if(volume->GetName()!="") VolumeName=volume->GetName();



  auto atrack=step->GetTrack();
  if(atrack->GetParticleDefinition()!=G4OpticalPhoton::Definition())  {

      ana->FillNtupleIColumn(1,0,run->GetCurrentEvent()->GetEventID());
      ana->FillNtupleIColumn(1,1,atrack->GetTrackID());
      ana->FillNtupleSColumn(1,2,atrack->GetParticleDefinition()->GetParticleName());
      ana->FillNtupleDColumn(1,3,atrack->GetPosition()[0]/mm);
      ana->FillNtupleDColumn(1,4,atrack->GetPosition()[1]/mm);
      ana->FillNtupleDColumn(1,5,atrack->GetPosition()[2]/mm);
      ana->FillNtupleDColumn(1,6,atrack->GetProperTime()/ns);
      ana->FillNtupleSColumn(1,7,VolumeName);
      ana->AddNtupleRow(1);
  }



        G4double edepStep = step->GetTotalEnergyDeposit();
        fEventAction->AddEdep(edepStep);
        // Opticks

        G4SteppingManager * sMg=G4EventManager::GetEventManager()->GetTrackingManager()->GetSteppingManager();
        G4StepStatus stepStatus=sMg->GetfStepStatus();
        if(stepStatus!=fAtRestDoItProc){
            G4ProcessVector * PostStepProc=sMg->GetfPostStepDoItVector();
            size_t MaxSteps=sMg->GetMAXofPostStepLoops();
            for (int stp=0;stp<MaxSteps;stp++){
                if((*PostStepProc)[stp]->GetProcessName()=="Scintillation"){
                    G4Scintillation *ScintProc= (G4Scintillation*) (*PostStepProc)[stp];
                    G4int num_photons=ScintProc->GetNumPhotons();


                    if(num_photons>0){
                        G4MaterialPropertiesTable * MPT = atrack->GetMaterial()->GetMaterialPropertiesTable();
                        G4double t1,t2=0;
                        G4int singlets,triplets=0;
                        t1=MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT1);
                        t2=MPT->GetConstProperty(kSCINTILLATIONTIMECONSTANT2);
                        singlets= floor(MPT->GetConstProperty(kSCINTILLATIONYIELD1)*num_photons);
                        triplets= ceil(MPT->GetConstProperty(kSCINTILLATIONYIELD2)*num_photons);
                        //std::cout << "Scintilation "<< num_photons <<" Amount of Singlets " <<singlets <<" Triplets " << triplets <<std::endl;
                        if(singlets>0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(atrack,step,singlets,0,t1);
                        if(triplets>0)
                            U4::CollectGenstep_DsG4Scintillation_r4695(atrack,step,triplets,1,t2);
                        //U4::CollectGenstep_DsG4Scintillation_r4695(atrack,step,num_photons,1,t2);
                    }

                }
            }
        }

        //U4::CollectGenstep_G4Cerenkov_modified(atrack,step,2000,1,0,3.14,1,0,2000,2000);

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
