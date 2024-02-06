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
/// \file B1/src/EventAction.cc
/// \brief Implementation of the B1::EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#  include "SEvt.hh"
#  include "NP.hh"
#  include "G4CXOpticks.hh"
#include "U4Recorder.hh"
#include "PhotonSD.hh"
#include "U4HitGet.h"
#include "G4SDManager.hh"
namespace {G4Mutex opticks_mt =G4MUTEX_INITIALIZER;}
namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: fRunAction(runAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  // Opticks
  auto fRecord=U4Recorder::Get();
    fRecord->BeginOfEventAction_(evt->GetEventID());
    auto fSEvt=SEvt::Get(evt->GetEventID());
    //fSEvt->beginOfEvent(evt->GetEventID());


  fEdep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
    // Opticks
    auto fSEvt=SEvt::Get(evt->GetEventID());
   // fSEvt->endOfEvent(evt->GetEventID());
    auto fRecord=U4Recorder::Get();
    fRecord->EndOfEventAction_(evt->GetEventID());

  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);

   G4cout<<" Opticks End of Event Action" <<G4endl;
    G4AutoLock lock(&opticks_mt);
    G4CXOpticks * g4cx=G4CXOpticks::Get();
    G4int eventID=evt->GetEventID();
    G4int ngenstep=SEvt::GetNumGenstepFromGenstep(eventID);
    G4int nphotons=SEvt::GetNumPhotonCollected(eventID);



    // Simulate the photons
      if(nphotons>0 and ngenstep>0){
          std::cout<<g4cx->desc()<<std::endl;
          std::cout<<"--- G4Optickx ---" << g4cx->descSimulate() <<std::endl;
          g4cx->simulate(eventID,0); // For Simulation
          cudaDeviceSynchronize();
          //g4cx->render();  // For Rendering
      }


    PhotonSD* SD = (PhotonSD*) G4SDManager::GetSDMpointer()->FindSensitiveDetector("PMT1");
    if(SEvt::GetNumHit(eventID)>0){
        SD->OpticksHits();

    }

    G4cout << "Number of Steps Generated " <<ngenstep << G4endl;
    G4cout << "Number of Photons Generated " <<nphotons << G4endl;
    G4cout << "Number of Hits Opticks  " <<SEvt::GetNumHit(eventID)<< G4endl;
    G4cout << "Number of Hits GEANT4  " <<SD->GetGEANT4Hits()<< G4endl;
    G4CXOpticks::Get()->reset(eventID);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
