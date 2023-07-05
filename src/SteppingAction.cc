//
// ********************************************************************
// * License and Disclaimer                                       	*
// *                                                              	*
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                         	*
// *                                                              	*
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.     	*
// *                                                              	*
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                  	*
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.      	*
// ********************************************************************
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "HistoManager.hh"

#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4HadronicProcess.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"




const int MAX_FILE_SIZE=3000000;

G4double ElectronEnergyArray[MAX_FILE_SIZE];
int NUMBEROFELECTRONS=0;
G4double PositronEnergyArray[MAX_FILE_SIZE];
int NUMBEROFPOSITRONS=0;
G4double NeutronEnergyArray[MAX_FILE_SIZE];
int NUMBEROFNEUTRONS=0;
G4double GammaEnergyArray[MAX_FILE_SIZE];
int NUMBEROFGAMMAS=0;

G4ThreeVector ElectronMomentumArray[MAX_FILE_SIZE];
G4ThreeVector PositronMomentumArray[MAX_FILE_SIZE];
G4ThreeVector NeutronMomentumArray[MAX_FILE_SIZE];
G4ThreeVector GammaMomentumArray[MAX_FILE_SIZE];

G4ThreeVector ElectronPositionArray[MAX_FILE_SIZE];
G4ThreeVector PositronPositionArray[MAX_FILE_SIZE];
G4ThreeVector NeutronPositionArray[MAX_FILE_SIZE];
G4ThreeVector GammaPositionArray[MAX_FILE_SIZE];

G4double ElectronEnergyArrayTarget[MAX_FILE_SIZE];
int NUMBEROFELECTRONSTARGET=0;
G4double PositronEnergyArrayTarget[MAX_FILE_SIZE];
int NUMBEROFPOSITRONSTARGET=0;
G4double NeutronEnergyArrayTarget[MAX_FILE_SIZE];
int NUMBEROFNEUTRONSTARGET=0;
G4double GammaEnergyArrayTarget[MAX_FILE_SIZE];
int NUMBEROFGAMMASTARGET=0;

G4ThreeVector ElectronMomentumArrayTarget[MAX_FILE_SIZE];
G4ThreeVector PositronMomentumArrayTarget[MAX_FILE_SIZE];
G4ThreeVector NeutronMomentumArrayTarget[MAX_FILE_SIZE];
G4ThreeVector GammaMomentumArrayTarget[MAX_FILE_SIZE];

G4ThreeVector ElectronPositionArrayTarget[MAX_FILE_SIZE];
G4double ElectronTimeArrayTarget[MAX_FILE_SIZE];
G4ThreeVector PositronPositionArrayTarget[MAX_FILE_SIZE];
G4double PositronTimeArrayTarget[MAX_FILE_SIZE];
G4ThreeVector NeutronPositionArrayTarget[MAX_FILE_SIZE];
G4double NeutronTimeArrayTarget[MAX_FILE_SIZE];
G4ThreeVector GammaPositionArrayTarget[MAX_FILE_SIZE];
G4double GammaTimeArrayTarget[MAX_FILE_SIZE];

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

 Run* run
   = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  // count processes
  //
  const G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4VProcess* process   =
               	const_cast<G4VProcess*>(endPoint->GetProcessDefinedStep());

	//G4cout << "Process Called:" << process->GetProcessName()<<G4endl;

  run->CountProcesses(process);


  // check that an real interaction occured (eg. not a transportation)
  G4StepStatus stepStatus = endPoint->GetStepStatus();
  G4bool transmit1 = stepStatus==fGeomBoundary;
  G4bool transmit2 = stepStatus==fWorldBoundary;


  bool kill=0;

  if(endPoint->GetGlobalTime()>0.011*ns)
  {
	kill=1;
	G4AnalysisManager* analysis = G4AnalysisManager::Instance();
	G4double energy = endPoint->GetKineticEnergy();
	G4double totenergy = endPoint->GetTotalEnergy();
	G4double charge = endPoint->GetCharge();
	G4ThreeVector momvec = endPoint->GetMomentum();
	G4ThreeVector position = endPoint->GetPosition();
	G4double beta = endPoint->GetBeta();
	G4double mom = sqrt(momvec[0]*momvec[0]+momvec[1]*momvec[1]+momvec[2]*momvec[2]);
	G4double posmag = sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]);
  	if(beta==1)
  	{
    	if(NUMBEROFGAMMASTARGET<MAX_FILE_SIZE)
    	{
      	GammaEnergyArrayTarget[NUMBEROFGAMMASTARGET]=energy;
      	GammaMomentumArrayTarget[NUMBEROFGAMMASTARGET]=momvec;
      	GammaPositionArrayTarget[NUMBEROFGAMMASTARGET]=position;
      	GammaTimeArrayTarget[NUMBEROFGAMMASTARGET]=endPoint->GetGlobalTime();
      	NUMBEROFGAMMASTARGET++;
    	}


  	}
  	else if(sqrt(totenergy*totenergy-mom*mom)>939&&sqrt(totenergy*totenergy-mom*mom)<940)
  	{
    	if(NUMBEROFNEUTRONSTARGET<MAX_FILE_SIZE)
    	{
      	NeutronEnergyArrayTarget[NUMBEROFNEUTRONSTARGET]=energy;
      	NeutronMomentumArrayTarget[NUMBEROFNEUTRONSTARGET]=momvec;
      	NeutronPositionArrayTarget[NUMBEROFNEUTRONSTARGET]=position;
      	NeutronTimeArrayTarget[NUMBEROFNEUTRONSTARGET]=endPoint->GetGlobalTime();
      	NUMBEROFNEUTRONSTARGET++;
    	}
  	}
  	else if(sqrt(totenergy*totenergy-mom*mom)>.5&&sqrt(totenergy*totenergy-mom*mom)<.52)
  	{
    	if(charge<0)
    	{
      	if(NUMBEROFELECTRONSTARGET<MAX_FILE_SIZE)
      	{
        	ElectronEnergyArrayTarget[NUMBEROFELECTRONSTARGET]=energy;
        	ElectronMomentumArrayTarget[NUMBEROFELECTRONSTARGET]=momvec;
        	ElectronPositionArrayTarget[NUMBEROFELECTRONSTARGET]=position;
        	ElectronTimeArrayTarget[NUMBEROFELECTRONSTARGET]=endPoint->GetGlobalTime();
        	NUMBEROFELECTRONSTARGET++;
      	}
  	}
    	else if(charge>0)
    	{
      	if(NUMBEROFPOSITRONSTARGET<MAX_FILE_SIZE)
      	{
        	PositronEnergyArrayTarget[NUMBEROFPOSITRONSTARGET]=energy;
        	PositronMomentumArrayTarget[NUMBEROFPOSITRONSTARGET]=momvec;
        	PositronPositionArrayTarget[NUMBEROFPOSITRONSTARGET]=position;
        	ElectronTimeArrayTarget[NUMBEROFPOSITRONSTARGET]=endPoint->GetGlobalTime();
        	NUMBEROFPOSITRONSTARGET++;
      	}
  	}
  	}

  	if(NUMBEROFGAMMASTARGET>.9*MAX_FILE_SIZE||NUMBEROFELECTRONSTARGET>.9*MAX_FILE_SIZE||NUMBEROFPOSITRONSTARGET>.9*MAX_FILE_SIZE||NUMBEROFNEUTRONSTARGET>.9*MAX_FILE_SIZE)
{
  PrintFiles();
}
  }
  //G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  if (transmit2)
  {
    if(kill==1)
    {

    fParticleFlag.clear();
    // kill event after first interaction
    //


    G4RunManager::GetRunManager()->AbortEvent();

    }
    return;
  }
  else if(transmit1)
  {
	G4AnalysisManager* analysis = G4AnalysisManager::Instance();
	G4double energy = endPoint->GetKineticEnergy();
	G4double totenergy = endPoint->GetTotalEnergy();
	G4double charge = endPoint->GetCharge();
	G4ThreeVector momvec = endPoint->GetMomentum();
	G4ThreeVector position = endPoint->GetPosition();
	G4double beta = endPoint->GetBeta();
	G4double mom = sqrt(momvec[0]*momvec[0]+momvec[1]*momvec[1]+momvec[2]*momvec[2]);
	G4double posmag = sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]);
  //     if(posmag<300*mm&&endPoint->GetMaterial()->GetName()=="interGalactic")
  //     {
    // 	if(beta==1)
    // 	{
  //   	// if(NUMBEROFGAMMASTARGET<MAX_FILE_SIZE)
  //   	// {
  //   	//   GammaEnergyArrayTarget[NUMBEROFGAMMASTARGET]=energy;
  //   	//   GammaMomentumArrayTarget[NUMBEROFGAMMASTARGET]=momvec;
  //   	//   GammaPositionArrayTarget[NUMBEROFGAMMASTARGET]=position;
  //   	//   NUMBEROFGAMMASTARGET++;
  //   	// }
  //
  //   	std::ofstream myfilegEn;
  //   	std::ofstream myfilegMom;
  //   	std::ofstream myfilegPos;
  //
  //   	myfilegEn.open("gammaenergytarget.txt",std::ios::app);
  //   	myfilegMom.open("gammamomentumtarget.txt",std::ios::app);
  //   	myfilegPos.open("gammapositiontarget.txt",std::ios::app);
  //
  //     	myfilegEn<< energy<<"\n";
  //     	myfilegMom<< momvec<<"\n";
  //     	myfilegPos<< position<<"\n";
  //
  // 	myfilegEn.close();
  // 	myfilegMom.close();
  // 	myfilegPos.close();
    // 	}
    // 	else if(sqrt(totenergy*totenergy-mom*mom)>939&&sqrt(totenergy*totenergy-mom*mom)<940)
    // 	{
  //   	if(NUMBEROFNEUTRONSTARGET<MAX_FILE_SIZE)
  //   	{
  //     	NeutronEnergyArrayTarget[NUMBEROFNEUTRONSTARGET]=energy;
  //     	NeutronMomentumArrayTarget[NUMBEROFNEUTRONSTARGET]=momvec;
  //     	NeutronPositionArrayTarget[NUMBEROFNEUTRONSTARGET]=position;
  //     	NeutronTimeArrayTarget[NUMBEROFNEUTRONSTARGET]=endPoint->GetGlobalTime();
  //     	NUMBEROFNEUTRONSTARGET++;
  //   	}
    // 	}
    // 	else if(sqrt(totenergy*totenergy-mom*mom)>.5&&sqrt(totenergy*totenergy-mom*mom)<.52)
    // 	{
    //		 if(charge<0)
    //		 {
  //     	if(NUMBEROFELECTRONSTARGET<MAX_FILE_SIZE)
  //     	{
  //       	ElectronEnergyArrayTarget[NUMBEROFELECTRONSTARGET]=energy;
  //       	ElectronMomentumArrayTarget[NUMBEROFELECTRONSTARGET]=momvec;
  //       	ElectronPositionArrayTarget[NUMBEROFELECTRONSTARGET]=position;
  //       	NUMBEROFELECTRONSTARGET++;
  //     	}
    //       }
    //		 else if(charge>0)
    //		 {
  //     	if(NUMBEROFPOSITRONSTARGET<MAX_FILE_SIZE)
  //     	{
  //       	PositronEnergyArrayTarget[NUMBEROFPOSITRONSTARGET]=energy;
  //       	PositronMomentumArrayTarget[NUMBEROFPOSITRONSTARGET]=momvec;
  //       	PositronPositionArrayTarget[NUMBEROFPOSITRONSTARGET]=position;
  //       	NUMBEROFPOSITRONSTARGET++;
  //     	}
    //       }
    // 	}
    // }
  if(kill==1)
  {

  fParticleFlag.clear();
  // kill event after first interaction
  //


  G4RunManager::GetRunManager()->AbortEvent();

}

  return;
  }






  //real processes : sum track length
  //
  G4double stepLength = aStep->GetStepLength();
  run->SumTrack(stepLength);

  //energy-momentum balance initialisation
  //
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4double Q         	= - prePoint->GetKineticEnergy();
  G4ThreeVector Pbalance = - prePoint->GetMomentum();

  //initialisation of the nuclear channel identification
  //
  G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
  G4String partName = particle->GetParticleName();


  G4String nuclearChannel = partName;
  G4HadronicProcess* hproc = dynamic_cast<G4HadronicProcess*>(process);
  const G4Isotope* target = NULL;
  if (hproc) target = hproc->GetTargetIsotope();
  G4String targetName = "XXXX";
  if (target) targetName = target->GetName();
  nuclearChannel += " + " + targetName + " --> ";
  if (targetName == "XXXX") run->SetTargetXXX(true);

  //scattered primary particle (if any)
  //
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  G4int ih = 1;
  if (aStep->GetTrack()->GetTrackStatus() == fAlive) {
	G4double energy = endPoint->GetKineticEnergy();
	analysis->FillH1(ih,energy);
	//
	G4ThreeVector momentum = endPoint->GetMomentum();
	Q    	+= energy;
	Pbalance += momentum;
	//
	nuclearChannel += partName + " + ";


  }

  //secondaries

  const std::vector<const G4Track*>* secondary
                                	= aStep->GetSecondaryInCurrentStep();



//G4cout << "Size:"<<secondary->size()<<G4endl;


if(secondary->size()>3)
{
  //kill=0;
  for(int i=0;i<secondary->size();i++)
  {
	if((*secondary)[i]->GetDefinition()->GetParticleName()=="neutron")
	{
  	//kill=0;
	}
  }
}
if(secondary->size()==0)
{
  //G4cout << "Particle:" << partName<<G4endl;
}
if(secondary->size()>0)
{
  //G4cout << "This is the reaction:" << partName << "--->" <<(*secondary)[0]->GetDefinition()->GetParticleName()<<G4endl;
  if(partName=="e-")
  {
	//kill=1;
  }
//kill=1;
}



  for (size_t lp=0; lp<(*secondary).size(); lp++) {
	particle = (*secondary)[lp]->GetDefinition();
	G4String name   = particle->GetParticleName();


	G4String type   = particle->GetParticleType();
	G4double energy = (*secondary)[lp]->GetKineticEnergy();
	run->ParticleCount(name,energy);
	//energy spectrum
	ih = 0;
     	if (particle == G4Gamma::Gamma())   	ih = 2;
	else if (particle == G4Electron::Electron()) ih = 3;
	else if (particle == G4Neutron::Neutron())   ih = 4;
	else if (particle == G4Proton::Proton()) 	ih = 5;
	else if (particle == G4Deuteron::Deuteron()) ih = 6;
	else if (particle == G4Alpha::Alpha())   	ih = 7;
	else if (type == "nucleus")              	ih = 8;
	else if (type == "meson")                	ih = 9;
	else if (type == "baryon")               	ih = 10;
	//if (ih == 4) analysis->FillH1(ih,energy);
	//atomic mass
	if (type == "nucleus") {
  	G4int A = particle->GetAtomicMass();
  	analysis->FillH1(13, A);
	}
	//energy-momentum balance
	G4ThreeVector momentum = (*secondary)[lp]->GetMomentum();
	Q    	+= energy;
	Pbalance += momentum;
	//count e- from internal conversion together with gamma
	if (particle == G4Electron::Electron()) particle = G4Gamma::Gamma();
	//particle flag
	fParticleFlag[particle]++;
  }


  //energy-momentum balance
  G4double Pbal = Pbalance.mag();
  run->Balance(Pbal);
  ih = 11;
  analysis->FillH1(ih,Q);
  ih = 12;
  analysis->FillH1(ih,Pbal);

  // nuclear channel
  const G4int kMax = 16;
  const G4String conver[] = {"0","","2 ","3 ","4 ","5 ","6 ","7 ","8 ","9 ",
                         	"10 ","11 ","12 ","13 ","14 ","15 ","16 "};
  std::map<G4ParticleDefinition*,G4int>::iterator ip;
  for (ip = fParticleFlag.begin(); ip != fParticleFlag.end(); ip++) {
	particle = ip->first;
	G4String name = particle->GetParticleName();
	G4int nb = ip->second;
	if (nb > kMax) nb = kMax;
	G4String Nb = conver[nb];
	if (particle == G4Gamma::Gamma()) {
 	run->CountGamma(nb);
 	name = "gamma or e-";


	}
	if (ip != fParticleFlag.begin()) nuclearChannel += " + ";
	nuclearChannel += Nb + name;
  }


  ///G4cout << "\n nuclear channel: " << nuclearChannel << G4endl;
  run->CountNuclearChannel(nuclearChannel, Q);


  G4double totenergy = endPoint->GetTotalEnergy();
  G4ThreeVector momvec = endPoint->GetMomentum();
  G4double mom = sqrt(momvec[0]*momvec[0]+momvec[1]*momvec[1]+momvec[2]*momvec[2]);
  G4double charge = endPoint->GetCharge();

  //if(sqrt(totenergy*totenergy-mom*mom)<.5||sqrt(totenergy*totenergy-mom*mom)>.52||charge>=0)
  if(kill==1)
  {

  fParticleFlag.clear();
  // kill event after first interaction
  //


  G4RunManager::GetRunManager()->AbortEvent();

}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::PrintFiles()
{

  std::ofstream myfileeEn;
  myfileeEn.open("electronenergy.txt",std::ios::app);
  std::ofstream myfileeMom;
  myfileeMom.open("electronmomentum.txt",std::ios::app);
  std::ofstream myfileePos;
  myfileePos.open("electronposition.txt",std::ios::app);
  std::ofstream myfileeTime;

  for (int i=0;i<NUMBEROFELECTRONS;i++)
  {
	myfileeEn<< ElectronEnergyArray[i]<<"\n";
	myfileeMom<< ElectronMomentumArray[i]<<"\n";
	myfileePos<< ElectronPositionArray[i]<<"\n";
  }

myfileeEn.close();
myfileeMom.close();
myfileePos.close();


  std::ofstream myfilepEn;
  myfilepEn.open("positronenergy.txt",std::ios::app);
  std::ofstream myfilepMom;
  myfilepMom.open("positronmomentum.txt",std::ios::app);
  std::ofstream myfilepPos;
  myfilepPos.open("positronposition.txt",std::ios::app);
std::ofstream myfilepTime;

  for (int i=0;i<NUMBEROFPOSITRONS;i++)
  {
	myfilepEn<< PositronEnergyArray[i]<<"\n";
	myfilepMom<< PositronMomentumArray[i]<<"\n";
	myfilepPos<< PositronPositionArray[i]<<"\n";
  }

myfilepEn.close();
myfilepMom.close();
myfilepPos.close();


  std::ofstream myfilenEn;
  myfilenEn.open("neutronenergy.txt",std::ios::app);
  std::ofstream myfilenMom;
  myfilenMom.open("neutronmomentum.txt",std::ios::app);
  std::ofstream myfilenPos;
  myfilenPos.open("neutronposition.txt",std::ios::app);
  std::ofstream myfilenTime;

  for (int i=0;i<NUMBEROFNEUTRONS;i++)
  {
	myfilenEn<< NeutronEnergyArray[i]<<"\n";
	myfilenMom<< NeutronMomentumArray[i]<<"\n";
	myfilenPos<< NeutronPositionArray[i]<<"\n";
  }

myfilenEn.close();
myfilenMom.close();
myfilenPos.close();


   std::ofstream myfileGEn;
//   myfileGEn.open("gammaenergy.txt",std::ios::app);
   std::ofstream myfileGMom;
//   myfileGMom.open("gammamomentum.txt",std::ios::app);
   std::ofstream myfileGPos;
//   myfileGPos.open("gammaposition.txt",std::ios::app);
   std::ofstream myfileGTime;
//
//   for (int i=0;i<NUMBEROFGAMMAS;i++)
//   {
// 	myfileGEn<< GammaEnergyArray[i]<<"\n";
// 	myfileGMom<< GammaMomentumArray[i]<<"\n";
// 	myfileGPos<< GammaPositionArray[i]<<"\n";
//   }
//
// myfileGEn.close();
// myfileGMom.close();
// myfileGPos.close();



  myfileeEn.open("electronenergytarget.txt",std::ios::app);
  myfileeMom.open("electronmomentumtarget.txt",std::ios::app);
  myfileePos.open("electronpositiontarget.txt",std::ios::app);
  myfileeTime.open("electrontimetarget.txt",std::ios::app);

  for (int i=0;i<NUMBEROFELECTRONSTARGET;i++)
  {
	myfileeEn<< ElectronEnergyArrayTarget[i]<<"\n";
	myfileeMom<< ElectronMomentumArrayTarget[i]<<"\n";
	myfileePos<< ElectronPositionArrayTarget[i]<<"\n";
	myfileeTime<< ElectronTimeArrayTarget[i]<<"\n";
  }

myfileeEn.close();
myfileeMom.close();
myfileePos.close();
myfileeTime.close();


  myfilepEn.open("positronenergytarget.txt",std::ios::app);
  myfilepMom.open("positronmomentumtarget.txt",std::ios::app);
  myfilepPos.open("positronpositiontarget.txt",std::ios::app);
  myfilepTime.open("positrontimetarget.txt",std::ios::app);

  for (int i=0;i<NUMBEROFPOSITRONSTARGET;i++)
  {
	myfilepEn<< PositronEnergyArrayTarget[i]<<"\n";
	myfilepMom<< PositronMomentumArrayTarget[i]<<"\n";
	myfilepPos<< PositronPositionArrayTarget[i]<<"\n";
	myfilepTime<< PositronTimeArrayTarget[i]<<"\n";

  }

myfilepEn.close();
myfilepMom.close();
myfilepPos.close();
myfilepTime.close();


  myfilenEn.open("neutronenergytarget.txt",std::ios::app);
  myfilenMom.open("neutronmomentumtarget.txt",std::ios::app);
  myfilenPos.open("neutronpositiontarget.txt",std::ios::app);
  myfilenTime.open("neutrontimetarget.txt",std::ios::app);

  for (int i=0;i<NUMBEROFNEUTRONSTARGET;i++)
  {
	myfilenEn<< NeutronEnergyArrayTarget[i]<<"\n";
	myfilenMom<< NeutronMomentumArrayTarget[i]<<"\n";
	myfilenPos<< NeutronPositionArrayTarget[i]<<"\n";
	myfilenTime<< NeutronTimeArrayTarget[i]<<"\n";
  }

myfilenEn.close();
myfilenMom.close();
myfilenPos.close();
myfilenTime.close();


  myfileGEn.open("gammaenergytarget.txt",std::ios::app);
  myfileGMom.open("gammamomentumtarget.txt",std::ios::app);
  myfileGPos.open("gammapositiontarget.txt",std::ios::app);
  myfileGTime.open("gammatimetarget.txt",std::ios::app);

  for (int i=0;i<NUMBEROFGAMMASTARGET;i++)
  {
	myfileGEn<< GammaEnergyArrayTarget[i]<<"\n";
	myfileGMom<< GammaMomentumArrayTarget[i]<<"\n";
	myfileGPos<< GammaPositionArrayTarget[i]<<"\n";
	myfileGTime<< GammaTimeArrayTarget[i]<<"\n";
  }

myfileGEn.close();
myfileGMom.close();
myfileGPos.close();
myfileGTime.close();


NUMBEROFELECTRONS=0;
NUMBEROFPOSITRONS=0;
NUMBEROFNEUTRONS=0;
NUMBEROFGAMMAS=0;

NUMBEROFELECTRONSTARGET=0;
NUMBEROFPOSITRONSTARGET=0;
NUMBEROFNEUTRONSTARGET=0;
NUMBEROFGAMMASTARGET=0;

}
