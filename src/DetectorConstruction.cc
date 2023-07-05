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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4SolidStore.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4RunManager.hh"


#include "G4FieldManager.hh"



#include "F02ElectricFieldSetup.hh"


#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4TransportationManager.hh"
#include "G4IntegrationDriver.hh"
#include "G4ChordFinder.hh"

G4ElectricField* pEMfield;
G4EqMagElectricField* pEquation;
G4ChordFinder* pChordFinder;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fPBox(0), fLBox(0), fMaterial(0), fDetectorMessenger(0)
{
  fBoxSize = 10*m;
  DefineMaterials();
  SetMaterial("Molybdenum98");
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
 // define a Material from isotopes
 //
 MaterialWithSingleIsotope("Molybdenum98", "Mo98",  10.28*g/cm3, 42, 98);

 //NE213
 G4Element* H  = new G4Element("Hydrogen" ,"H" , 1.,  1.01*g/mole);
 G4Element* C  = new G4Element("Carbon"   ,"C" , 6., 12.00*g/mole);
 G4Material* ne213 =
 new G4Material("NE213", 0.874*g/cm3, 2);
 ne213->AddElement(H,    9.2*perCent);
 ne213->AddElement(C,   90.8*perCent);

 G4Material* hydrogen =
 new G4Material("hydrogen", 1.0*g/cm3, 1);
 hydrogen->AddElement(H, 1);

 G4Material* carbon =
 new G4Material("carbon", 1.0*g/cm3, 1);
 carbon->AddElement(C, 1);

 G4Material* plastic =
 new G4Material("plastic", 1.0*g/cm3, 2);
 plastic->AddElement(H, 1);
 plastic->AddElement(C, 1);

 // or use G4-NIST materials data base
 //
 G4NistManager* man = G4NistManager::Instance();
 man->FindOrBuildMaterial("G4_B");

 G4Element*  O = man->FindOrBuildElement("O");
 G4Element* Hf = man->FindOrBuildElement("Hf");

 G4Material* HfO2 = new G4Material("HfO2", 9.68*g/cm3, 2);
 HfO2->AddElement(Hf, 1);
 HfO2->AddElement(O , 2);

 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
  G4Material* material;
 if(name=="Platinum")
 {
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotopePt190 = new G4Isotope("Pt", 78, 190);
 G4Isotope* isotopePt192 = new G4Isotope("Pt", 78, 192);
 G4Isotope* isotopePt194 = new G4Isotope("Pt", 78, 194);
 G4Isotope* isotopePt195 = new G4Isotope("Pt", 78, 195);
 G4Isotope* isotopePt196 = new G4Isotope("Pt", 78, 196);
 G4Isotope* isotopePt198 = new G4Isotope("Pt", 78, 198);

 G4Element* element  = new G4Element("Platinum", "Pt", ncomponents=6);
 element->AddIsotope(isotopePt190, abundance= .012*perCent);
 element->AddIsotope(isotopePt192, abundance= .782*perCent);
 element->AddIsotope(isotopePt194, abundance= 32.86*perCent);
 element->AddIsotope(isotopePt195, abundance= 33.78*perCent);
 element->AddIsotope(isotopePt196, abundance= 25.21*perCent);
 element->AddIsotope(isotopePt198, abundance= 7.36*perCent);

 material = new G4Material(name, 21.45*g/cm3, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);
 }

 else if(name=="Gold")
 {
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotopeAu197 = new G4Isotope("Au", 79, 197);

 G4Element* element  = new G4Element("Gold", "Au", ncomponents=1);
 element->AddIsotope(isotopeAu197, abundance= 100*perCent);

 material = new G4Material(name, 19.32*g/cm3, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);
 }

 else if(name=="Hydrogen")
 {
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotopeH1 = new G4Isotope("H", 1, 1);

 G4Element* element  = new G4Element("Hydrogen", "H", ncomponents=1);
 element->AddIsotope(isotopeH1, abundance= 100*perCent);

 material = new G4Material(name, 1000*g/cm3, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);
 }

 else if(name=="Aluminum")
 {
 G4int ncomponents;
 G4double abundance, massfraction;


 G4Isotope* isotopeAl27 = new G4Isotope("Al", 13, 27);

 G4Element* element  = new G4Element("Aluminum", "Al", ncomponents=1);
 element->AddIsotope(isotopeAl27, abundance= 100*perCent);

 material = new G4Material(name, 19.32*g/cm3, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);


 }

 else if(name=="Rhenium"||name=="Uranium")
 {
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotopeRe185 = new G4Isotope("Re", 75, 185);
 G4Isotope* isotopeRe187 = new G4Isotope("Re", 75, 187);

 G4Element* element  = new G4Element("Rhenium", "Re", ncomponents=2);
 element->AddIsotope(isotopeRe185, abundance= 37.40*perCent);
 element->AddIsotope(isotopeRe187, abundance= 62.60*perCent);

 material = new G4Material(name, 19.32*g/cm3, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);
 }



 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

	G4NistManager *nist=G4NistManager::Instance();
	G4double energy[2]={1.239841939*eV/.9,1.239841939*eV/.2};
	G4double rindexWorld[2]={1, 1};

	G4double atomicNumber=1;
	G4double massOfMole=1.008*g/mole;
	G4double density=1.e-25*g/cm3;
	G4double temperature=2.73*kelvin;
	G4double pressure=3.e-18*pascal;

  G4Material *worldMat=new G4Material("interGalactic", atomicNumber,massOfMole, density, kStateGas, temperature, pressure);
  G4Material *EMMat=new G4Material("EMField", atomicNumber,massOfMole, density, kStateGas, temperature, pressure);

	G4MaterialPropertiesTable *mptWorld=new G4MaterialPropertiesTable();
	mptWorld->AddProperty("RINDEX", energy, rindexWorld, 2);


	worldMat->SetMaterialPropertiesTable(mptWorld);







 G4int ncomponents1;
 G4double abundance1, massfraction1;

 G4Isotope* isotopeC12 = new G4Isotope("C", 6, 12);
 G4Isotope* isotopeC13 = new G4Isotope("C", 6, 13);

 G4Element* elementC  = new G4Element("Carbon", "C", ncomponents1=2);
 elementC->AddIsotope(isotopeC12, abundance1= 98.9*perCent);
 elementC->AddIsotope(isotopeC13, abundance1= 1.1*perCent);

 G4Isotope* isotopeN14 = new G4Isotope("N", 7, 14);
 G4Isotope* isotopeN15 = new G4Isotope("N", 7, 15);

 G4Element* elementN  = new G4Element("Nitrogen", "N", ncomponents1=2);
 elementN->AddIsotope(isotopeN14, abundance1= 99.65*perCent);
 elementN->AddIsotope(isotopeN15, abundance1= 0.35*perCent);


 G4Isotope* isotopeAl27 = new G4Isotope("Al", 13, 27);

 G4Element* elementAl  = new G4Element("Aluminum", "Al", ncomponents1=1);
 elementAl->AddIsotope(isotopeAl27, abundance1= 100*perCent);


 G4Isotope* isotopeSi28 = new G4Isotope("Si", 14, 28);
 G4Isotope* isotopeSi29 = new G4Isotope("Si", 14, 29);
 G4Isotope* isotopeSi30 = new G4Isotope("Si", 14, 30);

 G4Element* elementSi  = new G4Element("Silicon", "Si", ncomponents1=3);
 elementSi->AddIsotope(isotopeSi28, abundance1= 92.23*perCent);
 elementSi->AddIsotope(isotopeSi29, abundance1= 4.67*perCent);
 elementSi->AddIsotope(isotopeSi30, abundance1= 3.10*perCent);


 G4Isotope* isotopeP31 = new G4Isotope("P", 15, 31);

 G4Element* elementP  = new G4Element("Phosphorus", "P", ncomponents1=1);
 elementP->AddIsotope(isotopeP31, abundance1= 100*perCent);


 G4Isotope* isotopeS32 = new G4Isotope("S", 16, 32);
 G4Isotope* isotopeS33 = new G4Isotope("S", 16, 33);
 G4Isotope* isotopeS34 = new G4Isotope("S", 16, 34);
 G4Isotope* isotopeS36 = new G4Isotope("S", 16, 36);

 G4Element* elementS  = new G4Element("Sulfur", "S", ncomponents1=4);
 elementS->AddIsotope(isotopeS32, abundance1= 95.02*perCent);
 elementS->AddIsotope(isotopeS33, abundance1= 0.75*perCent);
 elementS->AddIsotope(isotopeS34, abundance1= 4.21*perCent);
 elementS->AddIsotope(isotopeS36, abundance1= 0.02*perCent);


 G4Isotope* isotopeCr50 = new G4Isotope("Cr", 24, 50);
 G4Isotope* isotopeCr52 = new G4Isotope("Cr", 24, 52);
 G4Isotope* isotopeCr53 = new G4Isotope("Cr", 24, 53);
 G4Isotope* isotopeCr54 = new G4Isotope("Cr", 24, 54);

 G4Element* elementCr  = new G4Element("Chromium", "Cr", ncomponents1=4);
 elementCr->AddIsotope(isotopeCr50, abundance1= 4.35*perCent);
 elementCr->AddIsotope(isotopeCr52, abundance1= 83.79*perCent);
 elementCr->AddIsotope(isotopeCr53, abundance1= 9.5*perCent);
 elementCr->AddIsotope(isotopeCr54, abundance1= 2.36*perCent);


 G4Isotope* isotopeMn55 = new G4Isotope("Mn", 25, 55);

 G4Element* elementMn  = new G4Element("Manganese", "Mn", ncomponents1=1);
 elementMn->AddIsotope(isotopeMn55, abundance1= 100*perCent);


 G4Isotope* isotopeFe54 = new G4Isotope("Fe", 26, 54);
 G4Isotope* isotopeFe56 = new G4Isotope("Fe", 26, 56);
 G4Isotope* isotopeFe57 = new G4Isotope("Fe", 26, 57);
 G4Isotope* isotopeFe58 = new G4Isotope("Fe", 26, 58);

 G4Element* elementFe  = new G4Element("Iron", "Fe", ncomponents1=4);
 elementFe->AddIsotope(isotopeFe54, abundance1= 5.8*perCent);
 elementFe->AddIsotope(isotopeFe56, abundance1= 91.7*perCent);
 elementFe->AddIsotope(isotopeFe57, abundance1= 2.2*perCent);
 elementFe->AddIsotope(isotopeFe58, abundance1= .3*perCent);


 G4Isotope* isotopeNi58 = new G4Isotope("Ni", 28, 58);
 G4Isotope* isotopeNi60 = new G4Isotope("Ni", 28, 60);
 G4Isotope* isotopeNi61 = new G4Isotope("Ni", 28, 61);
 G4Isotope* isotopeNi62 = new G4Isotope("Ni", 28, 62);
 G4Isotope* isotopeNi64 = new G4Isotope("Ni", 28, 64);

 G4Element* elementNi  = new G4Element("Nickel", "Ni", ncomponents1=5);
 elementNi->AddIsotope(isotopeNi58, abundance1= 68.0769*perCent);
 elementNi->AddIsotope(isotopeNi60, abundance1= 26.2231*perCent);
 elementNi->AddIsotope(isotopeNi61, abundance1= 1.1399*perCent);
 elementNi->AddIsotope(isotopeNi62, abundance1= 3.6345*perCent);
 elementNi->AddIsotope(isotopeNi64, abundance1= 0.9256*perCent);






 G4Material* materialSt = new G4Material("Steel", 7.93*g/cm3, ncomponents1=9);
 materialSt->AddElement(elementFe, massfraction1=70.995*perCent);
 materialSt->AddElement(elementCr, massfraction1=18*perCent);
 materialSt->AddElement(elementNi, massfraction1=8*perCent);
 materialSt->AddElement(elementNi, massfraction1=2*perCent);
 materialSt->AddElement(elementMn, massfraction1=0.1*perCent);
 materialSt->AddElement(elementN, massfraction1=.03*perCent);
 materialSt->AddElement(elementS, massfraction1=.08*perCent);
 materialSt->AddElement(elementC, massfraction1=.75*perCent);
 materialSt->AddElement(elementSi, massfraction1=.045*perCent);


 G4Material* materialAl = new G4Material("Aluminum", 2.7*g/cm3, ncomponents1=1);
 materialAl->AddElement(elementAl, massfraction1=100.*perCent);








G4Tubs *solidWorld=new G4Tubs("Container",0*m,796.93*mm+30*mm,.5*m,0*deg,360*deg);
G4LogicalVolume *logicWorld=new G4LogicalVolume(solidWorld, worldMat, "logicWorld");
G4VPhysicalVolume *physWorld=new 					G4PVPlacement(0,G4ThreeVector(0.,0.,0.),logicWorld,"physWorld", 0, false, 0, true);



  G4ThreeVector AxisOfRotation=G4ThreeVector(0,1,0);
  G4ThreeVector AxisOfRotation1=G4ThreeVector(1,0,0);
  G4RotationMatrix * RotMat=new G4RotationMatrix();
  RotMat->rotate(90*deg,AxisOfRotation);
  RotMat->rotate(17*deg,AxisOfRotation1);



  // pEMfield = new G4UniformElectricField(
  //              G4ThreeVector(-5.0/.45*megavolt/cm,0.0,0.0));
  //
  // // Create an equation of motion for this field
  // pEquation = new G4EqMagElectricField(pEMfield);
  //
  // G4int nvar = 8;
  //
  // // Create the Runge-Kutta 'stepper' using the efficient 'DoPri5' method
  // auto pStepper = new G4DormandPrince745( pEquation, nvar );
  //
  // // Get the global field manager
  // auto fieldManager= G4TransportationManager::GetTransportationManager()->
  //      GetFieldManager();
  // // Set this field to the global field manager
  // fieldManager->SetDetectorField( pEMfield );
  //
  // G4double minStep     = 0.010*mm ; // minimal step of 10 microns
  //
  // // The driver will ensure that integration is control to give
  // //   acceptable integration error
  // auto pIntgrationDriver =
  //    new G4IntegrationDriver<G4DormandPrince745>(minStep,
  //                                                pStepper,
  //                                                nvar);
  //
  // pChordFinder = new G4ChordFinder(pIntgrationDriver);
  // fieldManager->SetChordFinder( pChordFinder );




   if (fMaterial->GetName()=="Gold"||fMaterial->GetName()=="Platinum")
   {

       G4Tubs*
       sBox = new G4Tubs("Container",0,5*mm,fBoxSize/2,0*deg,360*deg);   //its dimensions

       fLBox = new G4LogicalVolume(sBox,                     //its shape
                                  fMaterial,                 //its material
                                  fMaterial->GetName());     //its name




       G4bool allLocal =true;
       //fLBox->SetFieldManager(fieldManager,allLocal);

       fPBox = new G4PVPlacement(RotMat,                          //no rotation
                                 G4ThreeVector(62.47*mm+fBoxSize/2/std::cos(17*3.141592653589793/180),92.6*mm,0*cm),           //at (0,0,0)
                                 fLBox,                      //its logical volume
                                 fMaterial->GetName(),       //its name
                                 logicWorld,                          //its mother  volume
                                 false,                      //no boolean operation
                                 0);


      // G4Box*
      // sBoxEM = new G4Box("EMField",1*cm/2,1*cm/2,1*cm/2);   //its dimensions
      //
      // G4LogicalVolume* fLBoxEM = new G4LogicalVolume(sBoxEM,                     //its shape
      //                            EMMat,                 //its material
      //                            EMMat->GetName());     //its name
      //
      //
      //
      //
      // G4bool allLocal =true;
      // fLBoxEM->SetFieldManager(fieldManager,allLocal);
      //
      // G4PVPlacement* fPBoxEM = new G4PVPlacement(0,                          //no rotation
      //                           G4ThreeVector(62.47*mm+1*cm/2+7*mm,92.6*mm,0*cm),           //at (0,0,0)
      //                           fLBoxEM,                      //its logical volume
      //                           EMMat->GetName(),       //its name
      //                           logicWorld,                          //its mother  volume
      //                           false,                      //no boolean operation
      //                           0);



   }
  if (fMaterial->GetName()=="Hydrogen")
  {

        G4Sphere*
        sBox = new G4Sphere("Container",0,fBoxSize/2,0*deg,360*deg,0*deg,180*deg);   //its dimensions

        fLBox = new G4LogicalVolume(sBox,                     //its shape
                                   fMaterial,                 //its material
                                   fMaterial->GetName());     //its name




        G4bool allLocal =true;
        //fLBox->SetFieldManager(fieldManager,allLocal);

        fPBox = new G4PVPlacement(RotMat,                          //no rotation
                                  G4ThreeVector(62.47*mm,92.6*mm,0*cm),           //at (0,0,0)
                                  fLBox,                      //its logical volume
                                  fMaterial->GetName(),       //its name
                                  logicWorld,                          //its mother  volume
                                  false,                      //no boolean operation
                                  0);


       // G4Box*
       // sBoxEM = new G4Box("EMField",1*cm/2,1*cm/2,1*cm/2);   //its dimensions
       //
       // G4LogicalVolume* fLBoxEM = new G4LogicalVolume(sBoxEM,                     //its shape
       //                            EMMat,                 //its material
       //                            EMMat->GetName());     //its name
       //
       //
       //
       //
       // G4bool allLocal =true;
       // fLBoxEM->SetFieldManager(fieldManager,allLocal);
       //
       // G4PVPlacement* fPBoxEM = new G4PVPlacement(0,                          //no rotation
       //                           G4ThreeVector(62.47*mm+1*cm/2+7*mm,92.6*mm,0*cm),           //at (0,0,0)
       //                           fLBoxEM,                      //its logical volume
       //                           EMMat->GetName(),       //its name
       //                           logicWorld,                          //its mother  volume
       //                           false,                      //no boolean operation
       //                           0);



    }


 else if (fMaterial->GetName()=="Rhenium")
 {


  G4Box*
  sBox = new G4Box("Container",.5*cm,.5*cm,fBoxSize/2);   //its dimensions

  fLBox = new G4LogicalVolume(sBox,                     //its shape
                             fMaterial,                 //its material
                             fMaterial->GetName());     //its name

  fPBox = new G4PVPlacement(RotMat,                          //no rotation
                            G4ThreeVector(62.47*mm+fBoxSize/2/std::cos(17*3.141592653589793/180),92.6*mm,0*cm),           //at (0,0,0)
                            fLBox,                      //its logical volume
                            fMaterial->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);






 	G4Isotope* isotopeW180 = new G4Isotope("W", 74, 180);
 	G4Isotope* isotopeW182 = new G4Isotope("W", 74, 182);
 	G4Isotope* isotopeW183 = new G4Isotope("W", 74, 183);
 	G4Isotope* isotopeW184 = new G4Isotope("W", 74, 184);
 	G4Isotope* isotopeW186 = new G4Isotope("W", 74, 186);

 	G4Element* elementW  = new G4Element("Tungsten", "W", ncomponents1=5);
 	elementW->AddIsotope(isotopeW180, abundance1=  0.12*perCent);
 	elementW->AddIsotope(isotopeW182, abundance1= 26.50*perCent);
 	elementW->AddIsotope(isotopeW183, abundance1= 14.31*perCent);
 	elementW->AddIsotope(isotopeW184, abundance1= 30.64*perCent);
 	elementW->AddIsotope(isotopeW186, abundance1= 28.43*perCent);


 	G4Material* materialW = new G4Material("Tungsten", 19.3*g/cm3, ncomponents1=1);
 	materialW->AddElement(elementW, massfraction1=100.*perCent);

  	G4Box*
  	sBoxW = new G4Box("Container",.5*cm,.5*cm,.5*cm);   //its dimensions

  	G4Tubs*
  	sBoxDivot = new G4Tubs("Divot",0,1.1*mm,2*mm,0*deg,360*deg);   //its dimensions

  	G4RotationMatrix* RotMatNull=new G4RotationMatrix();
  	G4ThreeVector* TransVec=new G4ThreeVector(0,0,3*mm);

  	G4SubtractionSolid* subBoxW=new G4SubtractionSolid("Containerw/Divot", sBoxW, sBoxDivot, RotMatNull, *TransVec);

  	G4LogicalVolume* fLBoxW = new G4LogicalVolume(subBoxW,                     //its shape
                             materialW,                 //its material
                             materialW->GetName());     //its name

  	G4PVPlacement* fPBoxW = new G4PVPlacement(RotMat,                          //no rotation
                G4ThreeVector(62.47*mm+(fBoxSize/2)/std::cos(17*3.141592653589793/180)+(fBoxSize/2+.5*cm)*std::cos(17*3.141592653589793/180),92.6*mm-(fBoxSize/2+.5*cm)*std::sin(17*3.141592653589793/180),0*cm),           //at (0,0,0)
                            fLBoxW,                      //its logical volume
                            materialW->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);
 }




 else if (fMaterial->GetName()=="Uranium")
 {


  G4Box*
  sBox = new G4Box("Container",.5*cm,.5*cm,fBoxSize/2);   //its dimensions

  fLBox = new G4LogicalVolume(sBox,                     //its shape
                             fMaterial,                 //its material
                             fMaterial->GetName());     //its name

  fPBox = new G4PVPlacement(RotMat,                          //no rotation
                            G4ThreeVector(62.47*mm+fBoxSize/2/std::cos(17*3.141592653589793/180),92.6*mm,0*cm),           //at (0,0,0)
                            fLBox,                      //its logical volume
                            fMaterial->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);






 	G4Isotope* isotopeW180 = new G4Isotope("W", 74, 180);
 	G4Isotope* isotopeW182 = new G4Isotope("W", 74, 182);
 	G4Isotope* isotopeW183 = new G4Isotope("W", 74, 183);
 	G4Isotope* isotopeW184 = new G4Isotope("W", 74, 184);
 	G4Isotope* isotopeW186 = new G4Isotope("W", 74, 186);

 	G4Element* elementW  = new G4Element("Tungsten", "W", ncomponents1=5);
 	elementW->AddIsotope(isotopeW180, abundance1=  0.12*perCent);
 	elementW->AddIsotope(isotopeW182, abundance1= 26.50*perCent);
 	elementW->AddIsotope(isotopeW183, abundance1= 14.31*perCent);
 	elementW->AddIsotope(isotopeW184, abundance1= 30.64*perCent);
 	elementW->AddIsotope(isotopeW186, abundance1= 28.43*perCent);


 	G4Material* materialW = new G4Material("Tungsten", 19.3*g/cm3, ncomponents1=1);
 	materialW->AddElement(elementW, massfraction1=100.*perCent);

  	G4Box*
  	sBoxW = new G4Box("Container",.5*cm,.5*cm,.5*cm);   //its dimensions

  	G4Tubs*
  	sBoxDivot = new G4Tubs("Divot",0,1.1*mm,2*mm,0*deg,360*deg);   //its dimensions

  	G4RotationMatrix* RotMatNull=new G4RotationMatrix();
  	G4ThreeVector* TransVec=new G4ThreeVector(0,0,3*mm);

  	G4SubtractionSolid* subBoxW=new G4SubtractionSolid("Containerw/Divot", sBoxW, sBoxDivot, RotMatNull, *TransVec);

  	G4LogicalVolume* fLBoxW = new G4LogicalVolume(subBoxW,                     //its shape
                             materialW,                 //its material
                             materialW->GetName());     //its name

  	G4PVPlacement* fPBoxW = new G4PVPlacement(RotMat,                          //no rotation
                G4ThreeVector(62.47*mm+(fBoxSize/2)/std::cos(17*3.141592653589793/180)+(fBoxSize/2+.5*cm)*std::cos(17*3.141592653589793/180),92.6*mm-(fBoxSize/2+.5*cm)*std::sin(17*3.141592653589793/180),0*cm),           //at (0,0,0)
                            fLBoxW,                      //its logical volume
                            materialW->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);




 	G4Isotope* isotopeW234 = new G4Isotope("U", 92, 234);
 	G4Isotope* isotopeW235 = new G4Isotope("U", 92, 235);
 	G4Isotope* isotopeW238 = new G4Isotope("U", 92, 238);

 	G4Element* elementU  = new G4Element("Uranium", "U", ncomponents1=3);
 	elementU->AddIsotope(isotopeW234, abundance1=  0.3000*perCent);
 	elementU->AddIsotope(isotopeW235, abundance1= 99.6946*perCent);
 	elementU->AddIsotope(isotopeW238, abundance1=  0.0054*perCent);


 	G4Material* materialU = new G4Material("Uranium", 19.05*g/cm3, ncomponents1=1);
 	materialU->AddElement(elementU, massfraction1=100.*perCent);

  	G4Tubs*
  	sBoxU = new G4Tubs("U",0,1.1*mm,2*mm,0*deg,360*deg);   //its dimensions

  	G4LogicalVolume* fLBoxU = new G4LogicalVolume(sBoxU,                     //its shape
                             materialU,                 //its material
                             materialU->GetName());     //its name

  	G4PVPlacement* fPBoxU = new G4PVPlacement(RotMat,                          //no rotation
                G4ThreeVector(62.47*mm+(fBoxSize/2)/std::cos(17*3.141592653589793/180)+(fBoxSize/2+.2*cm)*std::cos(17*3.141592653589793/180),92.6*mm-(fBoxSize/2+.2*cm)*std::sin(17*3.141592653589793/180),0*cm),           //at (0,0,0)
                            fLBoxU,                      //its logical volume
                            materialU->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);

 }




  G4Tubs*
  sBox1 = new G4Tubs("Container",787.4*mm,796.93*mm,781.22/2*mm,0*deg,360*deg);   //its dimensions

  G4LogicalVolume* fLBox1 = new G4LogicalVolume(sBox1,                     //its shape
                             materialSt,                 //its material
                             materialSt->GetName());     //its name

  G4Tubs*
  sBox2 = new G4Tubs("Container",0*m,796.93*mm,10*mm,0*deg,360*deg);   //its dimensions

  G4LogicalVolume* fLBox2 = new G4LogicalVolume(sBox2,                     //its shape
                             materialSt,                 //its material
                             materialSt->GetName());     //its name




  G4PVPlacement* fPBox1 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,0),           //at (0,0,0)
                            fLBox1,                      //its logical volume
                            materialSt->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);


  G4PVPlacement* fPBox2 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,(781.22/2+5)*mm),           //at (0,0,0)
                            fLBox2,                      //its logical volume
                            materialSt->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);

  G4PVPlacement* fPBox3 = new G4PVPlacement(0,                          //no rotation
                            G4ThreeVector(0,0,-(781.22/2+5)*mm),           //at (0,0,0)
                            fLBox2,                      //its logical volume
                            materialSt->GetName(),       //its name
                            logicWorld,                          //its mother  volume
                            false,                      //no boolean operation
                            0);



    G4cout << "\nThickness of Target: "
           << fBoxSize << G4endl;

  PrintParameters();

  //always return the root volume
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The Box is " << G4BestUnit(fBoxSize,"Length")
         << " of " << fMaterial->GetName()
         << "\n \n" << fMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial) {
    if(fMaterial != pttoMaterial) {
      fMaterial = pttoMaterial;
      if(fLBox) { fLBox->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
           << materialChoice << " not found" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
