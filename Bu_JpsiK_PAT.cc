#include <memory>
#include "Bu_JpsiK_PAT.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerAlgorithm.h"
#include "DataFormats/PatCandidates/interface/TriggerCondition.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"

typedef math::Error<3>::type CovarianceMatrix;

Bu_JpsiK_PAT::Bu_JpsiK_PAT(const edm::ParameterSet& iConfig):
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<std::vector<pat::GenericParticle>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),
  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),

  tree_(0), 
 
  nB(0), nMu(0),

  J_px1(0), J_py1(0), J_pz1(0),
  J_px2(0), J_py2(0), J_pz2(0),
  B_k_px(0), B_k_py(0), B_k_pz(0),

  pi1Q(0),
  pi1_trackerhits(0),
  pi1_pixelhits(0),
  pi1dz(0), pi1dzE(0),
  pi1dxy(0), pi1dxyE(0),

  B_mass(0), B_px(0), B_py(0), B_pz(0), B_charge(0),
  B_k_charge1(0), B_k_px_track(0), B_k_py_track(0), B_k_pz_track(0),
  
  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),

  B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),

  muon1Q(0), muon2Q(0), muon1Soft(0), muon2Soft(0),
  muon1StationTight(0),muon2StationTight(0),
  mu1_trackerhits(0), mu2_trackerhits(0),
  mu1_pixelhits(0), mu2_pixelhits(0),
  mu1dz(0), mu1dzE(0), mu2dz(0), mu2dzE(0), 
  mu1dxy(0), mu1dxyE(0), mu2dxy(0), mu2dxyE(0),
  
  nVtx(0),
  nTrk(0), nTrk_Vtx(0),
  nTrk_VtxIP(0), nTrk_VtxIP_Q(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),

  pVtxIPX(0),   pVtxIPY(0),   pVtxIPZ(0), pVtxIPXE(0),   pVtxIPYE(0),   pVtxIPZE(0), pVtxIPCL(0),
  pVtxIPXYE(0),   pVtxIPXZE(0),   pVtxIPYZE(0),

  priRfVtxX(0), priRfVtxY(0), priRfVtxZ(0), priRfVtxXE(0), priRfVtxYE(0), priRfVtxZE(0), priRfVtxCL(0),
  priRfVtxXYE(0), priRfVtxXZE(0), priRfVtxYZE(0),
  priRfNTrkDif(0),
 
  B_Prob(0), B_J_Prob(0), 
 
  B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
  B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
  B_DecayVtxXYE(0),   B_DecayVtxXZE(0),   B_DecayVtxYZE(0),

  B_J_DecayVtxX(0),   B_J_DecayVtxY(0),   B_J_DecayVtxZ(0),
  B_J_DecayVtxXE(0),  B_J_DecayVtxYE(0),  B_J_DecayVtxZE(0),
  B_J_DecayVtxXYE(0), B_J_DecayVtxXZE(0), B_J_DecayVtxYZE(0),

  trigger(0),
  tri_PAL1DoubleMu0_open(0),
  tri_PAL1DoubleMu0(0), tri_PAL2DoubleMu0(0), tri_PAL3DoubleMu0(0),
  
  run(0), event(0),
  lumiblock(0)

{}


Bu_JpsiK_PAT::~Bu_JpsiK_PAT(){}

void Bu_JpsiK_PAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  // Get event content information

  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle<std::vector<pat::GenericParticle> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<edm::TriggerResults> triggerResults_handle;
  iEvent.getByToken(triggerResults_Label, triggerResults_handle);

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_, pruned);

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion3_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_bc_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_bc_ct = -9999.;

  if ( (isMC_ || OnlyGen_) && pruned.isValid() ) {
    int foundit = 0;
    for (size_t i=0; i<pruned->size(); i++) {
      foundit = 0;
      const reco::Candidate *dau = &(*pruned)[i];
      if ( (abs(dau->pdgId()) == 521) ) { //&& (dau->status() == 2) ) {
            foundit++;
            gen_bc_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
            gen_bc_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
            for (size_t k=0; k<dau->numberOfDaughters(); k++) {
              const reco::Candidate *gdau = dau->daughter(k);
              if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {
                foundit++;
                gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
                gen_bc_ct = GetLifetime(gen_bc_p4,gen_bc_vtx,gen_jpsi_vtx);
                int nm=0;
                for (size_t l=0; l<gdau->numberOfDaughters(); l++) {
                  const reco::Candidate *mm = gdau->daughter(l);
                  if (mm->pdgId()==13) { foundit++;
                     if (mm->status()!=1) {
                        for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                           const reco::Candidate *mu = mm->daughter(m);
                           if (mu->pdgId()==13 ) { //&& mu->status()==1) {
                              nm++;
                              gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
                              break;
                           }
                        }
                     } else {
                       gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
                       nm++;
                     }
                  }
                  if (mm->pdgId()==-13) { foundit++;
                     if (mm->status()!=1) {
                        for (size_t m=0; m<mm->numberOfDaughters(); m++) {
                           const reco::Candidate *mu = mm->daughter(m);
                           if (mu->pdgId()==-13 ) { //&& mu->status()==1) {
                              nm++;
                              gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
                              break;
                           }
                        }
                     } else {
                       gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
                       nm++;
                     }
                  }
                }
                if (nm==2) gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                else foundit-=nm;
              }
              if ((abs(gdau->pdgId())==211 || abs(gdau->pdgId())==321) && gdau->status()==1) { foundit++;
                gen_pion3_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
              }
            } // for (size_t k
      }   // if (abs(dau->pdgId())==541 )
      if (foundit>=5) break;
    } // for i
    if (foundit!=5) {
       gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
       gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
       gen_bc_vtx.SetXYZ(0.,0.,0.);
       gen_jpsi_vtx.SetXYZ(0.,0.,0.);
       gen_bc_ct = -9999.;
       std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
    }
 }

 if ( OnlyGen_ ) { 
    tree_->Fill();
    return;
 }

 nB = 0; nMu = 0;
 trigger = 0;

 if ( triggerResults_handle.isValid()) {
   const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
   unsigned int NTRIGGERS = 7;     	          
   // para el 2016
   std::string TriggersToTest[NTRIGGERS] = {
     "HLT_PAL1DoubleMu0","HLT_PAL1DoubleMu10","HLT_PAL1DoubleMuOpen",
     "HLT_PAL2DoubleMu0","HLT_PAL2DoubleMu10",
     "HLT_PAL3DoubleMu0","HLT_PAL3DoubleMu10"};
   
   for (unsigned int i = 0; i < NTRIGGERS; i++) {
     for (int version = 1; version < 9; version++) {
       std::stringstream ss;
       ss << TriggersToTest[i] << "_v" << version;
       unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
       if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
	 trigger += (1<<i);
	 break;
       }
     }
   }      
 } else std::cout << "*** NO triggerResults found " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
 
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(BSLabel_, beamSpotHandle);
  //iEvent.getByLabel(BSLabel_, beamSpotHandle);
  if ( beamSpotHandle.isValid() ) beamSpot = *beamSpotHandle; 
  else std::cout << "No beam spot available from EventSetup" << endl;
 
  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  //edm::Handle<std::vector<reco::Vertex> > primaryVertices_handle;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  // get primary vertex
  bestVtx = *(primaryVertices_handle->begin());

  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  //priVtxXE = bestVtx.xError();
  //priVtxYE = bestVtx.yError();
  //priVtxZE = bestVtx.zError();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);

  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
  nVtx = primaryVertices_handle->size();
  nTrk = thePATTrackHandle->size();
  nTrk_Vtx = bestVtx.tracksSize();


  //lumiblock = iEvent.id().luminosityBlock();
  //run = iEvent.id().run();
  //event = iEvent.id().event();


  //*****************************************
  //Let's begin by looking for J/psi+pi^+

  unsigned int nMu_tmp = thePATMuonHandle->size();
  nMu = nMu_tmp;

  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  // **************** Prelimary cuts on muons *********************
	  if(iMuon1->track()->pt()<1.5) continue;// pt>1.5??
	  if(iMuon2->track()->pt()<1.5) continue;
	  if(fabs(iMuon1->eta())>2.5 || fabs(iMuon2->eta())>2.5) continue;


	  // "softmuons". If commented is because this condition is soo strong for HI?
	  //if( !(iMuon1->isSoftMuon(bestVtx)) ) continue;
	  //if( !(iMuon2->isSoftMuon(bestVtx)) ) continue;
	  //if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  //if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
	  
	  if(!(iMuon1->isTrackerMuon()))continue;
	  if(!(iMuon2->isTrackerMuon()))continue;
	  
	  if(!(iMuon1->isGlobalMuon()))continue;
	  if(!(iMuon2->isGlobalMuon()))continue;
	  
	  if(iMuon1->track()->hitPattern().numberOfValidPixelHits()<1)continue;	  
	  if(iMuon1->track()->hitPattern().numberOfValidHits()<6)continue;

	  if(iMuon2->track()->hitPattern().numberOfValidPixelHits()<1)continue;	  
	  if(iMuon2->track()->hitPattern().numberOfValidHits()<6)continue;

	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	 // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  // *****  end DCA for the 2 muons *********************
	  
	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;
	  
	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  VirtualKinematicParticleFactory vFactory;
		  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	  muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	 	  
	  KinematicParticleVertexFitter fitter;   	  
	  RefCountedKinematicTree psiVertexFitTree;
	  try{
	  psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }
	  //psiVertexFitTree = fitter.fit(muonParticles); 	  
	  if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }
	  
	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }
	  	  
	  //if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
	  if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue; //Juts Jpsi
	  //if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.8) continue;// consider Jpsi and psi(2S)

	  double Omb_J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	  if(Omb_J_Prob_tmp<0.01)continue;

	  //Now that we have a J/psi candidate, we look for K^+ candidates

	  //for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin(); 
	  //	   iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) {		   
	  for(std::vector<pat::GenericParticle>::const_iterator iTrack1 = thePATTrackHandle->begin();
	       iTrack1 != thePATTrackHandle->end(); ++iTrack1 ) 
	    {	     	      
	      pat::GenericParticle patTrack1 = *iTrack1;
	      
	      if(iTrack1->track()->charge()==0) continue;
	      if(iTrack1->track()->pt()<0.4) continue;
	      if(fabs(iTrack1->track()->eta())>2.5) continue;

	      //if(iTrack1->track()->hitPattern().numberOfValidPixelHits()<1)continue;
	      //if(iTrack1->track()->numberOfValidHits()<5)continue;
	      //if(!(iTrack1->track()->quality(reco::TrackBase::highPurity))) continue;
	      
	      //Now let's checks if our muons do not use the same tracks as we are using now
	      if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;		    
	      
	      reco::TransientTrack pionTT((*theB).build(iTrack1->track()));
	      
	      //ParticleMass pion_mass = 0.13957018;
	      //float pion_sigma = pion_mass*1.e-6;
	      ParticleMass kaon_mass = 0.493677;
	      float kaon_sigma = kaon_mass*1.e-6;
	      
	      float chi = 0.;
	      float ndf = 0.;
	      
	      // ***************************
	      // Jpsipion invariant mass (before kinematic vertex fit)
	      // ***************************
	      TLorentzVector pion14V, Jpsi4V; 
	      pion14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),kaon_mass);
	      
	      Jpsi4V.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());
	      
	      if ( (pion14V + Jpsi4V).M()<4.5 || (pion14V + Jpsi4V).M()>6.4 ) continue;
	      
	      //Now we are ready to combine!
	      // JPsi mass constraint is applied in the final Bu fit,
	      
	      vector<RefCountedKinematicParticle> vFitMCParticles;
	      vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	      vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	      vFitMCParticles.push_back(pFactory.particle(pionTT,kaon_mass ,chi,ndf,kaon_sigma));
	      
	      MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
	      KinematicConstrainedVertexFitter kcvFitter;
	      RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
	      if (!vertexFitTree->isValid()) {
		//std::cout << "caught an exception in the B vertex fit with MC" << std::endl;
		continue;
	      }
	      vertexFitTree->movePointerToTheTop();
	      
	      RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
	      RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();
	      if (!bDecayVertexMC->vertexIsValid()){
		// cout << "B MC fit vertex is not valid" << endl;
		continue;
	      }
	      
	      if ( (bCandMC->currentState().mass() < 4.9) || (bCandMC->currentState().mass() > 6.0) ) {
		// (debug) cout << "continue from bmass > 6.0 or < 4.9 = " << bCandMC->currentState().mass() << endl;
		continue;
	      }
	      
	      double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
	      if(B_Prob_tmp<0.01)continue;
	      
	      // get children from final Bu fit		    
	      vertexFitTree->movePointerToTheFirstChild();
	      RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();
	      
	      vertexFitTree->movePointerToTheNextChild();
	      RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();
	      
	      vertexFitTree->movePointerToTheNextChild();
	      RefCountedKinematicParticle T1CandMC = vertexFitTree->currentParticle();
	      
	      KinematicParameters Mu1KP = mu1CandMC->currentState().kinematicParameters();
	      KinematicParameters Mu2KP = mu2CandMC->currentState().kinematicParameters();		   
	      KinematicParameters Pi1KP = T1CandMC->currentState().kinematicParameters();
	      
	      // ********************* loop over all the primary vertices and we choose the one with the best pointing angle **************** 
	      reco::Vertex bestVtxIP;
	      
	      Double_t pVtxIPX_temp = -10000.0;
	      Double_t pVtxIPY_temp = -10000.0;
	      Double_t pVtxIPZ_temp = -10000.0;
	      Double_t pVtxIPXE_temp = -10000.0;
	      Double_t pVtxIPYE_temp = -10000.0;
	      Double_t pVtxIPZE_temp = -10000.0;
	      Double_t pVtxIPXYE_temp = -10000.0;
	      Double_t pVtxIPXZE_temp = -10000.0;
	      Double_t pVtxIPYZE_temp = -10000.0;
	      Double_t pVtxIPCL_temp = -10000.0;	
	      Double_t lip1 = -1000000.0;
	      for(size_t i = 0; i < primaryVertices_handle->size(); ++i) {
		const Vertex &vtx = (*primaryVertices_handle)[i];
		
		Double_t dx1 = (*bDecayVertexMC).position().x() - vtx.x(); 
		Double_t dy1 = (*bDecayVertexMC).position().y() - vtx.y();
		Double_t dz1 = (*bDecayVertexMC).position().z() - vtx.z();
		float cosAlphaXYb1 = ( bCandMC->currentState().globalMomentum().x() * dx1 + bCandMC->currentState().globalMomentum().y()*dy1 + bCandMC->currentState().globalMomentum().z()*dz1  )/( sqrt(dx1*dx1+dy1*dy1+dz1*dz1)* bCandMC->currentState().globalMomentum().mag() );
		
		if(cosAlphaXYb1>lip1)
		  {
		    lip1 = cosAlphaXYb1 ;
		    pVtxIPX_temp = vtx.x();
		    pVtxIPY_temp = vtx.y();
		    pVtxIPZ_temp = vtx.z();
		    pVtxIPXE_temp = vtx.covariance(0, 0);
		    pVtxIPYE_temp = vtx.covariance(1, 1);
		    pVtxIPZE_temp = vtx.covariance(2, 2);
		    pVtxIPXYE_temp = vtx.covariance(0, 1);
		    pVtxIPXZE_temp = vtx.covariance(0, 2);
		    pVtxIPYZE_temp = vtx.covariance(1, 2);
		    pVtxIPCL_temp = (TMath::Prob(vtx.chi2(),(int)vtx.ndof()) );
		    
		    bestVtxIP = vtx;
		    
		  }
		
	      }
	      
	      // try refitting the primary without the tracks in the B reco candidate		     		     
	      const reco::Muon *rmu1 = dynamic_cast<const reco::Muon *>(iMuon1->originalObject());
	      const reco::Muon *rmu2 = dynamic_cast<const reco::Muon *>(iMuon2->originalObject());
	      reco::TrackRef patTrack1_1 = iTrack1->track();
	      
	      // first get tracks from the original primary and count the tracks with pT>0.4GeV, |eta|<2.4 (, and highPurity?)
	      vector<reco::TransientTrack> vertexTracks;
	      int nTrk_VtxIP_Q_tmp = 0;
	      
	      for ( std::vector<TrackBaseRef >::const_iterator iTrack = bestVtxIP.tracks_begin();
		    iTrack != bestVtxIP.tracks_end(); ++iTrack) {

		reco::TrackRef trackRef = iTrack->castTo<TrackRef>();

		// count the tracks with pT>0.4GeV, |eta|<2.4, and highPurity
		//if( (trackRef->quality(reco::TrackBase::highPurity)) && (trackRef->pt()>0.4) && fabs(trackRef->eta())<2.4 ){ nTrk_VtxIP_Q_tmp++;}
		if( (trackRef->pt()>0.4) && fabs(trackRef->eta())<2.4 ){ nTrk_VtxIP_Q_tmp++;}
		//nTrk_VtxIP_Q_tmp++;

		// Now. Compare primary tracks to check for matches with B cand
		
		// the 3 tracks in the Bc candidate are  patTrack1_1 rmu1 and rmu2 
		if (  !( (patTrack1_1.key()==trackRef.key()) || (rmu1->track().key()==trackRef.key()) || (rmu2->track().key()==trackRef.key()) ) ) {
		  
		  //TransientTrack tt(trackRef, &(*bFieldHandle) );
		  reco::TransientTrack tt((*theB).build(trackRef));
		  vertexTracks.push_back(tt);
		}//else { std::cout << "found track match with primary" << endl;}
	      }
	      
	      // *** if no tracks in primary or no reco track included in primary then don't do anything ***
	      reco::Vertex bestVtxRf = bestVtxIP;
	      GlobalPoint PVRfP = GlobalPoint( bestVtxIP.x(), bestVtxIP.y(), bestVtxIP.z() );
	      
	      if (  vertexTracks.size()>0 && (bestVtxIP.tracksSize()!=vertexTracks.size()) ) {
		AdaptiveVertexFitter theFitter;
		TransientVertex v = theFitter.vertex(vertexTracks,PVRfP);
		if ( v.isValid() ) {		    
		  //set bestVtxRf as new best vertex to fill variables for refitting PV
		  bestVtxRf = reco::Vertex(v);
		  
		}
	      }

	      //cout << "bestVtxIP.tracksSize, nTrk_VtxIP_Q_tmp =  " << bestVtxIP.tracksSize() << ", " << nTrk_VtxIP_Q_tmp << endl;
	      
	      // ************ fill candidate variables now
	      
	      // You can get the momentum components from the final Bc childrens
	      J_px1->push_back( Mu1KP.momentum().x() );
	      J_py1->push_back( Mu1KP.momentum().y() );
	      J_pz1->push_back( Mu1KP.momentum().z() );
	      
	      J_px2->push_back(Mu2KP.momentum().x());
	      J_py2->push_back(Mu2KP.momentum().y());
	      J_pz2->push_back(Mu2KP.momentum().z());

	      B_k_px->push_back( Pi1KP.momentum().x() );
	      B_k_py->push_back( Pi1KP.momentum().y() );
	      B_k_pz->push_back( Pi1KP.momentum().z() );
	      
	      // now we gill fill our "nominal" variables
	      pi1Q->push_back(iTrack1->track()->quality(reco::TrackBase::highPurity));
	      pi1_trackerhits->push_back(iTrack1->track()->numberOfValidHits() );
	      pi1_pixelhits->push_back(iTrack1->track()->hitPattern().numberOfValidPixelHits() );
	      pi1dz->push_back(iTrack1->track()->dz(bestVtxIP.position()) );
	      pi1dzE->push_back(iTrack1->track()->dzError() );
	      pi1dxy->push_back(iTrack1->track()->dxy(bestVtxIP.position()) );
	      pi1dxyE->push_back(iTrack1->track()->dxyError() );
	      
	      B_mass->push_back(bCandMC->currentState().mass());
	      B_px->push_back(bCandMC->currentState().globalMomentum().x());
	      B_py->push_back(bCandMC->currentState().globalMomentum().y());
	      B_pz->push_back(bCandMC->currentState().globalMomentum().z());
	      B_charge->push_back(bCandMC->currentState().particleCharge());
	      
	      B_k_px_track->push_back(iTrack1->px() );
	      B_k_py_track->push_back(iTrack1->py() );
	      B_k_pz_track->push_back(iTrack1->pz() );
	      B_k_charge1->push_back(iTrack1->charge() );
	      
	      B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
	      B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
	      B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
	      B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
	      
	      B_J_px1->push_back(iMuon1->track()->px());
	      B_J_py1->push_back(iMuon1->track()->py());
	      B_J_pz1->push_back(iMuon1->track()->pz());		    
	      B_J_charge1->push_back(iMuon1->charge());
	      
	      B_J_px2->push_back(iMuon2->track()->px());
	      B_J_py2->push_back(iMuon2->track()->py());
	      B_J_pz2->push_back(iMuon2->track()->pz());
	      B_J_charge2->push_back(iMuon2->charge());	   
	      
	      B_Prob    ->push_back(B_Prob_tmp);
	      B_J_Prob  ->push_back(Omb_J_Prob_tmp);			    
	      
	      B_DecayVtxX ->push_back((*bDecayVertexMC).position().x());    
	      B_DecayVtxY ->push_back((*bDecayVertexMC).position().y());
	      B_DecayVtxZ ->push_back((*bDecayVertexMC).position().z());
	      
	      B_DecayVtxXE ->push_back(bDecayVertexMC->error().cxx());   
	      B_DecayVtxYE ->push_back(bDecayVertexMC->error().cyy());   
	      B_DecayVtxZE ->push_back(bDecayVertexMC->error().czz());
	      B_DecayVtxXYE ->push_back(bDecayVertexMC->error().cyx());
	      B_DecayVtxXZE ->push_back(bDecayVertexMC->error().czx());
	      B_DecayVtxYZE ->push_back(bDecayVertexMC->error().czy());		  
	      
	      B_J_DecayVtxX ->push_back( psi_vFit_vertex_noMC->position().x() );
	      B_J_DecayVtxY ->push_back( psi_vFit_vertex_noMC->position().y() );
	      B_J_DecayVtxZ ->push_back( psi_vFit_vertex_noMC->position().z() );
	      
	      B_J_DecayVtxXE ->push_back( psi_vFit_vertex_noMC->error().cxx() );
	      B_J_DecayVtxYE ->push_back( psi_vFit_vertex_noMC->error().cyy() );
	      B_J_DecayVtxZE ->push_back( psi_vFit_vertex_noMC->error().czz() );
	      B_J_DecayVtxXYE ->push_back( psi_vFit_vertex_noMC->error().cyx() );
	      B_J_DecayVtxXZE ->push_back( psi_vFit_vertex_noMC->error().czx() );
	      B_J_DecayVtxYZE ->push_back( psi_vFit_vertex_noMC->error().czy() );
	      
	      pVtxIPX->push_back( pVtxIPX_temp);
	      pVtxIPY->push_back(  pVtxIPY_temp);	    
	      pVtxIPZ->push_back(  pVtxIPZ_temp);
	      pVtxIPXE->push_back( pVtxIPXE_temp);
	      pVtxIPYE->push_back( pVtxIPYE_temp);	    
	      pVtxIPZE->push_back( pVtxIPZE_temp);
	      pVtxIPXYE->push_back( pVtxIPXYE_temp);
	      pVtxIPXZE->push_back( pVtxIPXZE_temp);	    
	      pVtxIPYZE->push_back( pVtxIPYZE_temp);
	      pVtxIPCL->push_back(  pVtxIPCL_temp);
	      nTrk_VtxIP->push_back( bestVtxIP.tracksSize() );
	      nTrk_VtxIP_Q->push_back( nTrk_VtxIP_Q_tmp );

	      priRfVtxX->push_back( bestVtxRf.x() );
	      priRfVtxY->push_back( bestVtxRf.y() );
	      priRfVtxZ->push_back( bestVtxRf.z() );
	      priRfVtxXE->push_back( bestVtxRf.covariance(0, 0) );
	      priRfVtxYE->push_back( bestVtxRf.covariance(1, 1) );
	      priRfVtxZE->push_back( bestVtxRf.covariance(2, 2) );
	      priRfVtxXYE->push_back( bestVtxRf.covariance(0, 1) );
	      priRfVtxXZE->push_back( bestVtxRf.covariance(0, 2) );
	      priRfVtxYZE->push_back( bestVtxRf.covariance(1, 2) );		  
	      priRfVtxCL->push_back( ChiSquaredProbability((double)(bestVtxRf.chi2()),(double)(bestVtxRf.ndof())) );
	      priRfNTrkDif->push_back( bestVtxIP.tracksSize() - vertexTracks.size() );
	      
	      // ******** quality for Jpsi muons *********************
	      //cout << "global muon1 " << iMuon1->isGlobalMuon() <<  ",  global muon2 " << iMuon2->isGlobalMuon() <<endl;
	      //if(!(iMuon1->isGlobalMuon()))continue;
	      //if(!(iMuon2->isGlobalMuon()))continue;
	      
	      muon1Q->push_back(iMuon1->track()->quality(reco::TrackBase::highPurity));
	      muon2Q->push_back(iMuon2->track()->quality(reco::TrackBase::highPurity));
	      
	      muon1Soft->push_back(iMuon1->isSoftMuon(bestVtxIP));
	      muon2Soft->push_back(iMuon2->isSoftMuon(bestVtxIP));

	      muon1StationTight->push_back(muon::isGoodMuon(*iMuon1,muon::TMOneStationTight));
	      muon2StationTight->push_back(muon::isGoodMuon(*iMuon2,muon::TMOneStationTight));

	      mu1_trackerhits->push_back(iMuon1->track()->hitPattern().numberOfValidHits());
	      mu2_trackerhits->push_back(iMuon2->track()->hitPattern().numberOfValidHits());
	      
	      mu1_pixelhits->push_back(iMuon1->track()->hitPattern().numberOfValidPixelHits());
	      mu2_pixelhits->push_back(iMuon2->track()->hitPattern().numberOfValidPixelHits());

	      mu1dz->push_back(iMuon1->track()->dz(bestVtxIP.position()) );
	      mu1dzE->push_back(iMuon1->track()->dzError() );
	      mu1dxy->push_back(iMuon1->track()->dxy(bestVtxIP.position()) );
	      mu1dxyE->push_back(iMuon1->track()->dxyError() );

	      mu2dz->push_back(iMuon2->track()->dz(bestVtxIP.position()) );
	      mu2dzE->push_back(iMuon2->track()->dzError() );
	      mu2dxy->push_back(iMuon2->track()->dxy(bestVtxIP.position()) );
	      mu2dxyE->push_back(iMuon2->track()->dxyError() );

	      // ******** here we will check for muon-trigger-machint  (iMuon1 and iMuon2) ********
	      const pat::TriggerObjectStandAloneCollection muHLTMatches1_PAL1 = iMuon1->triggerObjectMatchesByFilter("hltL1fL1sDoubleMu0BptxANDL1Filtered0");
	      const pat::TriggerObjectStandAloneCollection muHLTMatches2_PAL1 = iMuon2->triggerObjectMatchesByFilter("hltL1fL1sDoubleMu0BptxANDL1Filtered0");

	      const pat::TriggerObjectStandAloneCollection muHLTMatches1_PAL1O = iMuon1->triggerObjectMatchesByFilter("hltL1fL1sDoubleMuOpenBptxANDL1Filtered0");
	      const pat::TriggerObjectStandAloneCollection muHLTMatches2_PAL1O = iMuon2->triggerObjectMatchesByFilter("hltL1fL1sDoubleMuOpenBptxANDL1Filtered0");
	      
	      const pat::TriggerObjectStandAloneCollection muHLTMatches1_PAL2 = iMuon1->triggerObjectMatchesByFilter("hltL2fL1sDoubleMuOpenBptxANDL1f0L2Filtered0");
	      const pat::TriggerObjectStandAloneCollection muHLTMatches2_PAL2 = iMuon2->triggerObjectMatchesByFilter("hltL2fL1sDoubleMuOpenBptxANDL1f0L2Filtered0");
	      
	      const pat::TriggerObjectStandAloneCollection muHLTMatches1_PAL3 = iMuon1->triggerObjectMatchesByFilter("hltL3fL1sDoubleMuOpenBptxANDL1f0L2f0L3Filtered0");
	      const pat::TriggerObjectStandAloneCollection muHLTMatches2_PAL3 = iMuon2->triggerObjectMatchesByFilter("hltL3fL1sDoubleMuOpenBptxANDL1f0L2f0L3Filtered0");
	      
	      int tri_PAL1DoubleMu0_tmp = 0, tri_PAL1DoubleMu0_open_tmp = 0,tri_PAL2DoubleMu0_tmp = 0,  tri_PAL3DoubleMu0_tmp = 0;
	      
	      if (muHLTMatches1_PAL1.size() > 0 && muHLTMatches2_PAL1.size() > 0) tri_PAL1DoubleMu0_tmp = 1;
	      if (muHLTMatches1_PAL1O.size() > 0 && muHLTMatches2_PAL1O.size() > 0) tri_PAL1DoubleMu0_open_tmp = 1;	      
	      if (muHLTMatches1_PAL2.size() > 0 && muHLTMatches2_PAL2.size() > 0) tri_PAL2DoubleMu0_tmp = 1;
	      if (muHLTMatches1_PAL3.size() > 0 && muHLTMatches2_PAL3.size() > 0) tri_PAL3DoubleMu0_tmp = 1;
	      //cout << "TrkPal1Open, TrkPal1, TrkPal2, TrkPal3 = " << tri_PAL1DoubleMu0_open_tmp << ", "<< tri_PAL1DoubleMu0_tmp << ", "<< tri_PAL2DoubleMu0_tmp << ", "<< tri_PAL3DoubleMu0_tmp << endl;

	      tri_PAL1DoubleMu0_open->push_back(tri_PAL1DoubleMu0_open_tmp);	       
	      tri_PAL1DoubleMu0->push_back(tri_PAL1DoubleMu0_tmp);
	      tri_PAL2DoubleMu0->push_back(tri_PAL2DoubleMu0_tmp);
	      tri_PAL3DoubleMu0->push_back(tri_PAL3DoubleMu0_tmp);

	      
	      nB++;	       
	      muonParticles.clear();
	      vFitMCParticles.clear();
	    }
	}
    }
  
   
   if (nB > 0 ) tree_->Fill();

   nB = 0; nMu = 0;
   trigger = 0;
   tri_PAL1DoubleMu0_open->clear();
   tri_PAL1DoubleMu0->clear(); tri_PAL2DoubleMu0->clear(); tri_PAL3DoubleMu0->clear();

   
   J_px1->clear();  J_py1->clear();  J_pz1->clear(); 
   J_px2->clear();  J_py2->clear();  J_pz2->clear();
   B_k_px->clear(); B_k_py->clear(); B_k_pz->clear();

   pi1Q->clear();
   pi1_trackerhits->clear(); 
   pi1_pixelhits->clear();
   pi1dz->clear(); pi1dzE->clear(); 
   pi1dxy->clear(); pi1dxyE->clear();

   B_charge->clear();

   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear(); 
   B_k_charge1->clear(); B_k_px_track->clear(); B_k_py_track->clear(); B_k_pz_track->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();

   B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(); B_J_charge1->clear();
   B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(); B_J_charge2->clear();

   muon1Q->clear(); muon2Q->clear(); muon1Soft->clear(); muon2Soft->clear();
   muon1StationTight->clear(); muon2StationTight->clear(); 
   mu1_trackerhits->clear(); mu2_trackerhits->clear();
   mu1_pixelhits->clear(); mu2_pixelhits->clear();
   mu1dz->clear(); mu1dzE->clear();  mu2dz->clear(); mu2dzE->clear(); 
   mu1dxy->clear(); mu1dxyE->clear(); mu2dxy->clear(); mu2dxyE->clear(); 
   

   B_Prob->clear(); B_J_Prob->clear();

   B_DecayVtxX->clear();     B_DecayVtxY->clear();     B_DecayVtxZ->clear();
   B_DecayVtxXE->clear();    B_DecayVtxYE->clear();    B_DecayVtxZE->clear();
   B_DecayVtxXYE->clear();   B_DecayVtxXZE->clear();   B_DecayVtxYZE->clear();

   B_J_DecayVtxX->clear();   B_J_DecayVtxY->clear();   B_J_DecayVtxZ->clear();
   B_J_DecayVtxXE->clear();  B_J_DecayVtxYE->clear();  B_J_DecayVtxZE->clear();
   B_J_DecayVtxXYE->clear(); B_J_DecayVtxXZE->clear(); B_J_DecayVtxYZE->clear();

   nVtx = 0;
   nTrk = 0; nTrk_Vtx = 0;
   nTrk_VtxIP->clear(); nTrk_VtxIP_Q->clear();
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0;

   pVtxIPX->clear();  pVtxIPY->clear();  pVtxIPZ->clear();
   pVtxIPXE->clear();  pVtxIPYE->clear();  pVtxIPZE->clear();  pVtxIPCL->clear();
   pVtxIPXYE->clear();  pVtxIPXZE->clear();  pVtxIPYZE->clear();

   priRfVtxX->clear(); priRfVtxY->clear(); priRfVtxZ->clear(); priRfVtxXE->clear(); priRfVtxYE->clear(); 
   priRfVtxZE->clear(); priRfVtxXYE->clear(); priRfVtxXZE->clear(); priRfVtxYZE->clear(); priRfVtxCL->clear(); 
   priRfNTrkDif->clear(); 
  
}

bool Bu_JpsiK_PAT::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

// ------------ method to check the presence of pixel hits  ------------
bool Bu_JpsiK_PAT::hasFirstLayerPixelHits(const reco::TransientTrack& track)
{
  using namespace reco;
  const HitPattern& p = track.hitPattern();      
  for (int i=0; i<p.numberOfHits(HitPattern::TRACK_HITS); i++) {
    uint32_t pattern = p.getHitPattern(HitPattern::TRACK_HITS, i);   
    if (p.pixelBarrelHitFilter(pattern) || p.pixelEndcapHitFilter(pattern) ) {
      if (p.getLayer(pattern) == 1) {
    if (p.validHitFilter(pattern)) {
      return true;
    }
      }
    }
  }
  return false;
} 

// ------------ method called once each job just before starting event loop  ------------

void
Bu_JpsiK_PAT::beginJob()
{
  
  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;
  
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bc ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");

  tree_->Branch("J_px1", &J_px1);
  tree_->Branch("J_py1", &J_py1);
  tree_->Branch("J_pz1", &J_pz1);

  tree_->Branch("J_px2", &J_px2);
  tree_->Branch("J_py2", &J_py2);
  tree_->Branch("J_pz2", &J_pz2);

  tree_->Branch("B_k_px", &B_k_px);
  tree_->Branch("B_k_py", &B_k_py);
  tree_->Branch("B_k_pz", &B_k_pz);

  tree_->Branch("pi1Q", &pi1Q);
  tree_->Branch("pi1_trackerhits", &pi1_trackerhits);
  tree_->Branch("pi1_pixelhits", &pi1_pixelhits);
  tree_->Branch("pi1dz", &pi1dz);
  tree_->Branch("pi1dzE", &pi1dzE);
  tree_->Branch("pi1dxy", &pi1dxy);
  tree_->Branch("pi1dxyE", &pi1dxyE);

  tree_->Branch("B_charge", &B_charge);
  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);

  tree_->Branch("B_k_charge1", &B_k_charge1);
  tree_->Branch("B_k_px_track", &B_k_px_track);
  tree_->Branch("B_k_py_track", &B_k_py_track);
  tree_->Branch("B_k_pz_track", &B_k_pz_track);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);

  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("muon1Q", &muon1Q);
  tree_->Branch("muon2Q", &muon2Q);
  tree_->Branch("muon1Soft", &muon1Soft);
  tree_->Branch("muon2Soft", &muon2Soft);
  tree_->Branch("muon1StationTight", &muon1StationTight);
  tree_->Branch("muon2StationTight", &muon2StationTight);
  tree_->Branch("mu1_trackerhits", &mu1_trackerhits);
  tree_->Branch("mu2_trackerhits", &mu2_trackerhits);
  tree_->Branch("mu1_pixelhits", &mu1_pixelhits);
  tree_->Branch("mu2_pixelhits", &mu2_pixelhits);

  tree_->Branch("mu1dz", &mu1dz);
  tree_->Branch("mu1dzE", &mu1dzE);
  tree_->Branch("mu1dxy", &mu1dxy);
  tree_->Branch("mu1dxyE", &mu1dxyE);

  tree_->Branch("mu2dz", &mu2dz);
  tree_->Branch("mu2dzE", &mu2dzE);
  tree_->Branch("mu2dxy", &mu2dxy);
  tree_->Branch("mu2dxyE", &mu2dxyE);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
       
  tree_->Branch("B_DecayVtxX",     &B_DecayVtxX);
  tree_->Branch("B_DecayVtxY",     &B_DecayVtxY);
  tree_->Branch("B_DecayVtxZ",     &B_DecayVtxZ);
  tree_->Branch("B_DecayVtxXE",    &B_DecayVtxXE);
  tree_->Branch("B_DecayVtxYE",    &B_DecayVtxYE);
  tree_->Branch("B_DecayVtxZE",    &B_DecayVtxZE);
  tree_->Branch("B_DecayVtxXYE",    &B_DecayVtxXYE);
  tree_->Branch("B_DecayVtxXZE",    &B_DecayVtxXZE);
  tree_->Branch("B_DecayVtxYZE",    &B_DecayVtxYZE);
 
  tree_->Branch("B_J_DecayVtxX",   &B_J_DecayVtxX);
  tree_->Branch("B_J_DecayVtxY",   &B_J_DecayVtxY);
  tree_->Branch("B_J_DecayVtxZ",   &B_J_DecayVtxZ);
  tree_->Branch("B_J_DecayVtxXE",  &B_J_DecayVtxXE);
  tree_->Branch("B_J_DecayVtxYE",  &B_J_DecayVtxYE);
  tree_->Branch("B_J_DecayVtxZE",  &B_J_DecayVtxZE);
  tree_->Branch("B_J_DecayVtxXYE",  &B_J_DecayVtxXYE);
  tree_->Branch("B_J_DecayVtxXZE",  &B_J_DecayVtxXZE);
  tree_->Branch("B_J_DecayVtxYZE",  &B_J_DecayVtxYZE);

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");

  tree_->Branch("pVtxIPX",     &pVtxIPX);
  tree_->Branch("pVtxIPY",     &pVtxIPY);
  tree_->Branch("pVtxIPZ",     &pVtxIPZ);
  tree_->Branch("pVtxIPXE",     &pVtxIPXE);
  tree_->Branch("pVtxIPYE",     &pVtxIPYE);
  tree_->Branch("pVtxIPZE",     &pVtxIPZE);
  tree_->Branch("pVtxIPXYE",     &pVtxIPXYE);
  tree_->Branch("pVtxIPXZE",     &pVtxIPXZE);
  tree_->Branch("pVtxIPYZE",     &pVtxIPYZE);
  tree_->Branch("pVtxIPCL",     &pVtxIPCL);

  tree_->Branch("priRfVtxX",&priRfVtxX);
  tree_->Branch("priRfVtxY",&priRfVtxY);
  tree_->Branch("priRfVtxZ",&priRfVtxZ);
  tree_->Branch("priRfVtxXE",&priRfVtxXE);
  tree_->Branch("priRfVtxYE",&priRfVtxYE);
  tree_->Branch("priRfVtxZE",&priRfVtxZE);
  tree_->Branch("priRfVtxXYE",&priRfVtxXYE);
  tree_->Branch("priRfVtxXZE",&priRfVtxXZE);
  tree_->Branch("priRfVtxYZE",&priRfVtxYZE);
  tree_->Branch("priRfVtxCL",&priRfVtxCL);
  tree_->Branch("priRfNTrkDif",&priRfNTrkDif);
 
  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("nTrk",       &nTrk);
  tree_->Branch("nTrk_Vtx",   &nTrk_Vtx);
  tree_->Branch("nTrk_VtxIP",   &nTrk_VtxIP);
  tree_->Branch("nTrk_VtxIP_Q",   &nTrk_VtxIP_Q);  
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");

  tree_->Branch("trigger",  &trigger, "trigger/i");
  tree_->Branch("tri_PAL1DoubleMu0_open", &tri_PAL1DoubleMu0_open);
  tree_->Branch("tri_PAL1DoubleMu0", &tri_PAL1DoubleMu0);
  tree_->Branch("tri_PAL2DoubleMu0", &tri_PAL2DoubleMu0);
  tree_->Branch("tri_PAL3DoubleMu0", &tri_PAL3DoubleMu0);

  // *************************

// gen
  if (isMC_) {
     tree_->Branch("gen_bc_p4",     "TLorentzVector",  &gen_bc_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_pion3_p4",  "TLorentzVector",  &gen_pion3_p4);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_bc_vtx",    "TVector3",        &gen_bc_vtx);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_bc_ct",     &gen_bc_ct,        "gen_bc_ct/F");
  }

}

double Bu_JpsiK_PAT::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}

// ------------ method called once each job just after ending the event loop  ------------
void Bu_JpsiK_PAT::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Bu_JpsiK_PAT);
