#define EventTree_cxx
#include "EventTree.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <iostream>


#define mu1_pt_cut 25.0
#define mu2_pt_cut 12.0
#define mu_eta_cut 2.40
#define mu_iso     0.25

#define ele1_pt_cut 25.0
#define ele2_pt_cut 15.0
#define ele_eta_cut 2.50

#define jet_pt_cut  35.0
#define jet_eta_cut 2.4

#define ipsig_cut        1.15      
#define track_angle_cut -1.75
#define alpha_max_cut   0.75

struct lep
{
  int charge;
  int pdg;
  float energy;
  float pt;
  float eta;
  float phi;
};


void EventTree::Loop()
{
  //--------------------
  //Crete output TTree
  //---------------------
  CreateOutputTree();//initialize output tree


  //-------------------------------------
  //Normalization histograms
  //------------------------------------
  TH1F* NEvents = new TH1F("NEvents", "NEvents", 1, 0, 1);
  TH1F* NEvents_genweight = new TH1F("NEvents_genweight", "NEvents_genweight", 1, 0, 1);

  if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //std::cout << "reset" << std::endl;
      //-----------------------------------------------
      ResetTreeBranches();//reset output tree branches
      
      //-----------------------------
      //Fill Normalization histograms
      //-----------------------------
      NEvents->Fill( 0.5, 1.0 );
      NEvents_genweight->Fill(0.5, AODGenEventWeight);

      //-------------------------------------------
      //do muon selection and create good muon list
      //-------------------------------------------
      std::vector<lep> good_muons;
      for( int imu = 0; imu < nAODMu; imu++ )
	{
	  if( AOD_muPt->at(imu) < mu2_pt_cut )           continue;
	  if( fabs(AOD_muEta->at(imu)) > mu_eta_cut )    continue;
	  if( !AOD_muPassLooseID->at(imu) )              continue;
	  if( AOD_muPFdBetaIsolation->at(imu) > mu_iso ) continue;
	  lep mu_candidate;
	  mu_candidate.charge = AOD_muCharge->at(imu);
	  mu_candidate.pdg    = 13;
	  mu_candidate.energy = AOD_muEn->at(imu);
	  mu_candidate.pt     = AOD_muPt->at(imu);
	  mu_candidate.eta    = AOD_muEta->at(imu);
	  mu_candidate.phi    = AOD_muPhi->at(imu);
	  //push into good muon vector
	  good_muons.push_back(mu_candidate);
	}
      
      //---------------------------------------------------
      //do electron selection and create good electron list
      //---------------------------------------------------
      std::vector<lep> good_electrons;
      for( int iele = 0; iele < nAODEle; iele++ )
        {
	  if( AOD_elePt->at(iele) < ele2_pt_cut ) continue;
	  if( fabs(AOD_eleEta->at(iele)) > ele_eta_cut )    continue;
	  //id
	  bool pass_id = ( (AOD_eleIDbit->at(iele) >> 1) & 0x1 ) == 1;//Loose ID
	  if( !pass_id ) continue;//remove eletrons NOT passing ID
	  //check if electron is ecal in-crack
	  bool in_crack = ( fabs(AOD_eleEta->at(iele)) > 1.442 && fabs(AOD_eleEta->at(iele)) < 1.566 );
	  if( in_crack ) continue; //remove electron in crack
	  //obtain conversion veto
	  bool pass_convsersion_veto = (AOD_elePassConversionVeto->at(iele) > 0);
	  if( !pass_convsersion_veto ) continue;//remove electrons not passing conversion veto;
	  
	  //check overlap with muons
	  TLorentzVector current_electron;
	  current_electron.SetPtEtaPhiE( AOD_elePt->at(iele), AOD_eleEta->at(iele), AOD_elePhi->at(iele), AOD_eleEn->at(iele));
	  bool matched_muon = false;
	  for( auto & temp_mu : good_muons )
	    {
	      TLorentzVector current_muon;
	      current_muon.SetPtEtaPhiE(temp_mu.pt, temp_mu.eta, temp_mu.phi, temp_mu.energy);
	      if( current_electron.DeltaR(current_muon) < 0.4 ) matched_muon = true;
	      break;
	    }
	  //do not consider electrons matched to good muons
	  if( matched_muon ) continue;

	  lep ele_candidate;
          ele_candidate.charge = AOD_eleCharge->at(iele);
          ele_candidate.pdg    = 11;
          ele_candidate.energy = AOD_eleEn->at(iele);
          ele_candidate.pt     = AOD_elePt->at(iele);
          ele_candidate.eta    = AOD_eleEta->at(iele);
          ele_candidate.phi    = AOD_elePhi->at(iele);
          //push into good electron vector
          good_electrons.push_back(ele_candidate);
	  
	}


      //Z-mumu candidate
      float dilepton_sum_pt = 0.0;
      if( good_muons.size() >= 2 )
	{
	  for( int imu1 = 0; imu1 < good_muons.size(); imu1++ )
	    {
	      for( int imu2 = imu1+1; imu2 < good_muons.size(); imu2++ )
		{
		  if( good_muons.at(imu1).charge == good_muons.at(imu2).charge ) continue;//remove same charge dimuon candidates
		  float current_dilepton_sum_pt = good_muons.at(imu1).pt + good_muons.at(imu2).pt;
		  if( current_dilepton_sum_pt > dilepton_sum_pt )
		    {
		      dilepton_sum_pt = current_dilepton_sum_pt;//always select Z with scalar sum pt being the largest
		      //found better dilep candidate
		      TLorentzVector mu1;
		      mu1.SetPtEtaPhiE( good_muons.at(imu1).pt, good_muons.at(imu1).eta, good_muons.at(imu1).phi, good_muons.at(imu1).energy );
		      TLorentzVector mu2;
		      mu2.SetPtEtaPhiE( good_muons.at(imu2).pt, good_muons.at(imu2).eta, good_muons.at(imu2).phi, good_muons.at(imu2).energy );
		      TLorentzVector my_dilepton;
		      my_dilepton = mu1+mu2;
		      dilep_mass = my_dilepton.M();
		      dilep_pt   = my_dilepton.Pt();
		      dilep_eta  = my_dilepton.Eta();
		      dilep_phi  = my_dilepton.Phi();
		      //set lepton1 info
		      lep1_pdgID  = good_muons.at(imu1).pdg;
		      lep1_charge = good_muons.at(imu1).charge;
		      lep1_pt     = good_muons.at(imu1).pt;
		      lep1_eta    = good_muons.at(imu1).eta;
		      lep1_phi    = good_muons.at(imu1).phi;
		      //set lepton2 info
		      lep2_pdgID  = good_muons.at(imu2).pdg;
		      lep2_charge = good_muons.at(imu2).charge;
		      lep2_pt     = good_muons.at(imu2).pt;
		      lep2_eta    = good_muons.at(imu2).eta;
		      lep2_phi    = good_muons.at(imu2).phi;
		    }
		}
	    }
	  
	}
      else if ( good_electrons.size() >= 2 )
	{
	  dilepton_sum_pt = 0.0;
	  for( int iele1 = 0; iele1 < good_electrons.size(); iele1++ )
            {
	      for( int iele2 = iele1+1; iele2 < good_electrons.size(); iele2++ )
		{
		  if( good_electrons.at(iele1).charge == good_electrons.at(iele2).charge ) continue;//remove same charge dielectron candidates
                  float current_dilepton_sum_pt = good_electrons.at(iele1).pt + good_electrons.at(iele2).pt;
                  if( current_dilepton_sum_pt > dilepton_sum_pt )
                    {
		      dilepton_sum_pt = current_dilepton_sum_pt;//always select Z with scalar sum pt being the largest
		      //found better dilep candidate
                      TLorentzVector ele1;
		      ele1.SetPtEtaPhiE( good_electrons.at(iele1).pt, good_electrons.at(iele1).eta, good_electrons.at(iele1).phi, good_electrons.at(iele1).energy );
		      TLorentzVector ele2;
		      ele2.SetPtEtaPhiE( good_electrons.at(iele2).pt, good_electrons.at(iele2).eta, good_electrons.at(iele2).phi, good_electrons.at(iele2).energy );
		      TLorentzVector my_dilepton;
		      my_dilepton = ele1+ele2;
		      dilep_mass = my_dilepton.M();
		      dilep_pt = my_dilepton.Pt();
		      dilep_eta  = my_dilepton.Eta();
		      dilep_phi  = my_dilepton.Phi();
		      //set lepton1 info
		      lep1_pdgID  = good_electrons.at(iele1).pdg;
		      lep1_charge = good_electrons.at(iele1).charge;
		      lep1_pt     = good_electrons.at(iele1).pt;
		      lep1_eta    = good_electrons.at(iele1).eta;
		      lep1_phi    = good_electrons.at(iele1).phi;
		      //set lepton2 info
                      lep2_pdgID  = good_electrons.at(iele2).pdg;
                      lep2_charge = good_electrons.at(iele2).charge;
                      lep2_pt     = good_electrons.at(iele2).pt;
                      lep2_eta    = good_electrons.at(iele2).eta;
                      lep2_phi    = good_electrons.at(iele2).phi;
		    }
                }
            }
	}
      else if( good_muons.size() == 1 && good_electrons.size() == 1 )
	{
	  //this will be e-mu region -- TO BE DONE
	}


      //-----------------------------------------
      //do jet porting into tree_out variables
      //-----------------------------------------
      for( int ijet = 0; ijet < AODCaloJetPt->size(); ijet++ )
	{
	  if ( AODCaloJetPt->at(ijet) < jet_pt_cut ) continue;
	  if ( fabs(AODCaloJetEta->at(ijet)) > jet_eta_cut ) continue;

	  //------------------------------
	  //check overlap with muons
	  //------------------------------
          TLorentzVector current_jet;
          current_jet.SetPtEtaPhiE( AODCaloJetPt->at(ijet), AODCaloJetEta->at(ijet), AODCaloJetPhi->at(ijet), AODCaloJetEnergy->at(ijet));
          bool matched_muon = false;
          for( auto & temp_mu : good_muons )
            {
              TLorentzVector current_muon;
              current_muon.SetPtEtaPhiE(temp_mu.pt, temp_mu.eta, temp_mu.phi, temp_mu.energy);
              if( current_jet.DeltaR(current_muon) < 0.4 ) matched_muon = true;
              break;
            }
          //do not consider jet matched to good muons
          if( matched_muon ) continue;


	  //------------------------------
          //check overlap with electrons
          //------------------------------
          bool matched_electron = false;
          for( auto & temp_ele : good_electrons )
            {
              TLorentzVector current_ele;
              current_ele.SetPtEtaPhiE(temp_ele.pt, temp_ele.eta, temp_ele.phi, temp_ele.energy);
              if( current_jet.DeltaR(current_ele) < 0.4 ) matched_electron = true;
              break;
            }
          //do not consider jet matched to good muons
          if( matched_electron ) continue;

	  //--------------------------------------
	  //Save Jet Quantities in TTree variables
	  //--------------------------------------
	  jet_pt[n_jets]          = AODCaloJetPt->at(ijet);
	  jet_eta[n_jets]         = AODCaloJetEta->at(ijet);
	  jet_phi[n_jets]         = AODCaloJetPhi->at(ijet);
	  jet_ipsig[n_jets]       = AODCaloJetMedianLog10IPSig->at(ijet);
	  jet_track_angle[n_jets] = AODCaloJetMedianLog10TrackAngle->at(ijet);
	  jet_alpha_max[n_jets]   = AODCaloJetAlphaMax->at(ijet);
	  if( jet_ipsig[n_jets] > ipsig_cut && jet_track_angle[n_jets] > track_angle_cut && jet_alpha_max[n_jets] < alpha_max_cut ) ntags++;
	  n_jets++;
	}

      //----------------------------------
      //per-event requirements
      //at least 1 jet
      //zpt>0 and mll>50
      //----------------------------------
      if( n_jets == 0 )  continue;
      if( dilep_pt < 0 ) continue;
      if( dilep_mass < 50.0 ) continue;
      tree_out->Fill();
      
   }
   
   
   file_out = new TFile( this->fout_name.c_str(), "recreate");
   NEvents->Write();
   NEvents_genweight->Write();
   tree_out->Write();
   file_out->Close();
   
}


void EventTree::CreateOutputTree()
{
  tree_out = new TTree("displaced_jet", "displaced_jet");

  //define event output branches;
  tree_out->Branch("run",               &run,     "run/I");      //run number
  tree_out->Branch("lumis",   &lumis,     "lumis/I");      //lumi
  tree_out->Branch("event",             &event,     "event/l");      //event number
  //gen event-genWeight
  tree_out->Branch("gen_weight",         &AODGenEventWeight,     "gen_weight/F");      //event number

  //dilepton
  //tree_out->Branch("",      &,     "/");      //
  tree_out->Branch("dilep_mass",     &dilep_mass,    "dilep_mass/F");
  tree_out->Branch("dilep_pt",       &dilep_pt,      "dilep_pt/F");
  tree_out->Branch("dilep_eta",      &dilep_eta,     "dilep_eta/F");
  tree_out->Branch("dilep_phi",      &dilep_phi,     "dilep_phi/F");
  
  //leptons
  tree_out->Branch("lep1_pdgID",      &lep1_pdgID,    "lep1_pdgID/I");
  tree_out->Branch("lep1_charge",     &lep1_charge,   "lep1_charge/I");
  tree_out->Branch("lep1_pt",         &lep1_pt,       "lep1_pt/F");
  tree_out->Branch("lep1_eta",        &lep1_eta,      "lep1_eta/F");
  tree_out->Branch("lep1_phi",        &lep1_phi,      "lep1_phi/F");
  //
  tree_out->Branch("lep2_pdgID",      &lep2_pdgID,    "lep2_pdgID/I");
  tree_out->Branch("lep2_charge",     &lep2_charge,   "lep2_charge/I");
  tree_out->Branch("lep2_pt",         &lep2_pt,       "lep2_pt/F");
  tree_out->Branch("lep2_eta",        &lep2_eta,      "lep2_eta/F");
  tree_out->Branch("lep2_phi",        &lep2_phi,      "lep2_phi/F");

  //jets
  tree_out->Branch("n_jets",             &n_jets,             "n_jets/I");
  tree_out->Branch("jet_pt",             &jet_pt,             "jet_pt[n_jets]/F");
  tree_out->Branch("jet_eta",            &jet_eta,            "jet_eta[n_jets]/F");
  tree_out->Branch("jet_phi",            &jet_phi,            "jet_phi[n_jets]/F");
  tree_out->Branch("jet_ipsig",          &jet_ipsig,          "jet_ipsig[n_jets]/F");
  tree_out->Branch("jet_track_angle",    &jet_track_angle,    "jet_track_angle[n_jets]/F");
  tree_out->Branch("jet_alpha_max",      &jet_alpha_max,      "jet_alpha_max[n_jets]/F");

  //ntags
  tree_out->Branch("ntags",             &ntags,             "ntags/I");

};

void EventTree::ResetTreeBranches()
{
  dilep_mass   = -999.;
  dilep_pt     = -999.;
  dilep_eta    = -999.;
  dilep_phi    = -999.;
  lep1_pdgID   = -999;
  lep1_charge  = -999;
  lep1_pt      = -999.;
  lep1_eta     = -999.;
  lep1_phi     = -999.;
  lep2_pdgID   = -999;
  lep2_charge  = -999;
  lep2_pt      = -999.;
  lep2_eta     = -999.;
  lep2_phi     = -999.;
  n_jets = 0;
  ntags  = 0; 
  for(int i = 0; i < NJETS; i++)
    {
      jet_pt[NJETS]          = -999.;
      jet_eta[NJETS]         = -999.;
      jet_phi[NJETS]         = -999.;
      jet_ipsig[NJETS]       = -999.;
      jet_track_angle[NJETS] = -999.;
      jet_alpha_max[NJETS]   = -999.;
    }
};
