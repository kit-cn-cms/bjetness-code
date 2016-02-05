#include "Other/BJetness/interface/BJetness.hpp"
/////
//   Methods for the BJetness variables
/////
void BJetness::get_bjetness_vars(
                                           vector<pat::Jet> evtjets, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h,
                                           double& bjetnessFV_num_loosenoipnoiso_leps, double& bjetnessFV_numjettrksnopv, double& bjetnessFV_pvTrkOVcollTrk, double& bjetnessFV_avip3d_val, double& bjetnessFV_avip3d_sig, double& bjetnessFV_avsip3d_sig, double& bjetnessFV_avip1d_sig
                                          ){
  //Get BJetness trk info
  vector<Track> jetschtrks; jetschtrks.clear(); 
  double num_pvtrks              = 0;
  double num_npvtrks             = 0;
  double num_loosenoipnoiso_eles = 0;  
  double num_loose_mus           = 0;           
  vector<tuple<double, double, double> > jetsdir; jetsdir.clear(); 
  get_bjetness_trkinfos(evtjets, vtx, jetschtrks, num_pvtrks, num_npvtrks, electron_pat, muon_h, num_loosenoipnoiso_eles, num_loose_mus, jetsdir);
  bjetnessFV_num_loosenoipnoiso_leps = num_loosenoipnoiso_eles+num_loose_mus;
  bjetnessFV_numjettrksnopv          = num_npvtrks;
  if(jetschtrks.size()!=0){
    bjetnessFV_pvTrkOVcollTrk        = num_pvtrks/double(jetschtrks.size()); 
    //Get BJetness Impact Parameters
    double ip_valtemp = 0;
    //3D
    double jetchtrks_avip3d_val  = 0;
    double jetchtrks_avip3d_sig  = 0;
    double jetchtrks_avsip3d_sig = 0;
    get_avip3d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip3d_val,jetchtrks_avip3d_sig,jetchtrks_avsip3d_sig);
    ip_valtemp = jetchtrks_avip3d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_val = ip_valtemp;
    else                       bjetnessFV_avip3d_val = -996;
    ip_valtemp = jetchtrks_avip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_sig = ip_valtemp;
    else                       bjetnessFV_avip3d_sig = -996; 
    ip_valtemp = jetchtrks_avsip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip3d_sig = ip_valtemp;
    else                       bjetnessFV_avsip3d_sig = -996;
    //1D
    double jetchtrks_avip1d_sig  = 0;
    get_avip1d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip1d_sig);
    ip_valtemp = jetchtrks_avip1d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip1d_sig = ip_valtemp;
    else                       bjetnessFV_avip1d_sig = -996;    
  }else{
    bjetnessFV_pvTrkOVcollTrk        = -998;
    bjetnessFV_avip3d_val            = -998;
    bjetnessFV_avip3d_sig            = -998;
    bjetnessFV_avsip3d_sig           = -998;
    bjetnessFV_avip1d_sig            = -998;
  }
}
//Get the BJetness trk info 
void BJetness::get_bjetness_trkinfos(vector<pat::Jet> evtjets, const reco::Vertex& vtx, vector<Track>& jetchtrks, double& bjetness_num_pvtrks, double& bjetness_num_npvtrks, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_loosenoipnoiso_eles, double& bjetness_num_loose_mus, vector<tuple<double, double, double> >& jetsdir){
  //Loop over evt jet
  for(uint j=0; j<evtjets.size(); j++){
    pat::Jet jet = evtjets[j];
    //Access jet daughters
    vector<CandidatePtr> jdaus(jet.daughterPtrVector());
    sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
    for(uint jd=0; jd<jdaus.size(); jd++){
      const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
      //dR requirement
      if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
      Track trk = Track(jcand.pseudoTrack());
      bool isgoodtrk = is_goodtrk(trk,vtx);
      //Minimal conditions for a BJetness jet constituent 
      if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
        jetchtrks.push_back(trk);
        if(jcand.fromPV()==3) bjetness_num_pvtrks++;
        if(jcand.fromPV()==2) bjetness_num_npvtrks++;
        jetsdir.push_back(make_tuple(jet.px(),jet.py(),jet.pz()));
        if(fabs(jcand.pdgId())==13 && is_loosePOG_jetmuon(jcand,muon_h)) bjetness_num_loose_mus++;
        if(fabs(jcand.pdgId())==11 && is_loosePOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_loosenoipnoiso_eles++;
      }//Ch trks 
    }//Loop on jet daus 
  }//Loop on evt jet
}
//Check that the track is a good track
bool BJetness::is_goodtrk(Track trk,const reco::Vertex& vtx){
 bool isgoodtrk = false;
 if(trk.pt()>1 &&
   trk.hitPattern().numberOfValidHits()>=8 &&
   trk.hitPattern().numberOfValidPixelHits()>=2 &&
   trk.normalizedChi2()<5 &&
   std::abs(trk.dxy(vtx.position()))<0.2 &&
   std::abs(trk.dz(vtx.position()))<17
   ) isgoodtrk = true;
 return isgoodtrk;
}
//Look for loose muon (definition to look for candidates among jet daughters)
bool BJetness::is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h){
  bool ismu = false;
  for(const pat::Muon &mu : *muon_h){
    if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
     if(mu.isLooseMuon()) ismu = true;
     if(ismu) break;
    }
  }  
  return ismu;
}
//Look for loose electron ((definition to look for candidates among jet daughters)
bool BJetness::is_loosePOGNoIPNoIso_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx){
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele; 
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05){
      double ooEmooP = 999;
      if(lele.ecalEnergy()==0)                   ooEmooP = 1e30;
      else if(!std::isfinite(lele.ecalEnergy())) ooEmooP = 1e30;
      else                                       ooEmooP = fabs(1.0/lele.ecalEnergy() - lele.eSuperClusterOverP()/lele.ecalEnergy());
      if(!(fabs(lele.superCluster()->position().eta()) > 1.4442 && fabs(lele.superCluster()->position().eta()) < 1.5660)){
        if(fabs(lele.superCluster()->position().eta())<1.4442
          && (lele.full5x5_sigmaIetaIeta()<0.0103)
          && (lele.hcalOverEcal()<0.104)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.0105)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.115)
          && ooEmooP<0.102
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=2
          && lele.passConversionVeto()   
          ){
            isele = true;
        }
        if(fabs(lele.superCluster()->position().eta())>1.5660
          && (lele.full5x5_sigmaIetaIeta()<0.0301)
          && (lele.hcalOverEcal()<0.0897)
          && (fabs(lele.deltaEtaSuperClusterTrackAtVtx())<0.00814)
          && (fabs(lele.deltaPhiSuperClusterTrackAtVtx())<0.182)
          && ooEmooP<0.126
          && lele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<=1
          && lele.passConversionVeto()   
          ){
            isele = true;
        }
      } 
      if(isele) break;
    } 
  } 
  return isele; 
}
//Methods related to IP
void BJetness::get_avip3d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip3d_val  += valtemp;
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip3d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip3d_sig += valtemp;
  }
}

void BJetness::get_avip1d(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_sig){
  double valtemp = 0;
  vector<TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  SignedTransverseImpactParameter stip;
  for(uint t=0; t<ttrks.size(); t++){
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance());
    if(valtemp==valtemp) jetchtrks_avip1d_sig  += valtemp;
  }
}
vector<TransientTrack> BJetness::get_ttrks(vector<Track> trks, const TransientTrackBuilder& ttrkbuilder){
 vector<TransientTrack> ttrks;
 for(uint tr=0; tr<trks.size(); tr++){
  TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
  ttrks.push_back(ttrk);
 }
 return ttrks;
}
