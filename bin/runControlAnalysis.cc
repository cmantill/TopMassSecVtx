#include "UserCode/TopMassSecVtx/interface/MacroUtils.h"
#include "UserCode/TopMassSecVtx/interface/SmartSelectionMonitor.h"
#include "UserCode/TopMassSecVtx/interface/TopSelectionTools.h"
#include "UserCode/TopMassSecVtx/interface/DataEventSummaryHandler.h"
#include "UserCode/TopMassSecVtx/interface/LxyAnalysis.h"
#include "UserCode/TopMassSecVtx/interface/MuScleFitCorrector.h"
#include "UserCode/TopMassSecVtx/interface/LeptonEfficiencySF.h"
#include "UserCode/TopMassSecVtx/interface/BfragWeighter.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "TopQuarkAnalysis/TopTools/interface/MEzCalculator.h"

#include "TVectorD.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"

#include <fstream>
#include <iostream>

using namespace std;

FactorizedJetCorrector *fJesCor=0;
std::vector<JetCorrectionUncertainty *> fTotalJESUnc;
MuScleFitCorrector *fMuCor=0;
LeptonEfficiencySF fLepEff;
BfragWeighter *fBfragWgt=0;

//aggregates all the objects needed for the analysis
enum ControlBoxType { DIJETBOX=1, GAMMABOX=22, WBOX=24, ZBOX=23};
class ControlBox{
public:
  ControlBox() : cat(0), probeFlav(0), probeXb(0), probeFlavStr(""), selEffSF(1.0), probe(0), tag(0) { }

  void assignBox(Int_t evcat,data::PhysicsObject_t *itag, data::PhysicsObject_t *iprobe)
  {
    cat=evcat;
    tag=itag;
    probe=iprobe;
  }

  void assignProbeFlavour(data::PhysicsObjectCollection_t &gen)
  {
    const data::PhysicsObject_t &genJet=probe->getObject("genJet");
    if(genJet.pt()>0) {
      int jetFlav=genJet.info.find("id")->second;
      if(abs(jetFlav)==5)      { probeFlav=1; probeFlavStr="b";         }
      else if(abs(jetFlav)==4) { probeFlav=3; probeFlavStr="c";         }
      else                     { probeFlav=2; probeFlavStr="udsg";      }
    }
    else                       { probeFlav=4; probeFlavStr="unmatched"; }

    const data::PhysicsObject_t &genParton=probe->getObject("gen");
    int genId=genParton.info.find("id")->second;
    if(abs(genId)==5 && genParton.pt()>0){
      for(size_t imc=0; imc<gen.size(); imc++)
	{
	  int id=gen[imc].get("id");
	  if(abs(id)<500) continue;
	  if(deltaR(gen[imc],genParton)>0.5) continue;
	  probeXb=gen[imc].pt()/genParton.pt();
	  break;
	}
    }
  }

  ~ControlBox() { }
  Int_t cat;
  Int_t probeFlav;
  Float_t probeXb;
  TString probeFlavStr;
  Float_t selEffSF;
  data::PhysicsObject_t *probe, *tag;
};

ControlBox assignBox(int reqControlType,
		     data::PhysicsObjectCollection_t &leptons,
		     data::PhysicsObjectCollection_t &photons,
		     data::PhysicsObjectCollection_t &jets,
		     std::vector<LorentzVector> &met,
		     data::PhysicsObjectCollection_t &pf)
{
  ControlBox box;

  std::vector<int> dilCands, ljCands, vetoCands;
  for(size_t i=0; i<leptons.size(); i++){
    if(leptons[i].getFlag("passLL"))     dilCands.push_back(i);
    if(leptons[i].getFlag("passLJ"))     ljCands.push_back(i);
    else {
      if(leptons[i].getFlag("passLJveto")) vetoCands.push_back(i);
    }
  }
  if(reqControlType==DIJETBOX)
    {
      //at least two jets and no leptons
      if(jets.size()<2) return box;
      if(dilCands.size() || ljCands.size() || vetoCands.size()) return box;

      //require at least one soft muon
      std::vector<int> jetsWithSoftMuons;
      for(size_t ijet=0; ijet<jets.size(); ijet++){
	size_t pfstart=jets[ijet].get("pfstart");
	size_t pfend=jets[ijet].get("pfend");
	if(pf.size()<pfstart || pf.size()<pfend-1) continue;
	for(size_t ipfn=pfstart; ipfn<=pfend; ipfn++)
	  {
	    int id=pf[ipfn].get("id");
	    if(abs(id)!=13) continue;
	    if(pf[ipfn].pt()<5) continue;
	    jetsWithSoftMuons.push_back(ijet);
	    break;
	  }
      }
      if(jetsWithSoftMuons.size()==0) return box;

      //select one of the jets to be the probe
      size_t randomTagIdx( (size_t) gRandom->Uniform(0,jetsWithSoftMuons.size()) );
      Int_t randomProbeIdx(-1);
      for(size_t ijet=0; ijet<jets.size(); ijet++)
	{
	  if(ijet==randomTagIdx) continue;
	  float dR=deltaR(jets[ijet],jets[randomTagIdx]);
	  if(dR<2.09) continue;
	  if(randomProbeIdx<0) randomProbeIdx=ijet;
	  else if( dR>deltaR(jets[randomProbeIdx],jets[randomTagIdx])
		   && jets[ijet].pt()>jets[randomProbeIdx].pt() ) randomProbeIdx=ijet;
	}
      if(randomProbeIdx<0) return box;

      //require a tag on the tag jet
      bool tagHasCSVL(jets[randomTagIdx].getVal("csv")>0.405);
      if(!tagHasCSVL) return box;

      //all done here
      box.assignBox(reqControlType, &(jets[randomTagIdx]), &(jets[randomProbeIdx]));
    }
  else if(reqControlType==GAMMABOX)
    {
      //gamma+1 jet without leptons
      if(photons.size()!=1 || jets.size()==0) return box;
      if(dilCands.size() || ljCands.size() || vetoCands.size()) return box;

      //back to back configuration to the leading jet
      float dR=deltaR(photons[0],jets[0]);
      if(fabs(dR)<2.09) return box;

      //all done here
      box.assignBox(reqControlType, &(photons[0]), &(jets[0]));
    }
  else if(reqControlType==WBOX)
    {
      //1 lepton + 1 jet
      if(ljCands.size()!=1 || vetoCands.size()) return box;
      if(jets.size()==0) return box;

      //require significant MET
      float metsig=met[0].pt()/sqrt(jets[0].pt());
      if(metsig<3.5) return box;

      //require minimun transverse mass
      float mt(utils::cmssw::getMT<LorentzVector>( leptons[ ljCands[0] ], met[0]));
      if(mt<50) return box;

      //compute the full kinematics
      TLorentzVector metP4(met[0].px(),met[0].py(),0,met[0].pt());
      TLorentzVector leptonP4(leptons[ ljCands[0] ].px(),
			      leptons[ ljCands[0] ].py(),
			      leptons[ ljCands[0] ].pz(),
			      leptons[ ljCands[0] ].energy());
      MEzCalculator mEz;
      mEz.SetMET(metP4);
      mEz.SetLepton(leptonP4);
      double mEz_val( mEz.Calculate() );
      double vEnergy=sqrt(pow(met[0].pt(),2)+pow(mEz_val,2));
      LorentzVector v(met[0].px(),met[0].py(),mEz_val,vEnergy);
      LorentzVector lv=leptons[ ljCands[0] ]+ v;

      //back-to-back configuration
      float dphijj=deltaPhi(lv.phi(),jets[0].phi());
      if(fabs(dphijj)<2.09) return box;

      //balancing variable : better not MET smears all
      //float balance=jets[0].pt()/lv.pt();
      //if(balance<0.5 || balance>1.5) return box;

      data::PhysicsObject_t *w=new data::PhysicsObject_t(lv.px(),lv.py(),lv.pz(),lv.energy());
      w->set("id",24);
      box.assignBox(reqControlType, w, &(jets[0]));
      box.selEffSF=fLepEff.getSingleLeptonEfficiencySF(leptons[ ljCands[0] ].eta(),abs(leptons[ ljCands[0] ].get("id"))).first;
    }
  else if(reqControlType==ZBOX)
    {
      //2 leptons + 1 jet
      if(dilCands.size()<2) return box;
      if(jets.size()==0) return box;

      //require z window
      LorentzVector ll=leptons[ dilCands[0] ]+leptons[ dilCands[1] ];
      if(fabs(ll.mass()-91)>15) return box;

      //back to back configuration with leading pT jet

      float dR=deltaR(ll,jets[0]);
      if(dR<2.09) return box;

      data::PhysicsObject_t *z=new data::PhysicsObject_t(ll.px(),ll.py(),ll.pz(),ll.energy());
      z->set("id",23);
      box.assignBox(reqControlType, z, &(jets[0]));

      float eta1(fabs(leptons[ dilCands[0] ].eta())), eta2(fabs(leptons[ dilCands[1] ].eta()));
      float id1id2(abs(leptons[ dilCands[0] ].get("id"))*abs(leptons[ dilCands[1] ].get("id")));
      box.selEffSF=fLepEff.getDileptonEfficiencySF(eta1,eta2,id1id2).first;
    }

  //all done
  return box;
}

//
int main(int argc, char* argv[])
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
    return 0;
  }

  //
  // configure
  //
  const edm::ParameterSet &runProcess       = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  std::vector<std::string> urls=runProcess.getParameter<std::vector<std::string> >("input");
  TString url          = TString(urls[0]);
  TString baseDir      = runProcess.getParameter<std::string>("dirName");
  TString jecDir       = runProcess.getParameter<std::string>("jecDir");
  bool isMC            = runProcess.getParameter<bool>("isMC");
  int reqControlType   = runProcess.getParameter<int>("mctruthmode");
  double xsec          = runProcess.getParameter<double>("xsec");
  bool isV0JetsMC(isMC && (url.Contains("DYJetsToLL_50toInf") || url.Contains("TeV_WJets")));
  TString out          = runProcess.getParameter<std::string>("outdir");
  bool saveSummaryTree = runProcess.getParameter<bool>("saveSummaryTree");
  std::vector<string>  weightsFile = runProcess.getParameter<std::vector<string> >("weightsFile");

  //check running mode
  bool runDijetControl(reqControlType==DIJETBOX);
  bool runGammaControl(reqControlType==GAMMABOX);
  bool runWControl(reqControlType==WBOX);
  bool runZControl(reqControlType==ZBOX);

  //if given use weights to correct tag pt in MC
  //TGraphErrors *tagPtWeightGr=0;
  float maxPtForTagWgt(450);
  /*
  if(weightsFile.size() && isMC)
    {
      TFile *inWeightsF=TFile::Open((weightsFile[0]+"/ControlTagPtWeights.root").c_str());
      TString key("qcd");
      if(runWControl)     key="w";
      if(runZControl)     key="z";
      if(runGammaControl) key="photon";
      tagPtWeightGr=(TGraphErrors *)inWeightsF->Get(key);
      inWeightsF->Close();
      if(tagPtWeightGr) std::cout << "Tag pT will be reweighted" << std::endl;
    }
  */


  //jet energy scale uncertainties
  gSystem->ExpandPathName(jecDir);
  fJesCor = utils::cmssw::getJetCorrector(jecDir,isMC);
  TString srcnames[]={
    "Total",
    "AbsoluteMPFBias", //in-situ correlation group
    "RelativeFSR",     //JEC inter-calibration
    "PileUpDataMC", "PileUpPtBB", "PileUpPtEC", "PileUpPtHF",      //Pileup
    "AbsoluteStat", "AbsoluteScale","HighPtExtra", "SinglePionECAL", "SinglePionHCAL", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB", "RelativePtEC1",  "RelativePtEC2", "RelativePtHF", "RelativeStatEC2", "RelativeStatHF", //JEC uncorrelated
    "FlavorPureGluon", "FlavorPureQuark", "FlavorPureCharm", "FlavorPureBottom" //flavor JES
  };
  const int nsrc = sizeof(srcnames)/sizeof(TString);
  for (int isrc = 0; isrc < nsrc; isrc++) {
    JetCorrectorParameters *p = new JetCorrectorParameters((jecDir+"/DATA_UncertaintySources_AK5PFchs.txt").Data(), srcnames[isrc].Data());
    fTotalJESUnc.push_back( new JetCorrectionUncertainty(*p) );
  }

  //muon energy corrector
  fMuCor = getMuonCorrector(jecDir,url);

  //control the sec vtx analysis
  LxyAnalysis lxyAn;

  if(weightsFile.size()) {
    TString burl(weightsFile[0].c_str());
    burl += "/BfragWeights.root";
    fBfragWgt = new BfragWeighter( burl );
  }

  //prepare the output file
  TString proctag=gSystem->BaseName(url);
  Ssiz_t pos=proctag.Index(".root");
  proctag.Remove(pos,proctag.Length());
  TString outUrl(out);
  gSystem->ExpandPathName(outUrl);
  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += proctag;
  outUrl +=  "_filt"; outUrl += reqControlType;
  outUrl += ".root";
  TFile *spyFile=TFile::Open(outUrl, "recreate");
  spyFile->cd();
  TDirectory *spyDir=0;
  if(saveSummaryTree)
    {
      gSystem->Exec("mkdir -p " + out);
      gDirectory->SaveSelf();
      spyFile->rmdir(proctag);
      spyDir = spyFile->mkdir("dataAnalyzer");
      spyDir->cd();
      lxyAn.attachToDir(spyDir);
    }

  //
  // check input file
  //
  TFile *inF = TFile::Open(url);
  if(inF==0) return -1;
  if(inF->IsZombie()) return -1;

  //
  // pileup reweighter
  //
  std::vector<double> dataPileupDistributionDouble = runProcess.getParameter< std::vector<double> >("datapileup");
  std::vector<float> dataPileupDistribution; for(unsigned int i=0;i<dataPileupDistributionDouble.size();i++){dataPileupDistribution.push_back(dataPileupDistributionDouble[i]);}
  std::vector<float> mcPileupDistribution;
  if(isMC){
    TString puDist(baseDir+"/pileup");
    TH1F* histo = (TH1F *) inF->Get(puDist);
    if(!histo)std::cout<<"pileup histogram is null!!!\n";
    for(int i=1;i<=histo->GetNbinsX();i++){mcPileupDistribution.push_back(histo->GetBinContent(i));}
    delete histo;
  }
  while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
  while(mcPileupDistribution.size()>dataPileupDistribution.size())dataPileupDistribution.push_back(0.0);
  gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
  edm::LumiReWeighting *LumiWeights= isMC ? new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution): 0;

  //book histograms
  SmartSelectionMonitor controlHistos;
  TH1 *Hcutflow=(TH1 *)controlHistos.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;
  for(int ibin=1; ibin<=6; ibin++) Hcutflow->SetBinContent(ibin,1);
  controlHistos.addHistogram( new TH1F("pthat", ";#hat{p}_{T} [GeV]; Events",100,0,1500) );
  controlHistos.addHistogram( new TH1F ("nvertices", "; Vertex multiplicity; Events", 50, 0.,50.) );


  controlHistos.addHistogram( new TH1F ("balance", ";Balance=p_{T}(probe)/p_{T}(tag); Events", 50,0.,2.) );
  controlHistos.addHistogram( new TH1F ("mass", ";Mass; Events", 50,0.,150.) );
  for(size_t i=0; i<2; i++)
    {
      for(size_t j=0; j<2; j++)
	{
	  TString pf(i==0 ? "" : "raw");
	  pf+=(j==0 ? "tag" : "probe");
	  controlHistos.addHistogram( new TH1F (pf+"pt", "; Transverse momentum [GeV]; Events", 50,0.,maxPtForTagWgt) );
	  controlHistos.addHistogram( new TH1F (pf+"eta", "; Pseudo-rapidity; Events", 50, 0.,5) );
	}
    }

  //jet control
  TString jetFlavors[]={"","b","udsg","c","unmatched"};
  size_t nJetFlavors=sizeof(jetFlavors)/sizeof(TString);
  TH1 *hflav=controlHistos.addHistogram( new TH1F ("probeflav", ";Probe flavour; Jets", nJetFlavors, 0.,nJetFlavors) );
  for(int xbin=1; xbin<=hflav->GetXaxis()->GetNbins(); xbin++) hflav->GetXaxis()->SetBinLabel(xbin,jetFlavors[xbin-1]);

  ///
  // process events file
  //
  DataEventSummaryHandler evSummary;
  if( !evSummary.attach( (TTree *) inF->Get(baseDir+"/data") ) )  { inF->Close();  return -1; }
  const Int_t totalEntries=evSummary.getEntries();

  float cnorm=1.0;
  if(isMC){
    TH1F* cutflowH = (TH1F *) inF->Get(baseDir+"/cutflow");
    if(cutflowH) cnorm=cutflowH->GetBinContent(1);
  }
  Hcutflow->SetBinContent(1,cnorm);

  cout << "Processing: " << proctag << " @ " << url << endl
       << "Initial number of events: " << cnorm << endl
       << "Events in tree:           " << totalEntries << endl
       << " xSec x BR:               " << xsec << endl;

  //
  // analyze (puf...)
  //
  for (int inum=0; inum < totalEntries; ++inum)
    {
      if(inum%500==0) { printf("\r [ %d/100 ]",int(100*float(inum)/float(totalEntries))); cout << flush; }
      evSummary.getEntry(inum);
      DataEventSummary &ev = evSummary.getEvent();

      if(isV0JetsMC && ev.nup>5) continue;

      //pileup weight
      float weightNom(1.0);
      if(LumiWeights) 	weightNom  = LumiWeights->weight(ev.ngenITpu);

      //get the trigger
      float triggerPrescale(1.0),triggerThreshold(0);
      bool dijetTrigger(false), gammaTrigger(false), eTrigger(false), muTrigger(false), eeTrigger(false), mumuTrigger(false);
      if(runDijetControl)
	{
	  triggerThreshold=30;
	  for(int itrig=0; itrig<=4; itrig++)
	    {
	      if(!ev.t_bits[itrig]) continue;
	      dijetTrigger=true;
	      triggerPrescale=ev.t_prescale[itrig];
	      if(itrig==1) triggerThreshold=50;
	      if(itrig==2) triggerThreshold=80;
	      if(itrig==3) triggerThreshold=120;
	      if(itrig==4) triggerThreshold=320;
	      break;
	    }
	  //if(isMC) dijetTrigger=true;
	  if(!dijetTrigger) continue;
	}
      else if(runGammaControl)
	{
	  for(size_t itrig=10; itrig>=7; itrig--)
	    {
	      if(!ev.t_bits[itrig]) continue;
	      gammaTrigger=true;
	      triggerPrescale=ev.t_prescale[itrig];
	      if(itrig==10) triggerThreshold=92; //90
	      if(itrig==9)  triggerThreshold=77; //75
	      if(itrig==8)  triggerThreshold=50;
	      if(itrig==7)  triggerThreshold=36;
	      break;
	    }
	  //if(isMC) {gammaTrigger=true;}
	  if(!gammaTrigger) continue;
	}
      else if(runWControl)
	{
	  eTrigger  = ev.t_bits[13];
	  muTrigger = ev.t_bits[6];
	  if(!eTrigger && !muTrigger) continue;
	}
      else if(runZControl)
	{
	  eeTrigger   = ev.t_bits[0];
	  mumuTrigger = ev.t_bits[2] || ev.t_bits[3];
	  if(!eeTrigger && !mumuTrigger) continue;
	}


      //weight for the event
      float weight(weightNom*triggerPrescale);

      //select the objects
      data::PhysicsObjectCollection_t rawLeptons = evSummary.getPhysicsObject(DataEventSummaryHandler::LEPTONS);
      data::PhysicsObjectCollection_t leptons    = top::selectLeptons(rawLeptons, fMuCor, ev.rho, isMC);
      data::PhysicsObjectCollection_t rawPhotons = evSummary.getPhysicsObject(DataEventSummaryHandler::PHOTONS);
      data::PhysicsObjectCollection_t photons    = top::selectPhotons(rawPhotons, triggerThreshold, ev.rho);
      data::PhysicsObjectCollection_t rawJets    = evSummary.getPhysicsObject(DataEventSummaryHandler::JETS);
      LorentzVector jetDiff                      = utils::cmssw::updateJEC(rawJets, fJesCor, fTotalJESUnc, ev.rho, ev.nvtx, isMC);
      data::PhysicsObjectCollection_t jets       = top::selectJets(rawJets,leptons);
      if(runDijetControl)
	{
	  data::PhysicsObjectCollection_t dijetCtrlJets;
	  for(size_t ijet=0; ijet<jets.size(); ijet++)
	    {
	      if(jets[ijet].pt()<triggerThreshold || fabs(jets[ijet].eta())>2.4) continue;
	      dijetCtrlJets.push_back( jets[ijet] );
	    }
	}

      data::PhysicsObjectCollection_t recoMet    = evSummary.getPhysicsObject(DataEventSummaryHandler::MET);
      recoMet[2].SetPxPyPzE(recoMet[2].px()-jetDiff.px(),
			    recoMet[2].py()-jetDiff.py(),
			    0.,
			    sqrt(pow(recoMet[2].px()-jetDiff.px(),2)+pow(recoMet[2].py()-jetDiff.py(),2))
			    );
      std::vector<LorentzVector> met              = utils::cmssw::getMETvariations(recoMet[2], jets, leptons, isMC);
      data::PhysicsObjectCollection_t pf          = evSummary.getPhysicsObject(DataEventSummaryHandler::PFCANDIDATES);

      //build the analysis box
      ControlBox box = assignBox(reqControlType,leptons,photons,jets,met,pf);
      if(box.cat==0) continue;
      data::PhysicsObjectCollection_t gen;
      if(isMC){
	gen=evSummary.getPhysicsObject(DataEventSummaryHandler::GENPARTICLES);
	box.assignProbeFlavour(gen);
      }

      //efficiency weighting
      if(isMC)          weight *= box.selEffSF;

      //tag pT weighing
      float tagWeight(1.0);
      /*
      if(tagPtWeightGr) tagWeight=tagPtWeightGr->Eval(TMath::Min(box.tag->pt(),maxPtForTagWgt-5.0));
      weight *= tagWeight;
      */

      //gen control
      std::vector<TString> catsToFill(1,"all");
      if(isMC) controlHistos.fillHisto("pthat",catsToFill,ev.pthat,weight);

      //vertices control
      controlHistos.fillHisto("nvertices",catsToFill,ev.nvtx,weight);

      //tag/probe kinematics
      controlHistos.fillHisto("balance",   catsToFill, box.probe->pt()/box.tag->pt(),  weight);

      //before tag pT weighting, for control
      controlHistos.fillHisto("rawtagpt",     catsToFill, box.tag->pt(),               weight/tagWeight);
      controlHistos.fillHisto("rawtageta",    catsToFill, fabs(box.tag->eta()),        weight/tagWeight);
      controlHistos.fillHisto("rawprobept",   catsToFill, box.probe->pt(),             weight/tagWeight);
      controlHistos.fillHisto("rawprobeeta",  catsToFill, fabs(box.probe->eta()),      weight/tagWeight);

      //after tag pT weighting
      controlHistos.fillHisto("mass",      catsToFill, box.tag->mass(),                weight);
      controlHistos.fillHisto("tagpt",     catsToFill, box.tag->pt(),                  weight);
      controlHistos.fillHisto("tageta",    catsToFill, fabs(box.tag->eta()),           weight);
      controlHistos.fillHisto("probept",   catsToFill, box.probe->pt(),                weight);
      controlHistos.fillHisto("probeeta",  catsToFill, fabs(box.probe->eta()),         weight);
      controlHistos.fillHisto("probeflav", catsToFill, 0,                              weight);
      if(isMC) 	controlHistos.fillHisto("probeflav", catsToFill, box.probeFlav,        weight);

      //secondary vertex characteristics evaluated for the probe
      const data::PhysicsObject_t &svx = box.probe->getObject("svx");
      float lxy=svx.vals.find("lxy")->second;
      if(lxy)
	{
	  catsToFill.clear();
	  catsToFill.push_back("svx");

	  controlHistos.fillHisto("probeflav", catsToFill, 0,                              weight);
	  if(isMC) 	controlHistos.fillHisto("probeflav", catsToFill, box.probeFlav,        weight);
	}


      //save selected event
      if(!saveSummaryTree) continue;
      lxyAn.resetBeautyEvent();
      BeautyEvent_t &bev=lxyAn.getBeautyEvent();
      bev.evcat=box.cat;
      bev.gevcat=box.cat;
      bev.run=ev.run;
      bev.event=ev.event;
      bev.lumi=ev.lumi;
      bev.nvtx=ev.nvtx;
      bev.rho=ev.rho;
      bev.nw=1;
      bev.w[0]=weight;
      bev.qscale=ev.qscale;  bev.x1=ev.x1; bev.x2=ev.x2; bev.id1=ev.id1; bev.id2=ev.id2;
      std::vector<data::PhysicsObject_t *> boxtag,boxprobe;
      std::vector<data::PhysicsObject_t *> fwdjets;
      boxtag.push_back( box.tag );
      boxprobe.push_back( box.probe );
      lxyAn.analyze( boxtag, boxprobe, fwdjets, met, pf, gen);
      //add fragmentation weights using matched b's and B hadrons
      if(fBfragWgt) {
	for(Int_t ij=0; ij<bev.nj; ij++)
	  {
	    if(abs(bev.bid[ij])!=5 || bev.gjpt[ij]<=0 || bev.bhadpt[ij]<=0) continue;
	    std::vector<float> bfragWeights=fBfragWgt->getEventWeights( bev.bhadpt[ij]/bev.gjpt[ij] );
	    bev.bwgt[ij][0]=bfragWeights[0];
	    bev.bwgt[ij][1]=bfragWeights[1];
	    bev.bwgt[ij][2]=bfragWeights[2];
	    bev.bwgt[ij][3]=bfragWeights[3];
	    bev.bwgt[ij][4]=bfragWeights[4];
	    bev.bwgt[ij][5]=bfragWeights[5];
	  }
      }


      spyDir->cd();
      lxyAn.save();
    }

  inF->Close();

  //
  // save histos to local file
  //
  spyFile->cd();
  spyFile->Write();
  TVectorD constVals(2);
  constVals[0] = isMC ? cnorm : 1.0;
  constVals[1] = isMC ? xsec  : 1.0;
  constVals.Write("constVals");
  controlHistos.Write();
  spyFile->Close();
 }
