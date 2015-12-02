#! /usr/bin/env python
import ROOT
import utils
from array import array
import optparse
import numpy as np
import os,sys
from UserCode.TopMassSecVtx.storeTools_cff import fillFromStore
from makeSVLMassHistos import getMassTrees
from runDileptonUnfolding import *
from UserCode.TopMassSecVtx.PlotUtils import *


"""
Get first three Mellin moments from each of the distributions
"""
def getMoments(root,v,moments,syst,out):

    histos={}
    
    #open ROOT file
    fIn=ROOT.TFile.Open(root)

    var = str(v)
    #unfolded histograms generated level
    if moments == 'unf':
      unfolded       = getPathToObjects(fIn.Get('%s_unfolded'%var))
      key = 'unf'
      for p in unfolded:
         h=fIn.Get(p)
         hname=h.GetName()
         if 'data_unfolded' in hname : 
             histos['unf'] = h.Clone('data_unfolded')   

    if moments == 'gen':
      genList       = getPathToObjects(fIn.Get('%s_gen'%(var)))
      key = 'gen'
      for p in genList:
          h=fIn.Get(p)
          hname=h.GetName()
          if 'fromfiles' in out:
            if not 'signal_gen' in histos : 
                histos['gen']=h.Clone('gen')
            else : 
                histos['gen'].Add(h)
          else:
           if 'TTJets' in hname:
            if not 'gen' in histos : 
                histos['gen']=h.Clone('gen')
            else : 
                histos['gen'].Add(h)

    if moments == 'rec':
      recList       = getPathToObjects(fIn.Get('%s_rec'%(var)))
      key = 'signal'
      for p in recList:
          h=fIn.Get(p)
          hname=h.GetName()
          if 'fromfiles' in out:
            if not 'signal' in histos : 
                histos['signal']=h.Clone('signal')
            else : 
                histos['signal'].Add(h)
          else:
           if 'TTJets' in hname:
            if not 'signal' in histos : 
                histos['signal']=h.Clone('signal')
            else : 
                histos['signal'].Add(h)
           else:
            if not 'bkg' in histos : histos['bkg']=h.Clone('bkg')
            else : histos['bkg'].Add(h)


    #detach objects from file and close file
    for h in histos: histos[h].SetDirectory(0)
    fIn.Close()

    mean = histos[key].GetMean()
    meanerr = histos[key].GetMeanError()
    rms = histos[key].GetRMS()
    rmserr = histos[key].GetRMSError()

    u1 = mean
    u1err = meanerr
    u2 = mean**2 + rms**2
    u2err = sqrt((2*rms*rmserr)**2+(2*mean*meanerr)**2)

    #Create a text file with the three different moments for the corresponding mass
    outfileName = out+'moments_' + var + '_' + moments+'.txt'
    outfile = open(outfileName,'a')
    outfile.write('%s %4.6f %4.6f %4.6f %4.6f \n' %(syst,u1,u2))
    outfile.close()

"""
Loop over a tree and fill histograms you declared before
"""
def createHistosSyst(var,filename,isData,histos,syst):
    
    #Getting histograms labeling
    rec = var+'_rec'
    wgt = rec+'_wgt'
    gen = var+'_gen'
    wgt_gen = gen+'_wgt'
    mig = var+'_migration'
    
    #open file
    fIn=ROOT.TFile.Open(filename)
    
    #loop over events in the tree and fill histos
    tree=fIn.Get('DileptonInfo')
    print tree.GetEntriesFast()
    for i in xrange(0,tree.GetEntriesFast()):
        tree.GetEntry(i)

        #select only emu events

        if tree.EvCat != -11*13 : 
            continue
        if not isData: 
           if tree.GenLpPt == 0 or tree.GenLmPt == 0: continue

        baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4] 
        print baseWeight

        #base weight: BR fix for ttbar x pileup x lepton selection x xsec weight
        if syst == 'lesup' or syst == 'lesdn' or syst is None: baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4] 
        #Pileup up puup
        if syst == 'puup': baseWeight = tree.Weight[0]*tree.Weight[2]*tree.Weight[4] 
        #Pileup dpwn pudn
        if syst == 'pudn': baseWeight = tree.Weight[0]*tree.Weight[3]*tree.Weight[4]
        #Lepton selection up lepselup
        if syst == 'lepselup': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[5] 
        #Lepton selection down lepseldn
        if syst == 'lepseldn': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[6] 
        #Top p_{T} weight applied toppt
        if syst == 'toppt': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.Weight[10] 
        #Uncl. MET up umetup
        #if syst == 'umetup': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.METWeight[1] 
        #Uncl. MET down umetdn
        #if syst == 'umetdn': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.METWeight[2] 
        #b-tag eff up btagup
        # if syst == 'btagup': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.BtagWeight[1]
        #b-tag eff dn btagdn
        # if syst == 'btagdn': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.BtagWeight[2]
        #Jet energy scale up jesup
        # if syst == 'jesup': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.JESWeight[1]
        #Jet energy scale down jesdn
        # if syst == 'jesdn': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.JESWeight[2]
        #Jet energy resolution up jerup
        # if syst == 'jerup': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.JESWeight[3]
        #Jet energy resolution down jerdn
        # if syst == 'jerdn': baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4]*tree.JESWeight[4]

        #event weight
        weight = 1 if isData else baseWeight
        
        #positive lepton
        lp=ROOT.TLorentzVector()
        if syst == 'lesdn':
            #lesdn
            lp.SetPtEtaPhiM(tree.LpPt_sf[0],tree.LpEta,tree.LpPhi,0.)
        elif syst == 'lesup':
            #lesUp
            lp.SetPtEtaPhiM(tree.LpPt_sf[1],tree.LpEta,tree.LpPhi,0.)
        else:
            lp.SetPtEtaPhiM(tree.LpPt,tree.LpEta,tree.LpPhi,0.)
        glp=ROOT.TLorentzVector()
        glp.SetPtEtaPhiM(tree.GenLpPt,tree.GenLpEta,tree.GenLpPhi,0.)

        #negative lepton
        lm=ROOT.TLorentzVector()
        if syst == 'lesdn':
            #lesdn
            lm.SetPtEtaPhiM(tree.LmPt_sf[0],tree.LmEta,tree.LmPhi,0.) 
        elif syst == 'lesup':
            #lesUp
            lm.SetPtEtaPhiM(tree.LmPt_sf[1],tree.LmEta,tree.LmPhi,0.)
        else:
            lm.SetPtEtaPhiM(tree.LmPt,tree.LmEta,tree.LmPhi,0.) 
        glm=ROOT.TLorentzVector()
        glm.SetPtEtaPhiM(tree.GenLmPt,tree.GenLmEta,tree.GenLmPhi,0.)

        #charged lepton pair - pt
        ll=ROOT.TLorentzVector()
        ll = lp + lm
        gll=ROOT.TLorentzVector()
        gll = glp + glm

        #fill the histograms according to the distrubution variable
        #pT positive lepton
        if var == 'ptpos': 
            histos[rec].Fill(lp.Pt())
            binWidth = histos[wgt].GetXaxis().GetBinWidth(histos[wgt].GetXaxis().FindBin(lp.Pt() ) )
            histos[wgt].Fill(lp.Pt(),weight/binWidth)
            if not isData:
                    histos[gen].Fill(glp.Pt(),weight)
                    binWidthGen = histos[wgt_gen].GetXaxis().GetBinWidth(histos[wgt_gen].GetXaxis().FindBin(glp.Pt() ) )
                    histos[wgt_gen].Fill(glp.Pt(),weight/binWidthGen)
                    histos[mig].Fill(glp.Pt(),lp.Pt(),weight)

        #Second distribution: Pt(l+l-) = ll.Pt      
        if var == 'ptll': 
            histos[rec].Fill(ll.Pt(),weight)
            binWidth = histos[wgt].GetXaxis().GetBinWidth(histos[wgt].GetXaxis().FindBin(ll.Pt() ) )
            histos[wgt].Fill(ll.Pt(),weight/binWidth)
            if not isData:
                    histos[gen].Fill(gll.Pt(),weight)
                    binWidthGen = histos[wgt_gen].GetXaxis().GetBinWidth(histos[wgt_gen].GetXaxis().FindBin(gll.Pt() ) )
                    histos[wgt_gen].Fill(gll.Pt(),weight/binWidthGen)
                    histos[mig].Fill(gll.Pt(),ll.Pt(),weight)

        #Third distribution: M(l+l-) = ll.M
        if var == 'mll': 
            histos[rec].Fill(ll.M(),weight)
            binWidth = histos[wgt].GetXaxis().GetBinWidth(histos[wgt].GetXaxis().FindBin(ll.M() ) )
            histos[wgt].Fill(ll.M(),weight/binWidth)
            if not isData:
                    histos[gen].Fill(gll.M(),weight)
                    binWidthGen = histos[wgt_gen].GetXaxis().GetBinWidth(histos[wgt_gen].GetXaxis().FindBin(gll.M() ) )
                    histos[wgt_gen].Fill(gll.M(),weight/binWidthGen)
                    histos[mig].Fill(gll.M(),ll.M(),weight)

        #Fourth distribution: E(l+)+E(l-) = lp.E() + lm.E()
        if var == 'EposEm': 
            histos[rec].Fill(lp.E() + lm.E(),weight)
            binWidth = histos[wgt].GetXaxis().GetBinWidth(histos[wgt].GetXaxis().FindBin(lp.E() + lm.E() ) )
            histos[wgt].Fill(lp.E() + lm.E(),weight/binWidth)
            if not isData:
                    histos[gen].Fill(glp.E() + glm.E(),weight)
                    binWidthGen = histos[wgt_gen].GetXaxis().GetBinWidth(histos[wgt_gen].GetXaxis().FindBin(glp.E()+glm.E() ) )
                    histos[wgt_gen].Fill(glp.E()+glm.E(),weight/binWidthGen)
                    histos[mig].Fill(glp.E() + glm.E(),lp.E() + lm.E(),weight)

        #Fifth distribution: Pt(l+)+Pt(l-) = lp.Pt() + lm.Pt()
        if var == 'ptposptm': 
            histos[rec].Fill(lp.Pt() + lm.Pt(),weight)
            binWidth = histos[wgt].GetXaxis().GetBinWidth(histos[wgt].GetXaxis().FindBin(lp.Pt() + lm.Pt() ) )
            histos[wgt].Fill(lp.Pt() + lm.Pt(),weight/binWidth)
            if not isData:
                    histos[gen].Fill(glp.Pt() + glm.Pt(),weight)
                    binWidthGen = histos[wgt_gen].GetXaxis().GetBinWidth(histos[wgt_gen].GetXaxis().FindBin(glp.Pt()+glm.Pt() ) )
                    histos[wgt_gen].Fill(glp.Pt()+glm.Pt(),weight/binWidthGen)
                    histos[mig].Fill(glp.Pt() + glm.Pt(),lp.Pt() + lm.Pt(),weight)

    #close file
    fIn.Close()
      
"""
Run anaylsis for nominal sample
"""
def runWeights(opt):

    # Nominal for weight based systematics

    # define histograms 
    for filename in os.listdir(opt.nominal):
        if not os.path.splitext(filename)[1] == '.root': continue   
        isData = True if 'Data' in filename else False
        #if 'MC8TeV' in filename:
        if 'MC8TeV_TTJets_MSDecays_172v5' in filename: 
    
            filename = opt.nominal + filename
            print filename
            #define histograms
            histos = getAnalysisHistograms(opt.var)

            syst = opt.syst
            var = opt.var

            #filling histograms
            createHistosSyst(var,filename,isData,histos,syst)

            #dump histograms just for this file in the outDir
            fOut=ROOT.TFile.Open(os.path.join(opt.output,os.path.basename(filename)),'RECREATE')
            for h in histos: histos[h].Write()
            print 'Histograms saved in %s' % fOut.GetName()
            fOut.Close()

"""
Run anaylsis for systematics files
"""
def runFromFiles(opt):

    var = opt.var

    # Nominal for systematics from files

    # define histograms 
    if opt.syst == 'scaleup': files = ['MC8TeV_TTJets_MSDecays_scaleup.root']
    if opt.syst == 'scaledown': files = ['MC8TeV_TTJets_MSDecays_scaledown.root']
    if opt.syst == 'tchscaleup': files = ['MC8TeV_SingleT_t_scaleup.root','MC8TeV_SingleTbar_t_scaleup.root']
    if opt.syst == 'tchscaledown': files = ['MC8TeV_SingleT_t_scaledown.root','MC8TeV_SingleTbar_t_scaledown.root']
    if opt.syst == 'twchscaleup': files = ['MC8TeV_SingleT_tW_scaleup.root','MC8TeV_SingleTbar_tW_scaleup.root']
    if opt.syst == 'twchscaledown': files = ['MC8TeV_SingleT_tW_scaledown.root','MC8TeV_SingleTbar_tW_scaledown.root']
    if opt.syst == 'matchingup': files = ['MC8TeV_TTJets_MSDecays_matchingup.root']
    if opt.syst == 'matchingdown': files = ['MC8TeV_TTJets_MSDecays_matchingdown.root']
    if opt.syst == 'p11': files = ['MC8TeV_TTJets_TuneP11.root']
    if opt.syst == 'p11nocr': files =  ['MC8TeV_TTJets_TuneP11noCR.root', 'MC8TeV_TTJets_SemiLep_TuneP11noCR.root']
    if opt.syst == 'p11mpihi': files =  ['MC8TeV_TTJets_TuneP11mpiHi.root', 'MC8TeV_TTJets_SemiLep_TuneP11mpiHi.root']
    if opt.syst == 'powherw': files =  ['MC8TeV_TT_AUET2_powheg_herwig.root']


    for filename in files: 
        isData = True if 'Data' in filename else False

        filename = opt.fromfiles + filename
        print filename
        #define histograms
        histos = getAnalysisHistograms(var)

        syst = None

        #filling histograms
        createHistos(var,filename,isData,histos)

        #dump histograms just for this file in the outDir
        fOut=ROOT.TFile.Open(os.path.join(opt.output,os.path.basename(filename)),'RECREATE')
        for h in histos: histos[h].Write()
        print 'Histograms saved in %s' % fOut.GetName()
        fOut.Close()



"""
steer
"""
def main():
	usage = 'usage: %prog [options]'
	parser = optparse.OptionParser(usage)
	parser.add_option('-i', '--input',
                          dest='input',   
                          default='/store/cmst3/group/top/summer2015/treedir_bbbcb36/ttbar/', 
                          help='input directory with the files [default: %default]')
	parser.add_option('-r', '--root',
                          dest='root', 
                          default=None,
                          help='ROOT file with distributions to unfold [default: %default]')
	parser.add_option('-v', '--var',
                          dest='var', 
                          default='ptpos',
                          help='Variable to unfold (note requires var_rec,var_gen,var_migration plots stored in the ROOT file [default: %default]')
        parser.add_option('-g',
                          dest='migration', 
                          default=None,
                          help='Work with migration matrix from nominal sample [default: %default]')
	parser.add_option('--jobs',
                          dest='jobs', 
                          default=1,
                          type=int,
                          help='# of jobs to process in parallel the trees [default: %default]')
        parser.add_option('-m',
                          dest='ma', 
                          default='f',
                          help='Analysis for mass files [default: %default]')
	parser.add_option('-o', '--output',
                          dest='output', 
                          default='unfoldResults',                                                                       
                          help='Output directory [default: %default]')
        parser.add_option('-f',
                          dest='fromfiles', 
                          default=None,
                          help='Analysis syst from files [default: %default]')
        parser.add_option('-n',
                          dest='nominal', 
                          default=None,
                          help='Analysis for nominal sample [default: %default]')
        parser.add_option('-s',
                          dest='syst', 
                          default=None,
                          help='Analysis for systematic uncertainty [default: %default]')

	(opt, args) = parser.parse_args()

	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptTitle(0)
	ROOT.gROOT.SetBatch(True)
	setTDRStyle()
        ROOT.gSystem.Load("libUserCodeTopMassSecVtx")
	ROOT.AutoLibraryLoader.enable()
	os.system('mkdir -p %s' % opt.output)

	# Check if one needs to create a new workspace or run pseudo-experiments	
	if opt.root is None :
            if opt.fromfiles is None:
                if opt.nominal is None:
                    print 80*'-'
                    print 'Creating ROOT file with migration matrices, data and background distributions of %s from %s'%(opt.var,opt.input)
                    createSummaryTasks(opt)
                    print 80*'-'
                else:
                    print 80*'-'
                    print 'Creating ROOT file with migration matrices, data and background distributions of %s for nominal sample from %s'%(opt.var,opt.nominal)
                    runWeights(opt)
                    print 80*'-'
            else:
                print 80*'-'
                print 'Doing analysis of %s from %s'%(opt.var,opt.fromfiles)
                runFromFiles(opt)
                print 80*'-'
        else:
             print 80*'-'
             print 'Unfolding variable %s from %s'%(opt.var,opt.root)
             unfoldVariable(opt)
             print 80*'-'
        return 0

if __name__ == "__main__":
	sys.exit(main())
