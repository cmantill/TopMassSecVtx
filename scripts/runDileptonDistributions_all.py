#! /usr/bin/env python
import ROOT
from array import array
import optparse
import numpy as np
import os,sys
from UserCode.TopMassSecVtx.storeTools_cff import fillFromStore
from makeSVLMassHistos import addTopMassTreesFromDir
from makeSVLMassHistos import getMassTrees
from copy import deepcopy
from runPlotter import runPlotter, addPlotterOptions

from UserCode.TopMassSecVtx.PlotUtils import *

"""
Returns histograms to be filled in the loop, depending on the distribution variable you chose to work with
10 bins gen level
41 bins rec level
"""
def getAnalysisHistograms(var) :

    histos={}
  
    # pT charged-lepton pair
    if var == 'ptll': 
        bins_gen = [15.13,30,41.76,51.57,60.23,70.22,79.93,89.23,104.48,136.92] 
        bins_rec = [8,10.79,15.13,19.47,23.05,26.53,30,33.03,36.06,39.08,41.76,
                    44.27,46.79,49.31,51.57,53.73,55.88,58.04,60.23,62.73,65.23,
                    67.73,70.22,72.65,75.07,77.5,79.93,82.25,84.58,86.9,89.23,
                    92.32,95.79,99.26,104.48,110.19,116.54,125.0,136.92,159.29,250]    
        title = ';p_{T}(l^{+}l^{-}) [GeV];Events'

    # E positive lepton
    if var == 'Epos':
        bins_gen= [24.96,36.20,45.21,54.26,64.11,75.28,89.10,107.18,133.67,183.27]
        bins_rec = [24.96,28.42,31.33,33.76,36.20,38.54,40.76,42.99,45.21,
                    47.45,49.68,51.92,54.26,56.63,59.00,61.51,64.11,66.72,
                    69.51,72.38,75.28,78.55,81.83,85.42,89.10,93.15,97.33,
                    102.12,107.18,112.63,118.81,125.75,133.67,142.84,153.71,
                    166.59,183.27,205.40,238.91]
        title = ';E(l^{+}l^{-}) [GeV];Events'

    # pT positive lepton                                                                                                                                                                                                    
    if var == 'ptpos':
        bins_gen = [22.72,28.18,33.67,39.17,44.89,50.94,59.19,68.15,83.59,115.65]                                                                                                                                                
        bins_rec = [20,21.36,22.73,24.09,25.45,26.82,28.18,29.55,30.92,32.29,
                    33.67,35.04,36.41,37.79,39.17,40.57,42.01,43.45,44.89,46.33,
                    47.77,49.21,50.93,53,55.06,57.12,59.18,61.37,63.63,65.89,
                    68.15,70.73,74.76,78.74,83.59,88.75,95,102.85,115.65,138.46,250]                                                                                                                                            
        title = ';p_{T}(l^{+}) [GeV];Events'

    # M charged-lepton pair                                                                                                                                                                                                   
    if var == 'mll':
        bins_gen = [30.85,48.96,62.37,75.26,88.33,102.76,117.89,136.41,167.91,235.33]                                                                                                                                           
        bins_rec = [20,25.63,30.84,36.06,40.54,44.75,48.96,52.48,55.78,59.07,
                    62.37,65.6,68.82,72.03,75.26,78.51,81.75,84.99,88.34,91.98,
                    95.63,99.27,102.75,106.19,109.62,113.25,117.89,122.52,127.14,131.78,
                    136.41,142.28,148.52,157.31,167.91,180.29,193.4,214.25,235.33,289.1,450]                                                                                                                                    
        title = ';M(l^{+}l^{-}) [GeV];Events'

    # Scalar sum of E                                                                                                                                                                                                         
    if var == 'EposEm':
        bins_gen = [75.17,95.35,114.58,130.11,146.95,166.87,188.33,223.57,271.12,364.29]                                                                                                                                         
        bins_rec = [50,66.05,75.17,80.82,86.46,91.02,95.35,99.67,104.8,110.02,
                    114.6,118.6,122.55,126.42,130.11,133.81,137.5,142.22,146.95,151.8,
                    156.88,161.95,166.87,171.76,176.91,182.55,188.33,195.1,202.61,212.06,
                    223.57,233.28,244.76,257.62,271.12,288.75,307.35,328.98,364.29,420.14,500]                                                                                                                               
        title = ';E(l^{+})+E(l^{-}) [GeV];Events'

    # Scalar sum of Pt                                                                                                                                                                                                        
    if var == 'ptposptm':
        bins_gen = [55.73,68.94,80.52,89.01,97.64,107.1,118.8,134.43,154.04,204.56]                                                                                                                                              
        bins_rec = [46,50.15,55.73,60.65,63.41,66.18,68.94,71.83,74.79,77.76,
                    80.52,82.64,84.76,86.88,89.01,91.15,93.32,95.48,97.64,99.81,
                    102.22,104.66,107.1,109.53,112.53,115.67,118.8,122.54,126.65,130.7,
                    134.43,138.17,142.8,148.28,154.04,159.91,170.92,185.17,204.56,236.2,275]                                                                                                                                 
        title = ';p_{T}(l^{+})+p_{T}(l^{-}) [GeV];Events'

    # Labeling histograms
    rec = var+'_rec'
    gen = var+'_gen'

    # Declaring histos
    histos[rec]=ROOT.TH1F(rec,title,len(bins_rec)-1,array('d',bins_rec))
    histos[gen]=ROOT.TH1F(gen,title,len(bins_gen)-1,array('d',bins_gen))

    for h in histos:
        histos[h].Sumw2()
        histos[h].SetDirectory(0)

    return histos

"""
Loop over a tree and fill histograms you declared before
"""
def createHistos(filename,isData,histos,var):

    # Getting histograms labeling
    rec = var+'_rec'
    gen = var+'_gen'
    
    # Open file
    fIn=ROOT.TFile.Open(filename)
    
    # Loop over events in the tree and fill histos
    tree=fIn.Get('DileptonInfo')
    
    for i in xrange(0,tree.GetEntriesFast()):
        tree.GetEntry(i)

        # Select only emu events
        if tree.EvCat != -11*13 : continue
        if not isData: 
           if tree.GenLpPt == 0 or tree.GenLmPt == 0: continue

        # Base weight: BR fix for ttbar x pileup x lepton selection x xsec weight
        baseWeight = tree.Weight[0]*tree.Weight[1]*tree.Weight[4] #*tree.XSWeight
                        
        # Event weight
        weight = 1 if isData else baseWeight
        
        # Positive lepton
        lp=ROOT.TLorentzVector()
        lp.SetPtEtaPhiM(tree.LpPt,tree.LpEta,tree.LpPhi,0.)
        glp=ROOT.TLorentzVector()
        glp.SetPtEtaPhiM(tree.GenLpPt,tree.GenLpEta,tree.GenLpPhi,0.)

        # Negative lepton
        lm=ROOT.TLorentzVector()
        lm.SetPtEtaPhiM(tree.LmPt,tree.LmEta,tree.LmPhi,0.)       
        glm=ROOT.TLorentzVector()
        glm.SetPtEtaPhiM(tree.GenLmPt,tree.GenLmEta,tree.GenLmPhi,0.)

        # Charged lepton pair 
        ll=ROOT.TLorentzVector()
        ll = lp + lm
        gll=ROOT.TLorentzVector()
        gll = glp + glm

        # Fill the histograms according to the distrubution variable
        # pT of the lepton pair: Pt(l+l-) = ll.Pt      
        if var == 'ptll': 
            histos[rec].Fill(ll.Pt(),weight)
            if not isData:
                histos[gen].Fill(gll.Pt(),weight)

        # E positive lepton: E(l+) = lp.E()                                                                                                                                                                                   
        if var == 'Epos':
            histos[rec].Fill(lp.E(),weight)
            if not isData:
                histos[gen].Fill(glp.E(),weight)

        # pT positive lepton: Pt(l+) = lp.Pt()                                                                                                                                                                                  
        if var == 'ptpos':
            histos[rec].Fill(lp.Pt(),weight)                                                                                                                                                                                    
            if not isData:
                histos[gen].Fill(glp.Pt(),weight)

        # Invariant mass of the lepton pair: M(l+l-) = ll.M()                                                                                                                                                                    
        if var == 'mll':
            histos[rec].Fill(ll.M(),weight)                                                                                                                                                                                     
            if not isData:
                    histos[gen].Fill(gll.M(),weight)

        # Sum of energies: E(l+)+E(l-) = lp.E() + lm.E()                                                                                                                                                                    
        if var == 'EposEm':
            histos[rec].Fill(lp.E() + lm.E(),weight)                                                                                                                                                                            
            if not isData:
                histos[gen].Fill(glp.E() + glm.E(),weight)
                
        # Sum of Pt: Pt(l+)+Pt(l-) = lp.Pt() + lm.Pt()                                                                                                                                                                
        if var == 'ptposptm':
            histos[rec].Fill(lp.Pt() + lm.Pt(),weight)                                                                                                                                                                          
            if not isData:
                histos[gen].Fill(glp.Pt() + glm.Pt(),weight)

    # Close file
    fIn.Close()
 
"""
Wrapper for when the analysis is run in parallel
"""
def createSummaryPacked(args):
    filename,isData,outDir,filen,var = args
    try:
        # Define histograms                                                                  
        histos = getAnalysisHistograms(var)

        # Filling histograms with new binning                                           
        createHistos(filen,isData,histos,var)

        # Dump histograms in a file                                                                        
        fOut=ROOT.TFile.Open(os.path.join(outDir,filename),'RECREATE')
        for h in histos: histos[h].Write()
        fOut.Close()

    except ReferenceError:
        print 50*'<'
        print "  Problem with", filename, "continuing without"
        print 50*'<'

"""
Wrapper for mass files
"""        
def createSummaryTasksMass(opt):
    
    tasklist = []

    output = os.path.join(opt.outDir, opt.var)

    # Create list of massfiles
    masstrees, massfiles = getMassTrees(opt.input, verbose=True)
    masspoints = sorted(list(set([mass for mass,_ in masstrees.keys()])))

    mass = [166, 169, 171, 173, 175, 178]

    # Loop over mass files to process
    for m in mass:
        m = str(m)
        mass_out = 'mass_scan/' + m + '/'
        out = os.path.join(output, mass_out)
        os.system('mkdir -p %s' %out)
        
        var = opt.var 
        for filename in os.listdir(opt.input):
            if not os.path.splitext(filename)[1] == '.root': continue
            isData = True if 'Data' in filename else False
            if not m in filename: continue
            filen = os.path.join(opt.input,filename)

            # Select type: Process only TTJets not SingleTop, etc..

            if opt.type == 'mass_TTJets':
                if not 'TTJets_MSDecays_' in filename: continue
                print 'processing file%s'%filename
                tasklist.append((filename,isData,out,filen,var))

            if opt.type == 'mass_all':
                print 'processing file %s'%filename
                tasklist.append((filename,isData,out,filen,var))
        
        print ' Submitting jobs in %d threads for %s and %s' % (opt.jobs,m,var)
        import multiprocessing as MP
        pool = MP.Pool(opt.jobs)
        pool.map(createSummaryPacked,tasklist)

"""
Wrapper for files
"""
def createSummaryTasks(opt):

    tasklist=[]

    output = os.path.join(opt.outDir, opt.var)

    # Loop under files stored in eos
    if opt.input.find('/store')>=0:
        print 'going through files'
        for filename in os.listdir(opt.input):
            if not os.path.splitext(filename)[1] == '.root': continue   
            isData = True if 'Data' in filename else False
            filen = os.path.join(opt.input,filename)
            if opt.type == 'all':
                tasklist.append((filename,isData,output,filen,opt.var))
            if opt.type == 'TTJets':
                if not 'MC8TeV_TTJets_MSDecays_172v5' in filename: continue
                tasklist.append((filename,isData,output,filen,opt.var))
    else:
        for filename in os.listdir(args[0]):
            if not os.path.splitext(filename)[1] == '.root': continue   
            isData = True if 'Data' in filename else False
            filen = os.path.join(opt.input,filename)
            tasklist.append((filename,isData,output,filen,opt.var))

    # Loop over tasks
    if opt.jobs>0:
        print ' Submitting jobs in %d threads' % opt.jobs
        import multiprocessing as MP
        pool = MP.Pool(opt.jobs)
        pool.map(createSummaryPacked,tasklist)
    else:
        for filename,isData,outDir in tasklist:
            createSummaryPacked(filename=filename,isData=isData,outDir=output,var=opt.var)

	return 0

"""
steer
"""
def main():
	usage = 'usage: %prog [options]'
	parser = optparse.OptionParser(usage)
	parser.add_option('-i', '--input',
                          dest='input',   
                          default='/eos/uscms/store/user/cmantill/Dileptons2012/ttbar/', 
                          help='input directory with the files [default: %default]')
        parser.add_option('-o', '--outDir',
                          dest='outDir',
                          default='distributions',
                          help='Output directory [default: %default]')
	parser.add_option('-v', '--var',
                          dest='var', 
                          default='ptll',
                          help='Variable. Choose between: ptpos, Epos, ptll, mll, EposEm, ptposptm [default: %default]')
	parser.add_option('--jobs',
                          dest='jobs', 
                          default=10,
                          type=int,
                          help='# of jobs to process in parallel the trees [default: %default]')
        parser.add_option('-t',
                          dest='type', 
                          default='all',
                          help='Select if you want to process all the files or just the nominal samples. Options: all, TTJets, mass_TTJets, mass_all')
        parser.add_option('-l', '--lumi', dest='lumi', default=19701,
                          type='float',
                          help='Re-scale to integrated luminosity [pb]'
                          ' [default: %default]')
        parser.add_option('--normToData', dest='normToData', action="store_true",
                          help='Force normalization to data')
        parser.add_option('-s', '--silent', dest='silent', action="store_true",
                          help='Silent mode (no plots) [default: %default]')
        parser.add_option('-d', '--debug', dest='debug', action="store_true",
                          help='Dump the event yields table for each plot')


	(opt, args) = parser.parse_args()

	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetOptTitle(0)
	ROOT.gROOT.SetBatch(True)
	setTDRStyle()
        ROOT.gSystem.Load("libUserCodeTopMassSecVtx")
        ROOT.AutoLibraryLoader.enable()

        output = os.path.join(opt.outDir, opt.var)
	os.system('mkdir -p %s' % output)
       
	if 'mass' in opt.type:
                print 80*'-'
                print 'Creating ROOT file with data and background distributions of %s from %s'%(opt.var,opt.input)
                createSummaryTasksMass(opt)
                print 80*'-'
	else:
		print 80*'-'
               	print 'Creating ROOT file with data and background distributions of %s from %s'%(opt.var,opt.input)
         	createSummaryTasks(opt)
                
                print 80*'='
                print ' Producing plots to unfold for %s'%(opt.var)
                toUnfolddir = os.path.join(output, 'toUnfold')
                os.system('mkdir -p %s'% toUnfolddir)
                tUnfoptions = deepcopy(opt)
                tUnfoptions.outDir = toUnfolddir
                tUnfoptions.json = 'test/topss2014/samples.json'
                tUnfoptions.outFile = 'toUnfold.root' 
                tUnfoptions.filter = ""
                tUnfoptions.split = False
                tUnfoptions.verbose = 0
                tUnfoptions.plotMask = ''
                tUnfoptions.ratioRange = '0.5,1.5'
                tUnfoptions.excludeProcesses = ''
                tUnfoptions.cutUnderOverFlow = True

                runPlotter(args[0], tUnfoptions)                
        	print 80*'-'

        return 0

if __name__ == "__main__":
	sys.exit(main())
