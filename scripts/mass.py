#!/usr/bin/env python
from array import array
import ROOT
from ROOT import *
import optparse
import os,sys,numpy
import pprint

from runDileptonUnfolding import getPathToObjects
from UserCode.TopMassSecVtx.storeTools_cff import fillFromStore
from UserCode.TopMassSecVtx.PlotUtils import RatioPlot

"""
Get purity,stability and efficiency
"""
def getStability(root,v,m):

    histos={}
    
    #open ROOT file
    fIn=ROOT.TFile.Open(root)

    print 'analizing %s'%root

    var = str(v)

    genList       = getPathToObjects(fIn.Get('%s_gen'%(var)))
    for p in genList:
        h=fIn.Get(p)
        hname=h.GetName()
        # Signal + bkg
        if not 'gen' in histos : 
            histos['gen']=h.Clone('gen')
        else : 
            histos['gen'].Add(h)

    recList       = getPathToObjects(fIn.Get('%s_rec'%(var)))
    for p in recList:
        h=fIn.Get(p)
        hname=h.GetName()
        # Signal + bkg
        if not 'rec' in histos : 
            histos['rec']=h.Clone('rec')
        else : 
            histos['rec'].Add(h)

    migrationList = getPathToObjects(fIn.Get('%s_migration'%(var)))
    for p in migrationList:
        h=fIn.Get(p)
        hname=h.GetName()
        # Signal + bkg
        if not 'mig' in histos : 
            histos['mig']=h.Clone('mig')
        else :
            histos['mig'].Add(h)

    #pT positive lepton
    if var == 'ptpos':            
        bins_gen = [22.72,25.45,28.18,30.92,33.67,36.42,39.17,42.01,44.89,47.77,50.94,55.06,59.19,63.63,68.15,74.76,83.59,95,115.65]
    #pT charged-lepton pair
    if var == 'ptll': 
        bins_gen = [15.13,23.05,30,36.06,41.76,46.8,51.57,55.88,60.23,65.23,70.22,75.07,79.93,84.58,89.23,95.79,104.48,116.54,136.92]
    #M charged-lepton pair
    if var == 'mll': 
        bins_gen = [30.85,40.55,48.96,55.78,62.37,68.81,75.26,81.75,88.33,95.63,102.76,109.62,117.89,127.15,136.41,148.52,167.91,193.4,235.33]
    #Scalar sum of E
    if var == 'EposEm': 
        bins_gen = [75.17,86.46,95.35,104.8,114.58,122.55,130.11,137.5,146.95,156.88,166.87,176.91,188.33,202.62,223.57,244.76,271.12,307.35,364.29]
    #Scalar sum of Pt        
    if var == 'ptposptm': 
        bins_gen = [55.73,63.41,68.94,74.8,80.52,84.76,89.01,93.32,97.64,102.22,107.1,112.53,118.8,126.65,134.43,142.8,154.04,170.92,204.56]
    
    title_p = ';p_{T}(l^{+}) [GeV];Purity'
    title_s = ';p_{T}(l^{+}) [GeV];Stability'
    title_e = ';p_{T}(l^{+}) [GeV];Efficiency'

    nbins_gen = histos['gen'].GetXaxis().GetNbins()
    nbins_rec = histos['mig'].GetYaxis().GetNbins()

    histos['Purity']=ROOT.TH1F('Purity',title_p,len(bins_gen)-1,array('d',bins_gen))
    histos['Stability']=ROOT.TH1F('Stability',title_s,len(bins_gen)-1,array('d',bins_gen))
    histos['Efficiency']=ROOT.TH1F('Efficiency',title_e,len(bins_gen)-1,array('d',bins_gen))

    for i in xrange(1,len(bins_gen)):

        # Purity: # of events generated and correctly reconstructed in a given bin i relative to the # of events that are reconstructed in bin i but generated anywhere
        Nrec_gen_i = histos['mig'].ProjectionY("test",i,i+1,"e").Integral(2*i,2*i+1)
  
        Nrec_i = histos['rec'].GetBinContent(2*i+1) + histos['rec'].GetBinContent(2*i)
        pur = Nrec_gen_i/Nrec_i
        histos['Purity'].SetBinContent(i,pur)

        Ngen_i = histos['gen'].GetBinContent(i)
        Ngen_rec_i = histos['mig'].ProjectionX("test2",2*i,2*i+1,"e").Integral()

        stab = Nrec_gen_i/Ngen_i
        histos['Stability'].SetBinContent(i,stab)

        #eff = Nrec_gen_i/histos['rec'].Integral()
        histos['Efficiency'].SetBinContent(i, Nrec_gen_i)
         
    scale = 1.0/ histos['Efficiency'].Integral()
    histos['Efficiency'].Scale(scale)

    # Convert histograms to TGraphErrors
    graphs ={}

    graphs['Purity'] = ROOT.TGraphErrors(histos['Purity'])
    graphs['Efficiency'] = ROOT.TGraphErrors(histos['Efficiency'])
    graphs['Stability'] = ROOT.TGraphErrors(histos['Stability'])
   
    # Save histograms in canvas
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)

    #Create a canvas for plotting your graph
    c = ROOT.TCanvas('c','c')
    gStyle.SetCanvasDefH(600);
    gStyle.SetCanvasDefW(800);
    c.SetLeftMargin(0.15);
    c.SetRightMargin(0.25)
    c.SetBottomMargin(0.25);
    p1=ROOT.TPad('p1','p1',0.,0.,1.0,1.0)
    p1.Draw()

    #Create CMS labels
    tlat1 = TLatex()
    tlat1.SetNDC()
    tlat1.SetTextFont(42)
    tlat1.SetTextSize(0.04)
    tlat1.SetTextAlign(31)
    tlat1.DrawLatex(0.16, 0.91,'#bf{CMS}')

    tlat2 = TLatex()
    tlat2.SetNDC()
    tlat2.SetTextFont(42)
    tlat2.SetTextSize(0.04)
    tlat2.SetTextAlign(31)
    tlat2.DrawLatex(0.30, 0.91, '#it{Simulation}')
    tlat2.DrawLatex(0.90, 0.91,'8 TeV')

    tleg = ROOT.TLegend(0.70, 0.75, .89, 0.89)
    tleg.SetBorderSize(0)
    tleg.SetFillColor(0)
    tleg.SetFillStyle(0)
    tleg.SetShadowColor(0)
    tleg.SetTextFont(43)
    tleg.SetTextSize(20)

    #pT positive lepton
    if var == 'ptpos': title = 'p_{T}(l^{+}) [GeV]'
    #pT charged-lepton pair
    if var == 'ptll': title = 'p_{T}(l^{+}l^{-}) [GeV]'
    #M charged-lepton pair
    if var == 'mll': title = 'M(l^{+}l^{-}) [GeV]'
    #Scalar sum of E
    if var == 'EposEm': title = 'E(l^{+})+E(l^{-}) [GeV]'
    #Scalar sum of Pt
    if var == 'ptposptm': title = 'p_{T}(l^{+})+p_{T}(l^{-}) [GeV]'

    p1.cd()
    graphs['Purity'].SetTitle("")
    graphs['Purity'].SetLineColor(kBlue-3)
    graphs['Purity'].SetLineStyle(1)
    graphs['Purity'].SetLineWidth(1)
    graphs['Purity'].SetMarkerStyle(8)
    graphs['Purity'].SetMarkerColor(kBlue-3)
    graphs['Purity'].SetMarkerSize(0.9)
    graphs['Purity'].SetMaximum(1.0)
    graphs['Purity'].SetMinimum(0.0)
    graphs['Purity'].GetXaxis().SetLabelSize(0.037)
    graphs['Purity'].GetYaxis().SetTitleOffset(1.2)
    graphs['Purity'].GetYaxis().SetTitleSize(0.037)
    graphs['Purity'].GetYaxis().SetLabelSize(0.037)
    graphs['Purity'].GetXaxis().SetTitle(title)
    graphs['Purity'].GetYaxis().SetTitle()
    for i in xrange(0,len(bins_gen)-1):
      errx = graphs['Purity'].GetErrorX(i)
      graphs['Purity'].SetPointError(i,errx,0)
    graphs['Purity'].Draw("Ap")
    tleg.AddEntry(graphs['Purity'], 'Purity', 'pl')

    graphs['Efficiency'].SetMarkerStyle(23)
    graphs['Efficiency'].SetMarkerColor(kGreen-3)
    graphs['Efficiency'].SetMarkerSize(0.9)
    graphs['Efficiency'].SetLineColor(kGreen-3)    
    graphs['Efficiency'].SetLineWidth(1)
    for i in xrange(0,len(bins_gen)-1):
      errx = graphs['Efficiency'].GetErrorX(i)
      graphs['Efficiency'].SetPointError(i,errx,0)
    graphs['Efficiency'].Draw("Samep")
    tleg.AddEntry(graphs['Efficiency'], 'Efficiency', 'pl')

    graphs['Stability'].SetMarkerStyle(21)
    graphs['Stability'].SetMarkerColor(kRed)
    graphs['Stability'].SetMarkerSize(0.9)
    graphs['Stability'].SetLineColor(kRed)    
    graphs['Stability'].SetLineWidth(1)
    for i in xrange(0,len(bins_gen)-1):
      errx = graphs['Stability'].GetErrorX(i)
      graphs['Stability'].SetPointError(i,errx,0)
    graphs['Stability'].Draw("Samep")
    tleg.AddEntry(graphs['Stability'], 'Stability', 'pl')

    tleg.Draw()

    outputName = 'unfoldResults/'+var+'/pur_stab_eff'
    for ext in ['png','pdf']:
      c.SaveAs('%s.%s'%(outputName,ext))
    del c

"""
Get first three Mellin moments from each of the distributions
"""
def getMoments(root,v,moments,m):

    histos={}
    
    #open ROOT file
    fIn=ROOT.TFile.Open(root)

    var = str(v)

    unfolded       = getPathToObjects(fIn.Get('%s_unfolded_wgt'%var))
    for p in unfolded:
      h = fIn.Get(p)
      hname = h.GetName()
      if moments == 'unf': 
        if 'data_unfolded' in hname:
          histos['unf'] =  h.Clone('data_unfolded')
      # if moments == 'gen': 
      #   if '%s_gen'%var in hname:
      #     histos['gen'] = h.Clone('%s_gen'%var)
      if moments == 'gen': 
        if 'signal_gen' in hname:
          histos['gen'] = h.Clone('signal_gen')

    #detach objects from file and close file
    for h in histos: histos[h].SetDirectory(0)
    fIn.Close()

    mean = histos[moments].GetMean()
    meanerr = histos[moments].GetMeanError()
    rms = histos[moments].GetRMS()
    rmserr = histos[moments].GetRMSError()

    u1 = mean
    u1err = meanerr
    u2 = mean**2 + rms**2
    u2err = sqrt((2*rms*rmserr)**2+(2*mean*meanerr)**2)

    #Create a text file with the three different moments for the corresponding mass
    outfileName = 'Results/'+var+'/moments_' + var + '_' + moments +'nobkg.txt'
    m = float(m)
    mass = m + 0.5
    outfile = open(outfileName,'a')
    outfile.write('%4.6f %4.6f %4.6f %4.6f %4.6f\n' %(mass,u1,u1err,u2,u2err))
    outfile.close()

from UserCode.TopMassSecVtx.PlotUtils import *

"""
Make template fit
"""
def getFit(root,v):

    histos={}
    
    #open ROOT file
    fIn=ROOT.TFile.Open(root)

    var = str(v)

    unfolded       = getPathToObjects(fIn.Get('%s_unfolded_wgt'%var))
    for p in unfolded:
      h = fIn.Get(p)
      hname = h.GetName()
      if moments == 'unf': 
        if 'data_unfolded' in hname:
          histos['unf'] =  h.Clone('data_unfolded')
      # if moments == 'gen': 
      #   if '%s_gen'%var in hname:
      #     histos['gen'] = h.Clone('%s_gen'%var)
      if moments == 'gen': 
        if 'signal_gen' in hname:
          histos['gen'] = h.Clone('signal_gen')
          
    #detach objects from file and close file
    for h in histos: histos[h].SetDirectory(0)
    fIn.Close()

    #Create a canvas for plotting your graph
    c = ROOT.TCanvas('c','c')
    gStyle.SetCanvasDefH(600);
    gStyle.SetCanvasDefW(800);
    c.SetLeftMargin(0.15);
    c.SetRightMargin(0.25)
    c.SetBottomMargin(0.25);
    pad1=ROOT.TPad('p1','p1',0.,0.,1.0,1.0)
    pad1.Draw()

    #Create CMS labels
    tlat1 = TLatex()
    tlat1.SetNDC()
    tlat1.SetTextFont(42)
    tlat1.SetTextSize(0.04)
    tlat1.SetTextAlign(31)
    tlat1.DrawLatex(0.16, 0.91,'#bf{CMS}')

    tlat2 = TLatex()
    tlat2.SetNDC()
    tlat2.SetTextFont(42)
    tlat2.SetTextSize(0.04)
    tlat2.SetTextAlign(31)
    tlat2.DrawLatex(0.30, 0.91, '#it{Simulation}')
    tlat2.DrawLatex(0.90, 0.91,'8 TeV')

    tleg = ROOT.TLegend(0.12, 0.75, .50, 0.89)
    tleg.SetBorderSize(0)
    tleg.SetFillColor(0)
    tleg.SetFillStyle(0)
    tleg.SetShadowColor(0)
    tleg.SetTextFont(43)
    tleg.SetTextSize(13)

    #w.factory("Gaussian::gauss(mes[5.20,5.30],mean[5.28,5.2,5.3],width[0.0027,0.001,1])");

    #Define the fit function and make the fit 
    tf = ROOT.TF1("f1", "[0]*exp(-0.5*((x-[1])/(([2]+[3])+([3]-[2])*tanh())**2)", 160., 190.)
    tf.SetLineColor(kGreen+4)
    tf.SetLineWidth(2)
    tf.SetLineStyle(1)
    graphs['u1_unf'].Fit(tf, "WQ")

    tf = histos['unf'].GetFunction("tf")
    p0_unf=linear.GetParameter(0)
    p0Err_unf=linear.GetParError(0)
    p1_unf=linear.GetParameter(1)
    p1Err_unf=linear.GetParError(1)
    chi2_unf=linear.GetChisquare()

"""
Compare distributions from different masses
"""
def getRatioMass(compare,v):

    var = str(v)
      
    #Open root file containing histograms
    masspoints = [166,169,171,173,175,178]

    #outDir = 'results/%s'%var
    outDir = 'unfoldResults/%s'%var

    histos = {}

    if compare == 'unf':
        filename_nominal = 'unfoldResults/%s/nominal/plots/%s_unfolded_wgt.root'%(var,var)
        fIn_n=ROOT.TFile.Open(filename_nominal)
        unfList_nominal  = getPathToObjects(fIn_n.Get('%s_unfolded_wgt'%(var)))
        for p in unfList_nominal:
         h=fIn_n.Get(p)
         hname=h.GetName()
         if 'data_unfolded' in hname : 
             histos['172'] = h.Clone('data_unfolded')   

    if compare == 'gen' or compare == 'rec':
        filename_nominal = 'unfoldResults/%s/nominal/plots/plotter.root'%(var)
        fIn_n=ROOT.TFile.Open(filename_nominal)
        recList_nominal  = getPathToObjects(fIn_n.Get('%s_rec_wgt'%(var)))
        genList_nominal  = getPathToObjects(fIn_n.Get('%s_gen_wgt'%(var)))

    if compare == 'gen':
      for p in genList_nominal:
          h=fIn_n.Get(p)
          hname=h.GetName()
          if 'TTJets' in hname or 'SingleT' in hname:
            if not '172' in histos : 
                histos['172']=h.Clone('172')
            else : 
               histos['172'].Add(h)
    
    if compare == 'rec':
      for p in recList_nominal:
          h=fIn_n.Get(p)
          hname=h.GetName()
          if 'TTJets' in hname or 'SingleT' in hname:
            if not '172' in histos : 
                histos['172']=h.Clone('172')
            else : 
               histos['172'].Add(h)

    for m in masspoints:

        m = str(m)

        if compare == 'unf':
            filename = 'unfoldResults/%s/mass_scan/%s/plots/%s_unfolded_wgt.root'%(var,m,var)
            fIn=ROOT.TFile.Open(filename)
            unfList  = getPathToObjects(fIn.Get('%s_unfolded_wgt'%(var)))
            for p in unfList:
             h=fIn.Get(p)
             hname=h.GetName()
             if 'data_unfolded' in hname : 
                 histos[m] = h.Clone('data_unfolded')   

        if compare == 'gen' or compare == 'rec':
            filename = 'unfoldResults/%s/mass_scan/%s/plots/plotter.root'%(var,m)
            fIn=ROOT.TFile.Open(filename)
            recList       = getPathToObjects(fIn.Get('%s_rec_wgt'%(var)))
            genList       = getPathToObjects(fIn.Get('%s_gen_wgt'%(var)))

        if compare == 'gen':
          for p in genList:
            h=fIn.Get(p)
            hname=h.GetName()
            if 'TTJets' in hname:
              if not m in histos : 
                  histos[m]=h.Clone(m)
              else : 
                 histos[m].Add(h)

        if compare == 'rec':
          for p in recList:
            h=fIn.Get(p)
            hname=h.GetName()
            if 'TTJets' in hname:
              if not m in histos : 
                  histos[m]=h.Clone(m)
              else : 
                 histos[m].Add(h)

    ratplot = RatioPlot('ratioplot')
    ratplot.ratiotitle = "Ratio wrt 172.5 GeV"
    ratplot.ratiorange = (0.8, 1.2)

    reference = [histos['172']]
    ratplot.reference = reference

    colors = [ROOT.kViolet-6,
              ROOT.kBlue-3, 
              ROOT.kRed-4, 
              ROOT.kOrange-3,
              ROOT.kGreen-3,
              ROOT.kGreen+2,
              ROOT.kAzure-2,]
    
    ratplot.colors1=colors
    ratplot.legpos = (0.68, 0.40)
    ratplot.tagpos = (0.88,0.87)

    for h in sorted(histos.keys()):
        legentry = 'm_{t} = %s.5 GeV' % h
        histos[h].SetFillColor(0)
        try:
          histo = histos[h]
          ratplot.add(histo, legentry)
        except KeyError: pass
   
    if compare == 'unf':
      ratplot.tag = 'Unfolded'
      ratplot.show("%s_unf_comparisonNorm"%var, outDir)
    if compare == 'gen':
      ratplot.tag = 'Generated level'
      ratplot.show("%s_gen_comparisonNorm"%var, outDir)
    if compare == 'rec':
      ratplot.tag = 'Reconstructed level'
      ratplot.show("%s_rec_comparisonNorm"%var, outDir)
    
"""
Graph Mellin Moments
"""
def graphMoments(v):

     var = str(v)

     filename_unf = 'Results/'+var+'/moments_'+var+'_unfnobkg.txt'
     filename_gen = 'Results/'+var+'/moments_'+var+'_gennobkg.txt'

     #pT positive lepton
     if var == 'ptpos': title = 'p_{T}(l^{+})'
     #pT charged-lepton pair
     if var == 'ptll': title = 'p_{T}(l^{+}l^{-})'
     #M charged-lepton pair
     if var == 'mll': title = 'M(l^{+}l^{-})'       
     #Scalar sum of E
     if var == 'EposEm': title = 'E(l^{+})+E(l^{-})'
     #Scalar sum of Pt
     if var == 'ptposptm': title = 'p_{T}(l^{+})+p_{T}(l^{-})'

     m_unf,u1_unf,erru1_unf,u2_unf,erru2_unf = numpy.loadtxt(filename_unf,unpack =True)
     m_gen,u1_gen,erru1_gen,u2_gen,erru2_gen = numpy.loadtxt(filename_gen,unpack =True)

     mass_gen= array('d', m_gen)
     e0_gen = array('d', [0.0]*len(mass_gen))
     u1_gen = array('d', u1_gen)
     erru1_gen = array('d', erru1_gen)
     u2_gen = array('d', u2_gen)
     erru2_gen = array('d', erru2_gen)

     mass_unf= array('d', m_unf)
     e0_unf = array('d', [0.0]*len(mass_unf))
     u1_unf = array('d', u1_unf)
     erru1_unf = array('d', erru1_unf)
     u2_unf = array('d', u2_unf)
     erru2_unf = array('d', erru2_unf)

     #Create a canvas for plotting your graph
     c = ROOT.TCanvas('c','c')
     gStyle.SetCanvasDefH(600);
     gStyle.SetCanvasDefW(800);
     c.SetLeftMargin(0.15);
     c.SetRightMargin(0.25)
     c.SetBottomMargin(0.25);
     pad1=ROOT.TPad('p1','p1',0.,0.,1.0,1.0)
     pad1.Draw()

     #Create CMS labels
     tlat1 = TLatex()
     tlat1.SetNDC()
     tlat1.SetTextFont(42)
     tlat1.SetTextSize(0.04)
     tlat1.SetTextAlign(31)
     tlat1.DrawLatex(0.16, 0.91,'#bf{CMS}')

     tlat2 = TLatex()
     tlat2.SetNDC()
     tlat2.SetTextFont(42)
     tlat2.SetTextSize(0.04)
     tlat2.SetTextAlign(31)
     tlat2.DrawLatex(0.30, 0.91, '#it{Simulation}')
     tlat2.DrawLatex(0.90, 0.91,'8 TeV')

     tleg = ROOT.TLegend(0.12, 0.75, .50, 0.89)
     tleg.SetBorderSize(0)
     tleg.SetFillColor(0)
     tleg.SetFillStyle(0)
     tleg.SetShadowColor(0)
     tleg.SetTextFont(43)
     tleg.SetTextSize(13)

     graphs = {}
     multigraphs = {}

     graphs['u1_unf'] = ROOT.TGraphErrors(len(mass_unf),mass_unf,u1_unf,e0_unf,erru1_unf)
     graphs['u1_gen'] = ROOT.TGraphErrors(len(mass_gen),mass_gen,u1_gen,e0_gen,erru1_gen)

     graphs['u2_unf'] = ROOT.TGraphErrors(len(mass_unf),mass_unf,u2_unf,e0_unf,erru2_unf)
     graphs['u2_gen'] = ROOT.TGraphErrors(len(mass_gen),mass_gen,u2_gen,e0_gen,erru2_gen)

     multigraphs['u1'] = ROOT.TMultiGraph() # Comparison between unfolded and generated level distribution
     multigraphs['u2'] = ROOT.TMultiGraph() # Comparison between unfolded and generated level distribution

     colors = [kBlue+1,kViolet-6,kRed,kBlue-5,kGreen+3,kRed+3]
     markers = [8,21,22,23,24,25,26]
     lines = [2,3,4,5,6,7,1,6,5]
     i=0

     tlat3 = TLatex()
     tlat3.SetNDC()
     tlat3.SetTextFont(42)
     tlat3.SetTextSize(0.03)
     tlat3.SetTextAlign(31)
     tlat3.DrawLatex(0.8,0.24,'Fit Results')

     gStyle.SetOptStat(0)
     gStyle.SetOptTitle(0)

     graphs['u1_unf'].SetMarkerColor(kBlue+1)
     graphs['u1_unf'].SetMarkerStyle(22)
     graphs['u1_unf'].SetMarkerSize(0.9)
     graphs['u1_unf'].SetLineColor(kBlue+1)
     tlat3.SetTextColor(kBlue+1)
     tleg.AddEntry(graphs['u1_unf'], '#mu^{(1)}_{%s} Unfolded'%title, 'pl')

     #Define the fit function and make the linear fit 
     tf = ROOT.TF1("linear", "pol1", 160., 190.)
     tf.SetLineColor(kGreen+4)
     tf.SetLineWidth(2)
     tf.SetLineStyle(1)
     graphs['u1_unf'].Fit(tf, "WQ")

     tf = graphs['u1_unf'].GetFunction("tf")
     p0_unf=linear.GetParameter(0)
     p0Err_unf=linear.GetParError(0)
     p1_unf=linear.GetParameter(1)
     p1Err_unf=linear.GetParError(1)
     chi2_unf=linear.GetChisquare()

     outfileName = 'Results/%s/p1_unf.txt'%var
     outfile = open(outfileName,'a')
     outfile.write('%4.6f\n' %(p1_unf))
     outfile.close()
     print p0_unf 
     print p0Err_unf 
     print p1_unf 
     print p1Err_unf
     print chi2_unf
     tlat3.DrawLatex(0.8,0.19,'#mu^{(1)}_{%s}=%3.4f #pm %3.4f + (%3.4f #pm %3.4f)*m_{t}'%(title,p0_unf,p0Err_unf,p1_unf,p1Err_unf))
     multigraphs['u1'].Add(graphs['u1_unf'])

     graphs['u1_gen'].SetMarkerColor(kGreen-3)
     graphs['u1_gen'].SetMarkerStyle(22)
     graphs['u1_gen'].SetMarkerSize(0.9)
     graphs['u1_gen'].SetLineColor(kGreen-3)
     tlat3.SetTextColor(kGreen-3)
     tleg.AddEntry(graphs['u1_gen'], '#mu^{(1)}_{%s} Generated'%title, 'pl')

     tf1 = ROOT.TF1("linear_g", "pol1", 160., 190.)
     tf1.SetLineColor(kGreen-3)
     tf1.SetLineWidth(2)
     tf1.SetLineStyle(2)
     graphs['u1_gen'].Fit(tf1, "WQ")

     p0_gen=linear_g.GetParameter(0)
     p0Err_gen=linear_g.GetParError(0)
     p1_gen=linear_g.GetParameter(1)
     p1Err_gen=linear_g.GetParError(1)
     chi2=linear_g.GetChisquare()
     outfileName = 'Results/%s/p1_gen.txt'%var
     outfile = open(outfileName,'a')
     outfile.write('%4.6f\n' %(p1_gen))
     outfile.close()
     tlat3.DrawLatex(0.8,0.14,'#mu^{(1)}_{%s}=%3.4f #pm %3.4f + (%3.4f #pm %3.4f)*m_{t}'%(title,p0_gen,p0Err_gen,p1_gen,p1Err_gen))
     multigraphs['u1'].Add(graphs['u1_gen'])

     # th = ROOT.TF1("th", "0.1347*173+0.1939*x", 160., 190.)
     # th.SetLineColor(kRed-4)
     # th.SetLineWidth(4)
     # th.SetLineStyle(2)
     # tleg.AddEntry(th, '#mu^{(1)}_{%s} TH'%title, 'pl')

     pad1.cd()
     multigraphs['u1'].Draw("AP")
     multigraphs['u1'].GetXaxis().SetTitle("Top Mass [GeV]")
     multigraphs['u1'].GetYaxis().SetTitle("#mu^{(1)}_{%s} [GeV]"%title)
     multigraphs['u1'].GetXaxis().SetTitleSize(0.037)
     multigraphs['u1'].GetXaxis().SetLabelSize(0.037)
     multigraphs['u1'].GetYaxis().SetTitleOffset(1.2)
     multigraphs['u1'].GetYaxis().SetTitleSize(0.037)
     multigraphs['u1'].GetYaxis().SetLabelSize(0.037)
        
     tleg.Draw()

     outputName = 'Results/'+var+'/u1_'+var+'nobkg'
     for ext in ['png','pdf']:
      c.SaveAs('%s.%s'%(outputName,ext))
     del c

def plotter(opt):

    var = ['ptpos']
    #var = ['ptll','mll','EposEm','ptposptm']

    #Tau Test
    for v in var: os.system("python scripts/runDileptonUnfolding.py -r Results/%s/nominal/plots/plotter.root -v %s -o Results/%s/tautest -g f -w b"%(v,v,v)) 

    # for v in var:

    #     #pT positive lepton
    #     if v == 'ptpos': title = 'p_{T}(l^{+})'
  
     #    #pT charged-lepton pair
     #    if v == 'ptll': title = 'p_{T}(l^{+}l^{-})'
     #    #M charged-lepton pair
     #    if v == 'mll': title = 'M(l^{+}l^{-})'
     #    #Scalar sum of E
     #    if v == 'EposEm': title = 'E(l^{+})+E(l^{-})'
     #    #Scalar sum of Pt
     #    if v == 'ptposptm': title = 'p_{T}(l^{+})+p_{T}(l^{-})'

        #print 'doing something'
        #tau = [8e-6,1e-5,4e-5,8e-5,1e-4,2e-4,4e-4,8e-4,1e-3,2e-3,3e-3,4e-3,5e-3,6e-3,7e-3,8e-3,9e-3,1e-2,2e-2]
        #tau = numpy.arange(0.484,1,1e-2)
        #for i in tau:
        # print 'doing something'
        # print i
        # os.system("python scripts/runDileptonUnfolding.py -r Results/%s/nominal/plots/plotter.root -v %s -o Results/%s/tautest -g f -t %s"%(v,v,v,i))
  
        # c = ROOT.TCanvas('c','c')
        # gStyle.SetCanvasDefH(600);
        # gStyle.SetCanvasDefW(800);
        # c.SetLeftMargin(0.15);
        # c.SetRightMargin(0.25)
        # c.SetBottomMargin(0.25);
        # pad1=ROOT.TPad('p1','p1',0.,0.,1.0,1.0)
        # pad1.Draw()

        # filename = 'Results/%s/tau_%s.txt'%(v,v)

        # tau,rho = numpy.loadtxt(filename,unpack =True)
        # tau = array('d', tau)
        # rho = array('d', rho)
        # taugraph = ROOT.TGraph(len(tau),tau,rho)
        # taugraph.SetMarkerColor(kBlue+1)
        # taugraph.SetMarkerStyle(22)
        # taugraph.SetMarkerSize(0.4)
        # taugraph.SetLineColor(kBlue+1)

        # #taugraph.GetXaxis().SetRangeUser(0.0001,0.02)
        # # pol2 = ROOT.TF1("pol2", "pol2", 0.006, 0.014)
        # # pol2.SetLineColor(kGreen-3)
        # # pol2.SetLineWidth(2)
        # # pol2.SetLineStyle(2)
        # # taugraph.Fit("pol2","R")
        # # min = pol2.GetMinimumX(0.006,0.014)
        # # print min

        # pad1.cd()
        # taugraph.Draw("AP")

        # taugraph.GetXaxis().SetTitle("#tau")
        # taugraph.GetYaxis().SetTitle("#rho^{avg}_{%s}"%title)
        # taugraph.GetXaxis().SetTitleSize(0.037)
        # taugraph.GetXaxis().SetLabelSize(0.037)
        # taugraph.GetYaxis().SetTitleOffset(1.2)
        # taugraph.GetYaxis().SetTitleSize(0.037)
        # taugraph.GetYaxis().SetLabelSize(0.037)
        # #taugraph.GetXaxis().SetLimits(0.006,0.014)

        # outputName = 'Results/'+v+'/tau' 
        # for ext in ['png','pdf']:
        #  c.SaveAs('%s.%s'%(outputName,ext))
        # del c

      #os.system("python scripts/runDileptonUnfolding.py -n /store/cmst3/group/top/summer2015/treedir_bbbcb36/ttbar/ -o Results/%s/nominal/ --jobs 8 -v %s"%(v,v))
      #os.system("python scripts/runDileptonUnfolding.py -n templates/nominal/ -o Results/%s/nominal/ --jobs 8 -v %s"%(v,v))
      #os.system("python scripts/runPlotter.py -j test/topss2014/samples.json -o Results/%s/nominal/plots Results/%s/nominal/ --cutUnderOverFlow"%(v,v))
      #os.system("python scripts/runDileptonUnfolding.py -r Results/%s/nominal/MC8TeV_TTJets_MSDecays_172v5.root -v %s -o Results/%s/nominal/wgt/ -g store -w b"%(v,v,v)) 
      #os.system("python scripts/runDileptonUnfolding.py -r Results/%s/nominal/plots/plotter.root -v %s -o Results/%s/nominal/ -g store -w b "%(v,v,v)) 

      #m = 172
      #root = 'Results/%s/nominal/plots/%s_unfolded_wgt.root'%(v,v)
      #root =  'unfoldResults/%s/nominal/plots/%s_unfolded.root'%(v,v)
      #getMoments(root,v,'unf',m)
      #getMoments(root,v,'gen',m)
      #root = 'unfoldResults/%s/nominal/plots/plotter.root'%(v)
      #getMoments(root,v,'rec',m)

      # root = 'unfoldResults/%s/nominal/plots/plotter.root'%(v)
      # getStability(root,v,m)

      #os.system("python scripts/runDileptonUnfolding.py -i /store/cmst3/group/top/summer2015/treedir_bbbcb36/ttbar/mass_scan -o Results/%s/mass_scan/ --jobs 8 -m True -v %s"%(v,v))

      # Create an array with the different masses
      #mass = [166,169,171,173,175,178]
      #mass = [169]
                 
      #for m in mass:
      #  m = str(m)

      #   # Get all the plots
        #os.system("python scripts/runPlotter.py -j test/topss2014/mass_scan_samples.json -o Results/%s/mass_scan/%s/plots Results/%s/mass_scan/%s --cutUnderOverFlow"%(v,m,v,m))
        #os.system("python scripts/runDileptonUnfolding.py -r Results/%s/mass_scan/%s/MC8TeV_TTJets_MSDecays_%sv5.root -v %s -o  Results/%s/mass_scan/%s -g f -w b"%(v,m,m,v,v,m)) 
        #os.system("python scripts/runDileptonUnfolding.py -r Results/%s/mass_scan/%s/plots/plotter.root -v %s -o  Results/%s/mass_scan/%s -g f -w b"%(v,m,v,v,m)) 

        #root = 'Results/%s/mass_scan/%s/plots/%s_unfolded_wgt.root'%(v,m,v)
      #  root = 'unfoldResults/%s/mass_scan/%s/plots/%s_unfolded.root'%(v,m,v)
        #getMoments(root,v,'unf',m)
        #getMoments(root,v,'gen',m)
        #root = 'unfoldResults/%s/mass_scan/%s/plots/plotter.root'%(v,m)
        #getMoments(root,v,'rec',m)

      #graphMoments(v)

      # getRatioMass('unf',v)
      # getRatioMass('gen',v)
      # getRatioMass('rec',v)

def main():
     usage = 'usage: %prog [options]'
     parser = optparse.OptionParser(usage)
     parser.add_option('-v', '--var',
                          dest='var', 
                          default='ptpos',
                          help='Variable to unfold (note requires var_rec,var_gen,var_migration plots stored in the ROOT file [default: %default]')
     parser.add_option('-o', '--output',
                          dest='output', 
                          default='unfoldResults',                                                                       
                          help='Output directory [default: %default]')
     (opt, args) = parser.parse_args()

     ROOT.gStyle.SetOptStat(0)
     ROOT.gStyle.SetOptTitle(0)
     ROOT.gROOT.SetBatch(True)

     plotter(opt)

if __name__ == "__main__":
    sys.exit(main())
