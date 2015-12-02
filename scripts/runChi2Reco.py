#!/usr/bin/env python                                                                                       
from array import array
import os, sys,numpy
import ROOT
import optparse
import pprint
from ROOT import TMath

from UserCode.TopMassSecVtx.PlotUtils import *

"""                                                                                                                                                                                                                              
Chi square test
"""
def chitest():

    histos={}
    chisquares = {}
    MASS = ['166','169','171','172','173','175','178']
    #VARS = ['ptll','mll','ptpos','ptposptm','Epos','EposEm']
    #TAGS = ['lesup','lesdn','jesup','jesdn','jerup','jerdn','puup','pudn','btagup','btagdn','mistagup','mistagdn','lepselup','lepseldn']
    TAGS = ['powherw','matchingup','matchingdown','scaleup','scaledown','toppt']
    VARS = ['ptposptm']
    for v in VARS:

        FILES = [('Dileptons/%s/reco/MC8TeV_TTJets_MSDecays_166v5.root'%v,'166'),
                 ('Dileptons/%s/reco/MC8TeV_TTJets_MSDecays_169v5.root'%v,'169'),
                 ('Dileptons/%s/reco/MC8TeV_TTJets_MSDecays_171v5.root'%v,'171'),
                 ('Dileptons/%s/reco/MC8TeV_TTJets_MSDecays_172v5.root'%v,'172'),
                 ('Dileptons/%s/reco/MC8TeV_TTJets_MSDecays_173v5.root'%v,'173'),
                 ('Dileptons/%s/reco/MC8TeV_TTJets_MSDecays_175v5.root'%v,'175'),
                 ('Dileptons/%s/reco/MC8TeV_TTJets_MSDecays_178v5.root'%v,'178')]

        minx_172 = 0
        
        chi2calibration  = []
        err_1 = []

        for fil,mass in FILES:
            filename = ROOT.TFile.Open(fil)
            hist_ref = filename.Get('rec')
            hist_ref.SetDirectory(0)
            filename.Close()
            
            print 'Having %s fixed for %s'%(mass,v)
            chisquare = [] 
            for f,m in FILES:
                fIn = ROOT.TFile.Open(f)
                hist = fIn.Get('gen')
                hist.SetDirectory(0)
                fIn.Close()
                chi2 = hist.Chi2Test(hist_ref,"CHI2")
                chisquare.append(chi2)
                print '%s %3.4f'%(m,chi2)          
            outputName = 'Dileptons/%s/chi_test_%s'%(v,mass)
            minx,b,c1 = graphChi2(v,mass,chisquare,outputName)
        
            if mass != '178' and mass != '166':
                print 'm_{t}=%.2f + %.2f - %.2f'%(minx,b,c1)
                chi2calibration.append(minx)
                err_1.append(b)
            
            if mass == '172':
                minx_172 = minx
                print minx_172

        print chi2calibration
        print err_1

        graph(v,chi2calibration,err_1)

        # Calculate systematics for shape measurement at reco level
        systematics = [] 

        for tag in TAGS:
            fil = 'Dileptons/%s/reco/DileptonKin_%s.root'%(v,tag)
            filename = ROOT.TFile.Open(fil)
            hist_ref = filename.Get('rec')
            hist_ref.SetDirectory(0)
            filename.Close()

            print 'Having %s fixed for %s'%(tag,v)
            chisquare = []
            for f,m in FILES:
                fIn = ROOT.TFile.Open(f)
                hist = fIn.Get('gen')
                hist.SetDirectory(0)
                fIn.Close()
                chi2 = hist.Chi2Test(hist_ref,"CHI2")
                chisquare.append(chi2)
                print '%s %3.4f'%(m,chi2)
            mass = '172'
            outputName = 'Dileptons/%s/systs/chi_test_%s'%(v,tag)
            minx,b,c1 = graphChi2(v,mass,chisquare,outputName)

            print 'm_{t}=%.2f + %.2f - %.2f'%(minx,b,c1)
            print 'nominal value = %.2f'%minx_172
            diff = minx-minx_172
            print '%s %s %.2f %.2f %.2f'%(v,tag,minx,minx_172,diff)
            systematics.append((tag,diff))

        diffplus = 0
        diffmin = 0
        for tag,val in systematics:
            print '%s %2.3f'%(tag,val)
            if val >= 0: diffplus += val**2
            if val < 0: diffmin += val**2
        print 'diffplus %2.3f'%sqrt(diffplus)
        print 'diffminus %2.3f'%sqrt(diffmin)
        print err_1

"""                                                                                                                                                                                                                             
Graph Chi2 pol2 plot                                                                                                                                                                                                             
"""
def graphChi2(var,m,chi2,outputName):

    #pT positive lepton                                                                                                                                                                                                         
    if var == 'ptpos': title = 'p_{T}(l^{+})'
    #pT positive lepton                                                                                                                                                                                                         
    if var == 'Epos': title = 'E(l^{+})'
    #pT charged-lepton pair                                                                                                                                                                           
    if var == 'ptll': title = 'p_{T}(l^{+}l^{-})'
    #M charged-lepton pair                                                                                                                                                                                                      
    if var == 'mll': title = 'M(l^{+}l^{-})'
    #Scalar sum of E                                                                                                                                                                                                            
    if var == 'EposEm': title = 'E(l^{+})+E(l^{-})'
    #Scalar sum of Pt                                                                                                                                                                                                           
    if var == 'ptposptm': title = 'p_{T}(l^{+})+p_{T}(l^{-})'
    
    mass = [166.5, 169.5, 171.5, 172.5, 173.5, 175.5 , 178.5]
    # mass = [169.5, 171.5, 172.5, 173.5, 175.5]

    mass = array('d', mass)
    chi2 = array('d', chi2)
        
    #Create a canvas for plotting your graph                                                                                                                                                                                 
    c = ROOT.TCanvas('c','c')
    ROOT.gStyle.SetCanvasDefH(600);
    ROOT.gStyle.SetCanvasDefW(600);
    c.SetLeftMargin(0.15);
    c.SetRightMargin(0.25)
    c.SetBottomMargin(0.25);
    pad1=ROOT.TPad('p1','p1',0.,0.,1.0,1.0)
    pad1.Draw()

    graphs = {} 
    graphs['chi2'] = ROOT.TGraph(len(mass),mass,chi2)
       
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    graphs['chi2'].SetMarkerColor(ROOT.kBlack)
    graphs['chi2'].SetMarkerStyle(1)
    graphs['chi2'].SetMarkerSize(2)
    graphs['chi2'].SetLineColor(ROOT.kBlack)
    
    pol2 = ROOT.TF1("pol2", "pol2", 166.5, 178.5)
    pol2.SetLineWidth(4)
    pol2.SetLineStyle(1)
    pol2.SetLineColor(ROOT.kBlack)
    graphs['chi2'].Fit("pol2","R")
    minx = pol2.GetMinimumX(166.5,178.5)
    miny = pol2.GetMinimum(166.5,178.5)
    by = miny +1
    print 'minx %f'%minx
    print 'miny %f'%miny
    print 'by%f'%by

    pol2 = graphs['chi2'].GetFunction("pol2")
    p0 = pol2.GetParameter(0)
    p1 = pol2.GetParameter(1)
    p2 = pol2.GetParameter(2)
    p0 = p0- by
    b = (-p1-sqrt(p1**2-4*p2*p0))
    b = minx - b/(2*p2)
    c1 = (-p1+sqrt(p1**2-4*p2*p0))
    c1 = c1/(2*p2) - minx
    ROOT.gStyle.SetOptFit(0)
    tlat3 = TLatex()
    tlat3.SetNDC()
    tlat3.SetTextFont(42)
    tlat3.SetTextSize(0.03)
    tlat3.SetTextAlign(31)
    tlat3.DrawLatex(0.8,0.24,'Fit Results')
    tlat3.DrawLatex(0.8,0.19,'m_{t}=%.2f + %.2f - %.2f'%(minx,b,c1))
    
    graphs['chi2'].Draw("AP")
    graphs['chi2'].GetXaxis().SetTitle("Top Mass [GeV]")
    graphs['chi2'].GetYaxis().SetTitle("#chi^{2} %s"%title)
    graphs['chi2'].GetXaxis().SetTitleSize(0.037)
    graphs['chi2'].GetXaxis().SetLabelSize(0.037)
    graphs['chi2'].GetYaxis().SetTitleOffset(1.2)
    graphs['chi2'].GetYaxis().SetTitleSize(0.037)
    graphs['chi2'].GetYaxis().SetLabelSize(0.037)

    #for ext in ['png','pdf']:
    #    c.SaveAs('%s.%s'%(outputName,ext))
    #del c

    return minx,b,c1

"""                                                                                                                                                                                                                             
Graph Calibration chi2                                                                                        
"""
def graph(var,moments_unf,sigma_unf):

    #pT positive lepton                                                                                                                                                                                                         
    if var == 'ptpos': title = 'p_{T}(l^{+})'
    #pT positive lepton                                                                                                                                                                                                         
    if var == 'Epos': title = 'E(l^{+})'
    #pT charged-lepton pair                                                                                                                                                                                                  
    if var == 'ptll': title = 'p_{T}(l^{+}l^{-})'
    #M charged-lepton pair                                                                                                                                                                                                    
    if var == 'mll': title = 'M(l^{+}l^{-})'
    #Scalar sum of E                                                                                                                                                                                                           
    if var == 'EposEm': title = 'E(l^{+})+E(l^{-})'
    #Scalar sum of Pt                                                                                                                                                                                                          
    if var == 'ptposptm': title = 'p_{T}(l^{+})+p_{T}(l^{-})'

    mass = [169.5, 171.5, 172.5, 173.5, 175.5]

    mass = array('d', mass)
    moments_unf = array('d', moments_unf)
    sigma_mass = array('d', [0.0]*len(mass))
    sigma_unf = array('d', sigma_unf)

    c = ROOT.TCanvas('c','c')
    ROOT.gStyle.SetCanvasDefH(600);
    ROOT.gStyle.SetCanvasDefW(600);
    pad1=ROOT.TPad('p1','p1',0.,0,1.0,1.0)
    pad1.Draw()

    tlat = TLatex()
    tlat.SetNDC()
    tlat.SetTextFont(61)
    tlat.SetTextSize(0.04)
    tlat.SetTextAlign(31)
    prelim_text1 = 'CMS'
    tlat.DrawLatex(0.25, 0.95, prelim_text1)

    tlat1 = ROOT.TLatex()
    tlat1.SetNDC()
    tlat1.SetTextFont(42)
    tlat1.SetTextSize(0.040)
    tlat1.SetTextAlign(31)
    prelim_text2 ='#it{Preliminary}'
    tlat1.DrawLatex(0.46, 0.95, prelim_text2)
    prelim_text3 = '19.7 fb^{-1}'
    tlat1.DrawLatex(0.85, 0.95, prelim_text3)
    prelim_text4 ='(8 TeV)'
    tlat1.DrawLatex(0.97, 0.95, prelim_text4)

    graphs = {}
    multigraphs = {}

    graphs['moments_unf'] = ROOT.TGraphErrors(len(mass),mass,moments_unf,sigma_mass,sigma_unf)
    multigraphs['u1'] = ROOT.TMultiGraph()

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)

    graphs['moments_unf'].SetMarkerStyle(8)
    graphs['moments_unf'].SetMarkerSize(1)

    linear = ROOT.TF1("linear", "pol1")
    linear.SetLineColor(ROOT.kRed)
    linear.SetLineStyle(1)
    graphs['moments_unf'].Fit(linear, "WQ")

    linear = graphs['moments_unf'].GetFunction("linear")
    p0_unf=linear.GetParameter(0)
    p0Err_unf=linear.GetParError(0)
    p1_unf=linear.GetParameter(1)
    p1Err_unf=linear.GetParError(1)
    chi2_unf=linear.GetChisquare()

    ROOT.gStyle.SetOptFit(0)
    tlat3 = TLatex()
    tlat3.SetNDC()
    tlat3.SetTextFont(42)
    tlat3.SetTextSize(0.02)
    tlat3.SetTextAlign(31)
    tlat3.SetTextColor(ROOT.kBlack)
    tlat3.DrawLatex(0.60,0.81,'Unf. (%3.3f #pm %3.3f) + (%3.3f #pm %3.3f) m_{top}'%(p0_unf,p0Err_unf,p1_unf,p1Err_unf))
    print 'Unf. (%3.3f #pm %3.3f) + (%3.3f #pm %3.3f)*m_{top}'%(p0_unf,p0Err_unf,p1_unf,p1Err_unf)

    multigraphs['u1'].Add(graphs['moments_unf'])

    pad1.cd()
    multigraphs['u1'].Draw("AP")
    multigraphs['u1'].GetXaxis().SetTitle("Top Mass [GeV]")
    multigraphs['u1'].GetYaxis().SetTitle("#chi^{2} %s [GeV]"%(title))
    multigraphs['u1'].GetXaxis().SetTitleSize(0.037)
    multigraphs['u1'].GetXaxis().SetLabelSize(0.037)
    multigraphs['u1'].GetYaxis().SetTitleOffset(1.6)
    multigraphs['u1'].GetYaxis().SetTitleSize(0.037)
    multigraphs['u1'].GetYaxis().SetLabelSize(0.037)

    outputName = 'Dileptons/%s/calibration_chi2_reco'%var
    #for ext in ['png','pdf']:
    #    c.SaveAs('%s.%s'%(outputName,ext))
    #del c


"""                                                                                                      
steer                                                                                                    
"""
def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gROOT.SetBatch(True)
    setTDRStyle()
    ROOT.gSystem.Load("libUserCodeTopMassSecVtx")
    ROOT.gSystem.Load('libGenVector')
    ROOT.gSystem.Load('libSmatrix')
    ROOT.AutoLibraryLoader.enable()

    print 80*'='
    chitest()
    print 80*'='
    
    return 0

if __name__ == "__main__":
	sys.exit(main())
