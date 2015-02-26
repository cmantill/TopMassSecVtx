#! /usr/bin/env python
import os, sys, re
import ROOT
import pickle
from UserCode.TopMassSecVtx.PlotUtils import RatioPlot

# MASSES = [163.5, 166.5, 169.5, 171.5, 172.5, 173.5, 175.5, 178.5, 181.5]
NBINS = 100
XMIN = 5.
XMAX = 200.
FITRANGE = (20., 140.)
MASSXAXISTITLE = 'm(SV,lepton) [GeV]'

#NTRKBINS = [(2,3), (3,4), (4,5), (5,7) ,(7,1000)]
NTRKBINS = [(2,3), (3,4), (4,1000)]
COMMONWEIGHT = "Weight[1]*Weight[4]*JESWeight[0]"

TREENAME = 'SVLInfo'
SELECTIONS = [
	('inclusive', 'abs(EvCat)<200', '#geq 1 lepton'),
	('ee',        'EvCat==-121',    'ee'),
	('em',        'EvCat==-143',    'e#mu'),
	('mm',        'EvCat==-169',    '#mu#mu'),
	('e',         'abs(EvCat)==11', 'e'),
	('m',         'abs(EvCat)==13', '#mu'),
	('inclusive_mrank1', 'SVLMassRank==1&&SVLDeltaR<2.0&&CombCat%2!=0&&abs(EvCat)<200', '#geq 1 lepton'),
	('ee_mrank1',        'SVLMassRank==1&&SVLDeltaR<2.0&&CombCat%2!=0&&EvCat==-121', 'ee'),
	('em_mrank1',        'SVLMassRank==1&&SVLDeltaR<2.0&&CombCat%2!=0&&EvCat==-143', 'e#mu'),
	('mm_mrank1',        'SVLMassRank==1&&SVLDeltaR<2.0&&CombCat%2!=0&&EvCat==-169', '#mu#mu'),
	('e_mrank1',         'SVLMassRank==1&&SVLDeltaR<2.0&&CombCat%2!=0&&abs(EvCat)==11', 'e'),
	('m_mrank1',         'SVLMassRank==1&&SVLDeltaR<2.0&&CombCat%2!=0&&abs(EvCat)==13', '#mu'),
]

COMBINATIONS = [
	('tot', '1'),
	('cor', 'CombInfo==1'),
	('wro', 'CombInfo==0'),
	('unm', 'CombInfo==-1'),
]

CHANNAMES = {
	'tt':'TTbar',
	't':'SingleT_t',
	'tW':'SingleT_tW',
	's':'SingleT_s'
}

def addTopMassTreesFromDir(dirname, trees, files):
	for filename in os.listdir(dirname):
		if not os.path.splitext(filename)[1] == '.root': continue

		## Select top samples
		if not "TTJets" in filename and not "SingleT" in filename: continue

		## Extract the mass (and splitting index, if it's there)
		matchres = re.match(r'.*([0-9]{3})v5\_?([0-9]+)?\.root', filename)
		if matchres: mass = float(matchres.group(1)) + 0.5
		else: mass = 172.5

		treefile = os.path.join(dirname, filename)
		tfile = ROOT.TFile.Open(treefile,'READ')
		tree = tfile.Get(TREENAME)

		chan = 'tt'

		if 'SingleT' in filename and '_tW' in filename:
			chan = 'tW'
		elif 'SingleT' in filename and '_t' in filename:
			chan = 't'
		elif 'SingleT' in filename and '_s' in filename:
			continue

		if not (mass,chan) in trees:
			trees[(mass, chan)] = ROOT.TChain(TREENAME)
			files[(mass, chan)] = []
		trees[(mass,chan)].Add(treefile)
		files[(mass,chan)].append(treefile)
def getMassTrees(inputdir, verbose=True):
	alltrees = {} # (mass,chan) -> chain
	allfiles = {} # (mass,chan) -> [locations]

	addTopMassTreesFromDir(inputdir, alltrees, allfiles)
	addTopMassTreesFromDir(os.path.join(inputdir, 'mass_scan'),
									 alltrees, allfiles)

	allmasses = sorted(list(set([mass for mass,_ in alltrees.keys()])))

	if verbose:
		print 80*'-'
		print "Found the following mass points (#entries):"
		for mass in sorted(allmasses):
			line = '  %5.1f GeV '%mass
			if (mass, 'tt') in alltrees:
				line += ' (tt: %7d) '     % alltrees[(mass, 'tt')].GetEntries()
			if (mass, 't') in alltrees:
				line += ' (t: %5d) '  % alltrees[(mass, 't')].GetEntries()
			if (mass, 'tW') in alltrees:
				line += ' (tW: %6d) ' % alltrees[(mass, 'tW')].GetEntries()
			print line
		print 80*'-'

	return alltrees, allfiles

def runSVLInfoTreeAnalysis((treefiles, histos, outputfile)):
	taskname = os.path.basename(outputfile)[:-5]
	chain = ROOT.TChain(TREENAME)
	for filename in treefiles:
		if not os.path.exists(filename):
			print "ERROR: file %s does not exist! Aborting" % filename
			return -1
		chain.Add(filename)
	print ' ... processing %s for %d histos from %7d entries' %(
		              taskname, len(histos), chain.GetEntries())

	from ROOT import gSystem
	gSystem.Load('libUserCodeTopMassSecVtx.so')
	from ROOT import SVLInfoTreeAnalysis
	ana = SVLInfoTreeAnalysis(chain)

	for hname,var,sel,nbins,xmin,xmax,xtitle in histos:
		ana.AddPlot(hname, var, sel, nbins, xmin, xmax, xtitle)
	ana.RunJob(outputfile)
	print '         %s done' % taskname

def runTasks(inputfiles, tasklist, opt, subdir):
	assert set(inputfiles.keys()) == set(tasklist.keys())
	tasks = []

	os.system('mkdir -p %s' % os.path.join(opt.outDir, subdir))

	for key, task in tasklist.iteritems():
		for ifilep in inputfiles[key]:
			ofilen = os.path.basename(ifilep)
			ofilep = os.path.join(opt.outDir, subdir, ofilen)
			tasks.append(([ifilep], task, ofilep))

	if opt.jobs > 1:
		import multiprocessing as MP
		pool = MP.Pool(opt.jobs)
		pool.map(runSVLInfoTreeAnalysis, tasks)
	else:
		for task in tasks:
			runSVLInfoTreeAnalysis(task)

def gatherHistosFromFiles(tasklist, massfiles, dirname, hname_to_keys):
	## First extract a list of ALL histogram names from the tasklist
	hnames = [t[0] for tasks in tasklist.values() for t in tasks]

	histos = {}
	for ifilen in os.listdir(dirname):
		if not os.path.splitext(ifilen)[1] == '.root': continue
		ifile = ROOT.TFile.Open(os.path.join(dirname,ifilen), 'READ')

		for hname in hnames:
			keys = hname_to_keys[hname]
			tag = keys[0]
			chan = keys[1]
			mass = keys[2]

			## Skip histo names that are not in this file
			tfilens = [os.path.basename(x) for x in massfiles[(mass,chan)]]
			if not ifilen in tfilens: continue

			try:
				histo = ifile.Get(hname)
				try:
					histo.SetDirectory(0)

					## Add histograms with same names
					if not keys in histos:
						histos[keys] = histo
					else:
						histos[keys].Add(histo)
				except AttributeError: pass

			except ReferenceError:
				print hname, "not found in", ifilen
	return histos

def main(args, opt):
	os.system('mkdir -p %s'%opt.outDir)
	masstrees, massfiles = getMassTrees(args[0], verbose=True)
	masspoints = sorted(list(set([mass for mass,_ in masstrees.keys()])))

	hname_to_keys = {} # hname -> (tag, chan, mass, comb, [ntk1])

	tasklist = {} ## (mass,chan) -> tasks
	for (mass,chan) in masstrees.keys():
		tasks = []
		for tag,sel,_ in SELECTIONS:

			htag = ("%s_%5.1f"%(tag,mass)).replace('.','')
			if not chan == 'tt':
				htag = ("%s_%s_%5.1f"%(tag,chan,mass)).replace('.','')

			for comb,combsel in COMBINATIONS:
				hname = "SVLMass_%s_%s" % (comb, htag)
				finalsel = "%s*(%s&&%s)"%(COMMONWEIGHT,sel,combsel)
				tasks.append((hname, 'SVLMass', finalsel,
					          NBINS, XMIN, XMAX, MASSXAXISTITLE))
				hname_to_keys[hname] = (tag, chan, mass, comb)

				for ntk1,ntk2 in NTRKBINS:
					tksel = "(SVNtrk>=%d && SVNtrk<%d)"%(ntk1,ntk2)
					finalsel = "%s*(%s&&%s&&%s)"%(COMMONWEIGHT, sel,
						                          combsel,tksel)
					hname = "SVLMass_%s_%s_%d" % (comb, htag, ntk1)
					tasks.append((hname, 'SVLMass', finalsel,
						          NBINS, XMIN, XMAX, MASSXAXISTITLE))
					hname_to_keys[hname] = (tag, chan, mass, comb, ntk1)

		tasklist[(mass,chan)] = tasks

	if not opt.cache:
		runTasks(massfiles, tasklist, opt, 'mass_histos')

	## Retrieve the histograms from the individual files
	# (tag, chan, mass, comb)      -> histo
	# (tag, chan, mass, comb, ntk) -> histo
	masshistos = gatherHistosFromFiles(tasklist, massfiles,
		                              os.path.join(opt.outDir,
		                              'mass_histos'),
		                              hname_to_keys)

	cachefile = open(".svlmasshistos.pck", 'w')
	pickle.dump(masshistos, cachefile, pickle.HIGHEST_PROTOCOL)
	cachefile.close()

	# ofi = ROOT.TFile(os.path.join(opt.outDir,'masshistos.root'),
	# 													   'recreate')
	# ofi.cd()

	# for key in masshistos.keys():
	# 	tag, chan, mass = key[0],key[1],key[2]
	# 	if not ofi.cd(tag):
	# 		outDir = ofi.mkdir(tag)
	# 		outDir.cd()

	# 	for comb,_ in COMBINATIONS:
	# 		masshistos[(tag,chan,mass,comb)].Write()
	# 		for ntk,_ in NTRKBINS:
	# 			masshistos[(tag,chan,mass,comb,ntk)].Write()
	# ofi.Write()
	# ofi.Close()

	ROOT.gStyle.SetOptTitle(0)
	ROOT.gStyle.SetOptStat(0)
	ROOT.gROOT.SetBatch(1)

	errorGetters = {} # tag -> function(chi2) -> error
	systematics = {} # (seltag, systname) -> error
	mass_scan_dir = os.path.join(opt.outDir, 'mass_scans')
	for tag,sel,seltag in SELECTIONS:
		print "... processing %s"%tag

		for chan in ['tt', 't', 'tW']:
			# print "   %s channel" % chan
			## Skip some useless combinations:
			if chan == 't' and ('ee' in tag or 'mm' in tag or 'em' in tag):
				continue
			if chan == 'tW' and tag in ['e', 'm', 'e_mrank1', 'm_mrank1']:
				continue

			ratplot = RatioPlot('ratioplot')
			ratplot.ratiotitle = "Ratio wrt 172.5 GeV"
			ratplot.ratiorange = (0.5, 1.5)

			ratplot.reference = masshistos[(tag,chan,172.5,'tot')]

			for mass in masspoints:
				legentry = 'm_{t} = %5.1f GeV' % mass
				try:
					histo = masshistos[(tag,chan,mass,'tot')]
					# print histo.GetName(), histo.GetEntries()
					ratplot.add(histo, legentry)
				except KeyError: pass
					# print "Can't find ", (tag,chan,mass,'tot')


			ratplot.tag = 'All combinations'
			ratplot.subtag = '%s %s' % (seltag, chan)
			ratplot.show("massscan_%s_%s_tot"%(tag,chan), mass_scan_dir)

			# chi2s = ratplot.getChiSquares(rangex=FITRANGE)
			# chi2stofit = []
			# for legentry in sorted(chi2s.keys()):
			# 	chi2stofit.append((float(legentry[8:-4]), chi2s[legentry]))
			# errorGetters[tag] = fitChi2(chi2stofit,
			# 							tag=seltag,
			# 							oname=os.path.join(mass_scan_dir,
			# 						    "chi2_simple_fit_%s.pdf"%tag),
			# 							drawfit=False)
			ratplot.reset()

			ratplot.reference = masshistos[(tag,chan,172.5,'cor')]
			for mass in masspoints:
				legentry = 'm_{t} = %5.1f GeV' % mass
				try: ratplot.add(masshistos[(tag,chan,mass,'cor')], legentry)
				except KeyError: pass
			ratplot.tag = 'Correct combinations'
			ratplot.subtag = '%s %s' % (seltag, chan)
			ratplot.show("massscan_%s_%s_cor"%(tag,chan), mass_scan_dir)
			ratplot.reset()

			ratplot.reference = masshistos[(tag,chan,172.5,'wro')]
			for mass in masspoints:
				legentry = 'm_{t} = %5.1f GeV' % mass
				try: ratplot.add(masshistos[(tag,chan,mass,'wro')], legentry)
				except KeyError: pass
			ratplot.tag = 'Wrong combinations'
			ratplot.subtag = '%s %s' % (seltag, chan)
			ratplot.show("massscan_%s_%s_wro"%(tag,chan), mass_scan_dir)
			ratplot.reset()

			ratplot.reference = masshistos[(tag,chan,172.5,'unm')]
			for mass in masspoints:
				legentry = 'm_{t} = %5.1f GeV' % mass
				try: ratplot.add(masshistos[(tag,chan,mass,'unm')], legentry)
				except KeyError: pass
			ratplot.tag = 'Unmatched combinations'
			ratplot.subtag = '%s %s' % (seltag, chan)
			ratplot.show("massscan_%s_%s_unm"%(tag,chan), mass_scan_dir)
			ratplot.reset()

	return 0



if __name__ == "__main__":
	import sys
	tmpargv  = sys.argv[:]     # [:] for a copy, not reference
	sys.argv = []
	from ROOT import gROOT, gStyle
	sys.argv = tmpargv
	from optparse import OptionParser
	usage = """
	usage: %prog [options] input_directory
	"""
	parser = OptionParser(usage=usage)
	parser.add_option('-v', '--verbose', dest='verbose', action="store",
					  type='int', default=1,
					  help='Verbose mode [default: %default (semi-quiet)]')
	parser.add_option('-j', '--jobs', dest='jobs', action="store",
					  type='int', default=1,
					  help=('Number of jobs to run in parallel '
					  	    '[default: single]'))
	parser.add_option('-o', '--outDir', dest='outDir', default='svlplots',
					  help='Output directory [default: %default]')
	parser.add_option('-c', '--cache', dest='cache', action="store_true",
					  help='Read from cache')
	(opt, args) = parser.parse_args()

	exit(main(args, opt))




