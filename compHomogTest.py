#! /bin/env python


# compHomogTest.py: a script to perform compositional homogeineity test using both simulation and Chi-squared test.
# Please type python compHomogTest.py -h to see the usage.


#####################################################################
import sys
import getopt
import commands

from p4 import *


#####################################################################
infiles = {'tree':None, 'data':None}
nSims = 100


#####################################################################
class ModelArgu(object):
    def __init__(self):
        self.invar_prop = 0.0
        self.compModel = 'empirical'
        self.compIsFree = False
        self.rateModel = 'wag'
        self.rateIsFree = False


class ConstructTree(object):
    def __init__(self):
        pass


#####################################################################
def read_infiles(infiles, constructTree):
    read(infiles['data'])
    d = Data()

    var.alignments[0].checkForDuplicateSequences(removeDupes=False)

    if hasattr(constructTree, 'method'):
        if constructTree.method == 'bionj':
            status, output = commands.getstatusoutput('which ' + constructTree.method)
            if status != 0:
                print "%s%s%s%s" % ('Sorry! The method to construct the tree used in simulation ', constructTree.method, ' is not available.', ' Exiting ......\n')
                show_help()
            dm = var.alignments[0].pDistances()
            t = dm.bionj()
        else:
            print "Sorry. Only bionj is supported in tree construction now. Alternatively, you need to provide a tree ('-t')."
            print
            show_help()
    else:
        read(infiles['tree'])
        t = var.trees[0]

    return([t,d])


def examine_required_info(infiles, constructTree):
    for i in ['tree', 'data']:
        if i == 'tree' and constructTree:
            continue
        if not infiles[i] or not os.path.isfile(infiles[i]):
            print i, "has to be specified. Exiting ......"
            sys.exit()


def show_help():
    print "compHomogTest: a simple script to test compositional homogeneity."
    print "python", sys.argv[0], '<-d|--data datafile>', '<-t|--tree treefile>';
    print 
    print "Optional arguments:"
    print '%-30s%s' % ("[-h|--h|--help]", 'Help message')
    print '%-30s%s' % ("[--invar_prop]", 'Proportion of invariable sites (0.0)')
    print '%-30s%s' % ("[--nSims]", 'No. of simulations (100)')
    print '%-30s%s' % ("[--gamma]", 'Shape parameter of Gamma distribution (Off)')
    print '%-30s%s' % ("[-c|--constructTree]", 'The method of tree construction used in simulation: bionj (Off)')
    print '%-30s%s' % ("[--compFree|--compIsFree]", 'Compositional frequencies are set to be unfixed (No)')
    print '%-30s%s' % ("[--compModel]", "Substitution model: 'equal', 'empirical', 'specified', 'cpREV', 'd78', 'jtt', 'mtREV24', 'mtmam', 'wag', 'rtRev', 'tmjtt94', 'tmlg99', 'lg', 'blosum62', 'hivb', 'mtart', 'mtzoa', 'gcpREV', 'stmtREV' (empirical)")
    print '%-30s%s' % ("[--rateFree|--rateIsFree]", 'Substitution rates are set to be variable (No)')
    print '%-30s%s' % ("[--rateModel]", "Substitution model: 'ones', '2p', 'specified', 'optimized', 'cpREV', 'd78', 'jtt', 'mtREV24', 'mtmam', 'wag', 'rtRev', 'tmjtt94', 'tmlg99', 'lg', 'blosum62', 'hivb', 'mtart', 'mtzoa', 'gcpREV', 'stmtREV' (wag)")
    print
    print '%-30s%s' % ('Author', 'Sishuo Wang @ Luo Haiwei Lab, Chinese University of Hong Kong. Please do not hesitate to write to sishuowang@hotmail.ca or sishuowang@cuhk.edu.hk for any bug report or suggestion. Your help is highly appreciated.')
    print '%-30s%s' % ('License', 'BSD')
    print '\n'
    sys.exit()


#####################################################################
modelArgu = ModelArgu()
constructTree = ConstructTree()


try:
    opts, args = getopt.getopt(sys.argv[1:], 't:d:c:', ['data=', 'tree=', 'gamma=', 'invar_prop=', 'nSims=', 'constructTree=', \
        'compFree', 'compIsFree', 'compModel=', \
        'rateFree', 'rateIsFree', 'rateModel='])
except getopt.GetoptError:
    print "Illegal params!"
    show_help()

if not opts:
    show_help()

for opt, value in opts:
    if opt == '-t' or opt == '--tree':
        infiles['tree'] = value
    elif opt == '-d' or opt == '--data':
        infiles['data'] = value
    elif opt == '--gamma':
        modelArgu.gamma = float(value)
    elif opt == '--invar_prop':
        modelArgu.invar_proFp = float(value)
    elif opt == '--compFree' or opt == '--compIsFree':
        modelArgu.compIsFree = True
    elif opt == '--compModel':
        modelArgu.compModel = value
    elif opt == '--rateFree' or opt == '--rateIsFree':
        modelArgu.rateIsFree = True
    elif opt == '--rateModel':
        modelArgu.rateModel = value
    elif opt == '--nSims':
        nSims = int(value)
    elif opt == '-c' or opt == '--constructTree':
        constructTree.method = value
    elif opt == '-h' or opt == '--h' or opt == '--help':
        show_help()


examine_required_info(infiles, constructTree)


#####################################################################
t, d = read_infiles(infiles, constructTree)


# Attach the data and a model to the tree.
t.data = d
t.newComp(free=modelArgu.compIsFree, spec=modelArgu.compModel)
t.newRMatrix(free=modelArgu.rateIsFree, spec=modelArgu.rateModel)
t.setPInvar(free=0, val=modelArgu.invar_prop)

if hasattr(modelArgu, 'gamma'):
    t.setNGammaCat(nGammaCat=4) # No. of GammaCat is set to 4.
    t.newGdasrv(free=0, val=modelArgu.gamma)

# heterogeneous model
# t.setModelThingsRandomly()

# optimize
t.optLogLike()

t.compoTestUsingSimulations(nSims=nSims, doIndividualSequences=0, doChiSquare=True, verbose=1)

