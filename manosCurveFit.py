''' manosCurveFit.py - Hosea Siu 2014
    Takes light curve data as input and fits it to a function according to Harris et al. 1989
'''

from lmfit import minimize, Parameters, Parameter, report_fit, fit_report
from operator import itemgetter
from os import listdir
import matplotlib.pyplot as plt, numpy as np, os.path, sys, string, cmd
from time import clock
from uncertainties import ufloat
from uncertainties.umath import *

basepath = os.path.abspath(os.path.dirname(sys.argv[0]))        # current directory path
objectName = ''

class RunOptionsShell(cmd.Cmd):
    intro = 'Welcome to manosCurveFit. Type help or ? to list commands.\n'
    prompt = '(manosCurveFit) '

    # default options, can be overridden
    fitOptions = {'minOrder':2, 'maxOrder':6, 'timer':False}
    outputOptions = {'printReport':True, 'saveReport':True, 'plotFullPeriod':True, 'plotErrorbars':True, \
                    'phaseFoldData':True, 'plotResiduals':True, 'plotPeriodErrors':True, 'showPlots':False}
   
    def do_setFitOptions(self, arg):       # sets options and returns other arguments as needed
        'Set the fitting options: minOrder, maxOrder, timer'
        args = arg.split()
        for i in range(len(self.fitOptions.keys())):
            key = self.fitOptions.keys()[i]
            if key in arg:
                if key == 'minOrder':
                    if int(args[args.index(key)+1]) < 0 or int(args[args.index(key)+1]) > self.fitOptions['maxOrder']:
                        print 'Error: order less than zero or greater than maximum'
                    else: 
                        self.fitOptions[key] = int(args[args.index(key)+1])                    
                elif key == 'maxOrder':
                    if int(args[args.index(key)+1]) < 0 or int(args[args.index(key)+1]) < self.fitOptions['minOrder']:
                        print 'Error: order less than zero or less than minimum'
                    else:
                        self.fitOptions[key] = int(args[args.index(key)+1])
                elif args[args.index(key)+1].lower() == 'true' or args[args.index(key)+1].lower() == 'false':       # check that the argument is a boolean
                    self.fitOptions[key] = args[args.index(key)+1].lower() == 'true'
                else:
                    print 'Error: unrecognized argument for ' + key
                print key + ' = ' + str(self.fitOptions[key])

    def do_setOutputOptions(self, arg):       # sets options and returns other arguments as needed
        'Set the output options: printReport, saveReport, plotFullPeriod, plotErrorbars, phaseFoldData, plotResiduals, plotPeriodErrors, showPlots'
        args = arg.split()
        for i in range(len(self.outputOptions.keys())):
            key = self.outputOptions.keys()[i]
            if key in arg:
                if args[args.index(key)+1].lower() == 'true' or args[args.index(key)+1].lower() == 'false':       # check that the argument is a boolean
                    self.outputOptions[key] = args[args.index(key)+1].lower() == 'true'
                else:
                    print 'Error: unrecognized argument for ' + key
                print key + ' = ' + str(self.outputOptions[key])

    def do_showObjects(self, arg):
        'Displays the list of all objects in \'Data\' folder'
        objects = lookInFolder('dir')[1]
        for i in range(len(objects)):
            print objects[i]

    def do_fitAll(self, arg):
        'Scans and attempts to fit all objects in \'Data\' folder; ignores previously fitted data by default, \'fitAll redo\' fits all objects regardless of any existing fits'
        objects = lookInFolder('dir')[1]
        args = arg.split()
        if len(args) != 0:
            if args[0] == 'redo':
                self.runFitting(objects, ignore = False)        # redo all objects regardless of whether or not they have a fit
        else:
            self.runFitting(objects)                        # only fit objects that haven't already been done

    def do_fit(self, arg):      # will fit all objects specified, even if they have already been fit before
        'Scans and attempts to fit all objects in \'Data\' folder that follow as arguments: fit object1_name object2_name...'
        self.runFitting(arg.split(), ignore = False)
        
    def do_exit(self, arg):
        'Exit the program'
        return True

    def runFitting(self, objects, ignore = True):      # ignore checks whether or not to attempt objects that have already been fit
        plt.close('all')     # close any remaining plots
        print '\nAttempting to fit ' + str(objects)
        for i in range(len(objects)):
            global objectName
            objectName = objects[i]
            directory, objectFiles = lookInFolder('file', objectName)
            if objectFiles is not None:     #check that the folder exists
                # check if the object has already been fit if ignore == True
                if ignore and objectName + 'LightCurve.png' in objectFiles:
                    print objectName + ' already processed- skipping to next object'
                # if it hasn't been fit, or is specified, run the fit routine
                else:
                    print '\nfitting ' + objectName

                    fileNamesAndFormat, offsets, guessMethod, T0, hardMinP, hardMaxP = extractRunOptions(objectName)

                    # note that these functions are all overloaded (there are extra options you can set- see the function definitions or documentation for examples)
                    lcd = LightCurveData(objectName, fileNamesAndFormat, offsetsList = offsets)         # read in data from the text file and create a LightCurveData object
                    (bestFit, m, periodsTested, periodErrors) = fitData(lcd, self.fitOptions, method = guessMethod, periodGuess = T0, hardMinPeriod = hardMinP, hardMaxPeriod = hardMaxP)       # fit the data to a model (can also add min and max periods in JD)
                    outputResults(bestFit, m, lcd, self.outputOptions, periodErrors = [periodsTested, periodErrors])             # plot and print the results
        print '\nDone\n'


''' Class used to store and manipulate MANOS light curve data
'''
class LightCurveData:
    ''' Gets the txt data files and converts them into a LightCurveData object
    '''
    def getData(self, fileName, formatSpec, offsetsList, passNumber):
        global basepath
        # TODO - have to change this for the MANOS file structure
        filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', fileName))
        textFile = np.loadtxt(filepath,'a12')          # extract data from text file
            
        # separate data into appropriately-named arrays - assumes that everything is a float (generalize?)
        tempData = {}

        for i in range(len(formatSpec)):
            tempData[formatSpec[i][0]] = textFile[:,formatSpec[i][1]].astype(float)

        if len(offsetsList) != 0:     # if offset values exist, then perform an offset
            tempData = self.offsetMags(tempData, offsetsList[passNumber])
            
        for i in range(len(formatSpec)):
            if passNumber == 0:
                self.data[formatSpec[i][0]] = tempData[formatSpec[i][0]]
            else:
                self.data[formatSpec[i][0]] = np.append(self.data[formatSpec[i][0]],tempData[formatSpec[i][0]])

    ''' Given a data type (i.e. 'jd', 'diffMag', etc), sorts the light curve data by that data type
        For fitting plots, we always sort by jd, though it may be switched to check other patterns
    '''
    def sortByDataType(self, key):
        if key is None:
            return
        
        dict = self.data
        keyIndex = dict.keys().index(key)
        numKeys = len(dict.keys())
        numEl = -1                  # checking to make sure that all dict entries have the same number of elements
        for i in range(numKeys):
            if numEl == -1:
                numEl = len(dict[dict.keys()[i]])
            else:
                if numEl != len(dict[dict.keys()[i]]):
                    print 'Error: invalid dictionary- nonrectangular'
                    return

        keyList = dict.keys()
        dictArray = np.zeros((numKeys,numEl))          # form a np array to use the transpose function
        for i in range(numKeys):
            dictArray[i] = dict[dict.keys()[i]]

        sortedList = sorted(dictArray.T.tolist(), key = itemgetter(keyIndex))

        temp = np.array(sortedList).T   # temporary variable for transpose
        for i in range(numKeys):
            dict[keyList[i]] = temp[i]
        self.data = dict

    ''' Given a dict of nights and corresponding magnitude offsets, offset the magnitude data by that much
        Also center the data by a weighted average
        "dataFile" refers to a data set from a specific file
    '''
    def offsetMags(self, dataFile, offsets = None):
        
        if offsets is None:
            return
        if len(np.unique(dataFile['night'])) == len(offsets):  # check if each night has an associated offset value
            for i in range(len(dataFile['diffMag'])):
                finishedOffset = False
                for j in range(len(offsets)):
                    if finishedOffset:
                        break
                    elif dataFile['night'][i] == offsets.keys()[j]:
                        dataFile['diffMag'][i] += offsets[offsets.keys()[j]]
            return dataFile
        else:
            print 'Error: number of nights does not match the number of offsets'
            return
            
    ''' LightCurveData initializer
        requires an object name, file names; can optionally take offsets and sorting
        (sorting should always be by jd for light curve plots)
    '''
    def __init__(self, objectName, fileNamesAndFormat, offsetsList = None):
        if fileNamesAndFormat is None:
            print 'Error: file name(s) and format(s) not specified'
            return

        # check for errors in the format specification
        for i in range(len(fileNamesAndFormat)):
            fileName = fileNamesAndFormat.keys()[i]
            formatSpec = fileNamesAndFormat[fileName]

            if len(formatSpec) <= 0:
                print 'Error: formatSpec input incorrect'
                return
            for j in range(len(formatSpec)):
                if len(formatSpec[j]) !=2:
                    print 'Error: formatSpec input incorrect'
                    return
        
        self.name = objectName
        self.data = {}      # initialize an empty dictionary with the keys in formatSpec
        for i in range(len(formatSpec)):
            self.data[formatSpec[i][0]] = {}

        # check for missing (required) specifications
        if 'jd' not in self.data.keys() or 'diffMag' not in self.data.keys() or 'magErr' not in self.data.keys():
            print 'Error: Insufficient data for light curve'
            print '       Did you include \'jd\', \'diffMag\', and \'magErr\'?'
            print '       Right now, I have ' + str(self.data.keys())
        else:
            for j in range(len(fileNamesAndFormat)):
                fileName = fileNamesAndFormat.keys()[j]
                formatSpec = fileNamesAndFormat[fileName]
                self.getData(fileName, formatSpec, offsetsList, j)
            self.sortByDataType('jd')
            self.data['diffMag'] -= np.average(self.data['diffMag'], weights = 1.0/self.data['magErr'])        # center the data around zero with a weighted average


''' Given a filepath to a directory, returns the filepath and list of all the files or subdirectories in that directory
'''
def lookInFolder(type, name = None):
    global basepath
    # TODO - have to change for MANOS file structure
    dirpath = os.path.abspath(os.path.join(basepath, "..", 'Data'))
    if type == 'file':
        if name in listdir(dirpath):            
            filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', name))
            return filepath, [ f for f in listdir(filepath) if os.path.isfile(os.path.join(filepath,f)) ]
        else:
            print name + ' does not exist in the directory'
            return None, None
    elif type == 'dir':
        return dirpath, [ f for f in listdir(dirpath) if os.path.isdir(os.path.join(dirpath,f)) ]
    else:
        print 'Error: invalid lookInFolder type'
        return None, None


''' Creates range lists of floats
    Has some error handling, but still has machine precision issues (shouldn't be an actual problem)
'''
def floatRange(start, stop, stepSize=0):    
    if start == stop:
        return []
    elif stepSize == 0:
        return [start, stop]
    elif (start < stop and stepSize > 0) or (start > stop and stepSize < 0):
        numbers = []
        start = float(start)
        stop = float(stop)
        current = start
        while True:
            numbers.append(current)
            if (start < stop):                  # allows for forward and backwards listing
                current = current + stepSize
                if current > stop:
                    break
            else:
                current = current - stepSize
                if current < stop:
                    break
        return numbers
    else:
        print 'Error: floatRange input error'
        return []


''' Generate the light curve model using Parameters
'''
def makeModel(params, t, mag=None, err=None):
    P = params['P'].value
    H = np.zeros(len(t))  
    for L in range(1,len(params.keys())/2):
        trigTerm = 2.0*np.pi*float(L)/float(P)
        offset = t-t[0]
        A = params['A'+str(L)].value
        B = params['B'+str(L)].value
        H += A*np.sin(trigTerm*offset)+B*np.cos(trigTerm*offset)
    H += params['y'].value              # add the y-shift offset
    if mag is None and err is None:     # no data given- the model values are returned
        return H
    elif err is None:                   # magnitude data is given, but no errors- the residuals are returned
        return H-mag
    else:                               # both magnitude and error are given- the normalized residuals are returned
        return (H-mag)/err


''' Generate the light curve model using uncertain values (for error propagation)
'''
def makeModelUncertainties(params, t):
    P = params['P']
    H = []
    for i in range(len(t)):
        H.append(ufloat(0,0))
    for L in range(1,len(params.keys())/2):
        trigTerm = 2.0*np.pi*float(L)/P
        offset = t-t[0]
        A = params['A'+str(L)]
        B = params['B'+str(L)]
        for i in range(len(H)):
            H[i] += A*sin(trigTerm*offset[i])+B*cos(trigTerm*offset[i])   # not np.math- using umath instead
    for i in range(len(H)):
        H[i] += params['y']                         # add the y-shift offset
    return H


''' Find the maximum recoverable period given a particular data set
    Use 5x the total observing window
'''
def findRecoverablePeriod(time):
    timeRange = np.max(time)-np.min(time)
#    sampleRate = len(time)/timeRange            # samples per JD
#    period = sampleRate/2.2                     # 2.2 used as a buffer, per Harris and Lupishko (Asteroids II)
    period = timeRange*5.0
    return period                               # in JD


''' Fit the data to a model using scipy.optimize.leastsq.
    Takes the dataset(a LightCurveData object), the range of orders to check (a two-element list)
    Three period guess cases:
        - no guess (Nyquist sampling criterion used with a 15 minute step size)
        - range of periods provided (three-element list of [minGuess, maxGuess, step])
        - single initial guess provided (three-element list of [min, max, guess]) where hardMinPeriod and hardMaxPeriod are hard limits
    The range of orders to check over may also be provided
'''
def fitData(LightCurveData, fitOptions, method = None, periodGuess = None, hardMinPeriod = None, hardMaxPeriod = None):
    if fitOptions['timer']:
        startTime = clock()
    if method is None:         # Nyquist sampling criterion is used at an interval of 15 minutes
         maxRecoverablePeriod = findRecoverablePeriod(LightCurveData.data['jd'])*24.0       # converted to hours
         initPeriodList = floatRange(0.25, maxRecoverablePeriod, 0.25)
         print 'Checking up to the maximum recoverable period (0.25 hours to ' + str(max(initPeriodList)) + ' hours)'
    if method == 'single':
        if type(periodGuess) is not float and type(periodGuess) is not int:
            print 'Error: \'single\' period guess method requires a guess input of type int or float'
            return
        else:
            initPeriodList = [float(periodGuess)]
    elif method == 'range':
        if type(periodGuess) is not list:
            print 'Error: \'range\' period guess method requires a guess input of type [min, max, stepSize] in floats'
            return
        elif len(periodGuess) != 3:
            print 'Error: \'range\' period guess method requires a guess input of type [min, max, stepSize] in floats'
            return            
        else:
            initPeriodList = floatRange(periodGuess[0], periodGuess[1], periodGuess[2])
    elif method is not None:
        print 'Error: invalid period guess method'
        return

    # convert periods from hours to JD
    initPeriodList = np.array(initPeriodList)/24.0
    if hardMaxPeriod is not None and (type(hardMaxPeriod) is float or type(hardMaxPeriod) is int):
        hardMaxPeriod /= 24.0
    if hardMinPeriod is not None and (type(hardMinPeriod) is float or type(hardMinPeriod) is int):
        hardMinPeriod /= 24.0
    
    # information from light curve for plotting
    time = LightCurveData.data['jd']-LightCurveData.data['jd'][0]       # use the minimum JD as a reference point of 0
    mag = LightCurveData.data['diffMag']
    err = LightCurveData.data['magErr']

    # storing the best fitting result
    [bestFit, bestParams, bestVar, bestOrder, bestPeriod] = [None, None, None, None, None]

    # store the periods and errors
    periodsTested = []
    periodErrors = []
    
    for p in range(len(initPeriodList)):
        initPeriod = initPeriodList[p]
        print 'attempting P = ' + str(initPeriod*24.0) + ' hours'
        for m in range(fitOptions['minOrder'],fitOptions['maxOrder']+1):

            params = Parameters()   # create the Parameters object and add P, A, and B values to it
            params.add('P', value = initPeriod, min = hardMinPeriod, max = hardMaxPeriod)
            params.add('y', value = 0)
            for i in range(m):
                params.add('A' + str(i+1), value = 0)
                params.add('B' + str(i+1), value = 0)
            
            fitResult = minimize(makeModel, params, args = (time, mag, err))
            periodsTested.append(params['P'].value*24.0)       # add the period to a list for the error plot
            residuals = makeModel(fitResult.params, time, mag)
            periodErrors.append(np.sqrt(sum(residuals**2)/float(len(residuals))))

            if fitResult.success:
                # check amplitude
                modelTime = np.linspace(time.min(), time.min()+fitResult.params['P'].value,100)
                modelMags = makeModel(fitResult.params, modelTime)
                amp = max(modelMags)-min(modelMags)
                if amp < 2.0:
                    n = float(len(time))       # number of observations
                    k = float(2*m+1)           # total free parameters
                    var = 1/(n-k)*sum(makeModel(params, time, mag, err)**2)
                    
                    if bestFit is None:
                        [bestFit, bestParams, bestVar, bestOrder, bestPeriod] = [fitResult, params, var, m, initPeriod]
                    elif bestVar > var:
                        [bestFit, bestParams, bestVar, bestOrder, bestPeriod] = [fitResult, params, var, m, initPeriod]
                else:
                    print 'model rejected because amplitude > 2.0 (amplitude = ' + str(amp) + ')'
            else:
                print 'optimization failed for P = ' + str(initPeriod) + ', m = ' + str(m)
    if fitOptions['timer']:
        print 'fitting runtime = ' + str(clock()) + ' s for ' + str(len(initPeriodList)) + ' initial period guesses'
    return (bestFit, bestOrder, periodsTested, periodErrors)


''' Phase folds the magnitude data given a set that has the phase offset at the first data point
'''
def phaseFold(time, period):
    foldedTime = []
    for i in range(len(time)):
        periodOffset = int(time[i]/period)
        if periodOffset >= 1:
            foldedTime.append(time[i]-period*periodOffset)
        else:
            foldedTime.append(time[i])

    print 'fitted period: ' + str(period*24.0) + ' h'
    return np.array(foldedTime)

''' Calculates amplitude with uncertainties using the uncertainties package
'''
def ampUncertainties(params, time):
    uncParams = {}                          # convert parameters into uncertainty values
    nominalValues = []
    stdDev = []
    for i in range(len(params.keys())):
        key = params.keys()[i]
        uncParams[key] = ufloat(params[key].value,params[key].stderr)
    
    H = makeModelUncertainties(uncParams, time)
    for i in range(len(H)):
        nominalValues.append(H[i].nominal_value)
        stdDev.append(H[i].std_dev)
    maxIndex = nominalValues.index(max(nominalValues))
    minIndex = nominalValues.index(min(nominalValues))
    return H[maxIndex]-H[minIndex]


''' Takes the best fit output and the data and plots them
    plotting error bars is optional (default is True)
    plotting the full phase is optional (default is True)
'''
def outputResults(fit, m, LightCurveData, outputOptions, periodErrors = None):
    # used for saving figures
    global basepath, objectName
    # TODO - have to change this for the MANOS file structure
    filepath = os.path.abspath(os.path.join(basepath, "..", 'Data', LightCurveData.name))
        
    time = LightCurveData.data['jd']-LightCurveData.data['jd'][0]    
    mag = LightCurveData.data['diffMag']
    err = LightCurveData.data['magErr']

    fittedPeriod = fit.params['P'].value    
    foldedTime = phaseFold(time, fittedPeriod)      # phase folded time

    timeRange = np.ptp(foldedTime)
    if fittedPeriod > timeRange*1.25:         # check if the fitted period is significantly longer than the data
        print 'Warning: fitted period is ' + str(int(((fittedPeriod/timeRange)-1)*100)) + '% longer than the range of the dataset'

    if outputOptions['phaseFoldData']:                               # shift plot for phase folded time
        time = foldedTime    
    
    if outputOptions['plotFullPeriod'] and time.min()+fittedPeriod < time.max():
        modelTime = np.linspace(time.min(), time.max(),100)
    elif outputOptions['plotFullPeriod']:
        modelTime = np.linspace(time.min(), time.min()+fittedPeriod,100)
    else:
        modelTime = np.linspace(time.min(), time.max(),100)

    plt.figure()
    if outputOptions['plotResiduals']:
        plt.subplot(211)
    modelMags = makeModel(fit.params, modelTime)
    plt.plot(modelTime, modelMags, 'b-')     # plot the model fit

    amp = ampUncertainties(fit.params,modelTime)    # get the amplitude and amplitude uncertainty

    print 'amplitude = ' + str(amp)
    if (outputOptions['plotErrorbars']):                                                    # plot the data with error bars (default)
        plt.errorbar(time, mag, err, fmt = 'rx')
    else:                                                               # plot the data without error plots
        plt.plot(time, mag, 'rx')

    if outputOptions['phaseFoldData']:
        plt.xlabel('delta JD')
    else:
        plt.xlabel('JD + ' + str(int(min(LightCurveData.data['jd']))))    # the subtraction simplifies display of the JD
    plt.ylabel('Differential Magnitude')
    plt.gca().invert_yaxis()        # flip the y axis

    periodUnc = ufloat(fit.params['P'].value,fit.params['P'].stderr)*24.0

    # shows up to 6 decimal places in the period uncertainty
    plt.title(objectName + ', P = ' + str(periodUnc) + ' h, a = ' + str(amp) + ', m = ' + str(m))
    lightCurveAxis = plt.axis()     # used to make sure that the residuals plots the x limits the same way
    
    if outputOptions['plotResiduals']:
        plt.subplot(212)
        residuals = makeModel(fit.params, time, mag)
        plt.plot(time, residuals, 'rx')
        plt.plot([0, max(time)], [0, 0], 'k--')      # plot a zero line
        plt.title('Residuals for ' + objectName + ' Fit')
        plt.ylabel('Residual Magnitude')
        plt.xlim(lightCurveAxis[0:2])
        if outputOptions['phaseFoldData']:
            plt.xlabel('delta JD')
        else:
            plt.xlabel('JD + ' + str(int(min(LightCurveData.data['jd']))))

    plt.subplots_adjust(hspace = 0.5)
    plt.savefig(os.path.join(filepath,(LightCurveData.name + 'LightCurve')))               # save the light curve plot
    print 'RMS of residuals = ' + str(np.sqrt(sum(residuals**2)/float(len(residuals))))

    if periodErrors is not None and outputOptions['plotPeriodErrors']:
        if type(periodErrors) is not list:
            print 'Error: period errors must be given as a list'
        else:
            if len(periodErrors) != 2:
                print 'Error: invalid period errors format'
            else:
                plt.figure()
                plt.plot(periodErrors[0], periodErrors[1], 'x')
                plt.title('RMS of Residuals for ' + objectName + ' Fit')
                plt.ylabel('RMS (magnitude)')
                plt.xlabel('period (h)')
                plt.savefig(os.path.join(filepath, (LightCurveData.name + 'MeanResiduals')))               # save the light curve plot

    if outputOptions['printReport']:
        print '\nFIT RESULTS:'
        report_fit(fit.params, show_correl = False)
    if outputOptions['saveReport']:
        f = open(os.path.join(filepath, (LightCurveData.name + 'FitReport.txt')), 'w')
        f.write(fit_report(fit.params,show_correl=False))
        f.write('\nmodel amplitude = ' + str(amp))
        f.close()

    if outputOptions['showPlots']:
        plt.show()


''' Given an object name, this will look in the appropriate folder and get all the data files, as well as the fitInfo specifications file
    The combined information is returned as the information necessary to create a LightCurveData object, fit the data, and print the results
'''
def extractRunOptions(objectName):
    filePath, objectFiles = lookInFolder('file', objectName)
    fileNamesAndFormat = {}
    offsets, method, T0, hardMinPeriod, hardMaxPeriod = None, None, None, None, None
    fitInfoFound = False    
    expectedNumFiles = 0
    readingOffsets = None
    
    for o in range(len(objectFiles)):
        if 'standard' in objectFiles[o]:        # all standardized data files will have 'standard' in their names
            with open(os.path.join(filePath,objectFiles[o])) as f:      # check if we're looking at a 3-col or 4-col file
                content = f.readlines()
                if len(string.split(content[0])) == 3:
                    fileNamesAndFormat[os.path.join(objectName,objectFiles[o])] = [['jd',0],['diffMag',1],['magErr',2]]
                elif len(string.split(content[0])) == 4:
                    fileNamesAndFormat[os.path.join(objectName,objectFiles[o])] = [['jd',0],['diffMag',1],['magErr',2],['night',3]]
                elif len(string.split(content[0])) == 8:    # this should be the actual MANOS format
                    fileNamesAndFormat[os.path.join(objectName,objectFiles[o])] = [['jd',3],['diffMag',6],['magErr',7]]
        elif objectFiles[o] == objectName + '_fitInfo.txt':
            if fitInfoFound == False:
                fitInfoFound = True
                with open(os.path.join(filePath,objectFiles[o])) as f:
                    content = f.readlines()
                    for i in range(len(content)):
                        lineWords = string.split(content[i])

                        # check for ignore cases
                        if lineWords is None:           # ignore lines that are blank 
                            continue
                        elif len(lineWords) == 0:       # ignore lines that are blank
                            continue
                        elif lineWords[0] is '\n' or lineWords[0][0] is '#':        # ignore lines that are newlines or have a '#' at the beginning
                            continue

                        # cases with useful information
                        else:
                            if readingOffsets is not None:
                                if 'ENDOFFSETS' in lineWords[0]:
                                    readingOffsets = None
                                else:
                                    if len(lineWords) != 2:
                                        print 'Error: invalid offset specification'
                                    else:
                                        offsets[readingOffsets][int(lineWords[0])] = float(lineWords[1])
                                
                            elif lineWords[0] == 'FILES':     # checking the number of files
                                expectedNumFiles = int(lineWords[1])        # the fitInfo file has a specified number of files
                            elif lineWords[0] == 'GUESS':     # reading in the period guess value(s)
                                numWords = len(lineWords)   # check that the number of parameters is right (if the line exists, it should be of length 3, 5, or 7
                                if numWords != 3 and numWords != 5 and numWords != 7:
                                    print 'Error: invalid guess format'
                                elif numWords == 3:         # the only case where there are three words on the 'guess' is if method == 'single'
                                    if lineWords[1] != 'single':
                                        print 'Error: invalid guess format'
                                    else:
                                        method = lineWords[1]
                                        T0 = float(lineWords[2])
                                elif numWords == 5:         # both the 'single' and 'range' cases can have 5 words, depending on whether hard mins and maxes are there
                                    if lineWords[1] == 'single':
                                        method = lineWords[1]
                                        T0 = lineWords[2]
                                        minPeriod = float(lineWords[3])
                                        maxPeriod = float(lineWords[4])
                                    elif lineWords[1] == 'range':
                                        method = lineWords[1]
                                        T0 = [float(lineWords[2]),float(lineWords[3]),float(lineWords[4])]
                                    else:
                                        print 'Error: invalid guess method'
                                elif numWords == 7:         # only 'range' case can have 7 words in the 'guess' line
                                    if lineWords[1] == 'range':
                                        method = lineWords[1]
                                        T0 = [float(lineWords[2]),float(lineWords[3]),float(lineWords[4])]
                                        minPeriod = float(lineWords[5])
                                        maxPeriod = float(lineWords[6])
                                    else:
                                        print 'Error: invalid guess format'
                            elif lineWords[0] == 'HARDMINPERIOD':           # the hard minimum period (no guess will be attempted for values below this)
                                if len(lineWords) != 2:
                                    print 'Error: invalid specification of a hard minimum guess'
                                else:
                                    hardMinPeriod = float(lineWords[1])
                            elif lineWords[0] == 'HARDMAXPERIOD':           # the hard maximum period (no guess will be attempted for values above this)
                                if len(lineWords) != 2:
                                    print 'Error: invalid specification of a hard maximum guess'
                                else:
                                    hardMaxPeriod = float(lineWords[1])
                            elif lineWords[0] == 'OFFSET':
                                if offsets is None:
                                    offsets = {}
                                offsets[lineWords[1]] = {}          # create a new sub-dictionary for a file's offset
                                readingOffsets = lineWords[1]
                                
            else:       # duplicate fitInfo file
                print 'Error: duplicate fitInfo file found'

    # turn the offset dictionary into a list, because dictionaries are not guarenteed to order keys in the same order
    offsetsList = []
    if offsets is not None:
        for i in range(len(fileNamesAndFormat.keys())):
            if fileNamesAndFormat.keys()[i] in offsets.keys():
                offsetsList.append(offsets[fileNamesAndFormat.keys()[i]])
            else:
                print 'Error: missing offset data for ' + fileNamesAndFormat.keys()[i]

    # error checking
    if len(fileNamesAndFormat) == 0:
        print 'Error: no data files found'
    if fitInfoFound == False:
        print 'Error: fitInfo file not found'
    else:
        if expectedNumFiles != len(fileNamesAndFormat):
            print 'Error: the expected and actual number of data files do not match'

    return fileNamesAndFormat, offsetsList, method, T0, hardMinPeriod, hardMaxPeriod
    

'''         start the command line script        '''
RunOptionsShell().cmdloop()
