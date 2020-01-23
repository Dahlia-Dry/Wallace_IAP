import os
from scipy import stats
from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy.io import fits
import ccdproc as ccdp
from ccdproc import Combiner
import matplotlib.pyplot as plt
import numpy as np


class CCDStats(object):
    def __init__(self, directory):
        self.directory = directory

    def create_master(self, frame):
        if frame == "Bias Frame":
            biases = []
            for filename in os.listdir(self.directory):
                print(filename)
                if filename.endswith('.fit') and fits.open(self.directory + filename, ignore_missing_end=True)[0].header['IMAGETYP'] == frame:
                    biases.append(CCDData.read(self.directory + filename, unit= 'adu'))
            print(len(biases))
            c = Combiner(biases)
            masterbias = c.median_combine()
            masterbias.write(self.directory + 'masterbias.fits')
            print('masterbias written to ' + self.directory)
        elif frame == "Dark Frame":
            darks = []
            for filename in os.listdir(self.directory):
                print(filename)
                if filename.endswith('.fit') and fits.open(self.directory + filename, ignore_missing_end=True)[0].header['IMAGETYP'] == frame:
                    darks.append(CCDData.read(self.directory + filename, unit='adu'))
            print(len(darks))
            c = Combiner(darks)
            masterdark = c.median_combine()
            masterdark.write(self.directory + 'masterdark.fits')
            print('masterdark written to ' + self.directory)
        else:
            raise TypeError('Frame should be Dark Frame or Bias Frame')

    def combine(self, data, mode):
        frames = []
        for file in data:
            if file.endswith('.fit'):
                frames.append(CCDData.read(self.directory + file, unit='adu'))
        c = Combiner(frames)
        if mode == 'mean':
            master = c.average_combine()
        elif mode == 'median':
            master= c.median_combine()
        else:
            raise ValueError('mode must be mean or median')
        master = np.asarray(master)
        master = np.subtract(master, 100)  # PEDESTAL
        return master

    def plot(self, xvar, yvar, frame, constraint=None, compile=None, fit=1, adu='mean'):
        i = 0
        colors = ['r','o','g']
        if constraint is not None:
            cdict = {}
            for filename in os.listdir(self.directory):
                if not filename.endswith('.fit') or fits.open(self.directory+ filename,ignore_missing_end=True)[0].header['IMAGETYP'] != frame:
                    pass
                else:
                    if float(fits.open(self.directory +filename, ignore_missing_end=True)[0].header[constraint]) not in cdict:
                        cdict[float(fits.open(self.directory +filename, ignore_missing_end=True)[0].header[constraint])] = [filename]
                    else:
                        cdict[float(fits.open(self.directory + filename, ignore_missing_end=True)[0].header[constraint])].append(filename)
        for key in cdict.keys():
            xdict = {}
            for filename in cdict[key]:
                if not filename.endswith('.fit') or fits.open(self.directory + filename, ignore_missing_end=True)[0].header['IMAGETYP'] != frame:
                        break
                else:
                    if float(fits.open(self.directory + filename, ignore_missing_end=True)[0].header[xvar]) not in xdict:
                        xdict[float(fits.open(self.directory + filename, ignore_missing_end=True)[0].header[xvar])] = [filename]
                    else:
                        xdict[float(fits.open(self.directory + filename, ignore_missing_end=True)[0].header[xvar])].append(filename)
            cdict[key] = xdict
            if compile is not None:
                for x in xdict.keys():
                    data = xdict[x]
                    if yvar == 'pixel count':
                        yarr = self.combine(data, compile)
                        if adu == 'mean':
                            y = np.mean(yarr)
                        elif adu =='median':
                            y = np.median(yarr)
                        else:
                            raise ValueError('adu must be mean or median')
                    else:
                        values = []
                        for file in data:
                            values.append(float(fits.open(self.directory + file, ignore_missing_end=True)[0].header[yvar]))
                        if compile == 'mean':
                            y = np.mean(values)
                        elif compile == 'median':
                            y = np.median(values)
                        else:
                            raise ValueError('compile must be mean or median if not None')
                    xdict[x] = [y]
            else:
                for x in xdict.keys():
                    yvalues = []
                    data = xdict[x]
                    for file in data:
                        if yvar == 'pixel count':
                            yarr = CCDData.read(self.directory + file, unit= 'adu')
                            yarr= np.asarray(yarr)
                            yarr = np.subtract(yarr, 100)  # PEDESTAL
                            if adu == 'mean':
                                y = np.mean(yarr)
                            elif adu =='median':
                                y = np.median(yarr)
                            else:
                                raise ValueError('adu must be mean or median')
                        else:
                            y = float(fits.open(self.directory + file, ignore_missing_end=True)[0].header[yvar])
                        yvalues.append(y)
                    xdict[x] = yvalues
            x = []
            y = []
            for val in xdict:
                y += xdict[val]
                for yval in xdict[val]:
                    x.append(val)
            plt.scatter(x, y, label=key)
            i = i + 1
            slope, intercept, r, p, std_err = stats.linregress(x, y)
            coefs = [slope, intercept]
            plt.plot(np.linspace(0, max(x)), [coefs[0] * i + coefs[1] for i in np.linspace(min(x), max(x))])

        title = xvar + ' vs ' + yvar + ' for ' + self.directory[-3:-1]
        legend = plt.legend(loc='lower right', shadow=True)
        plt.title(title)
        plt.xlabel(xvar)
        plt.ylabel(yvar)

        plt.xlim(0, np.max(np.array(x)) + 10)
        plt.ylim(coefs[1] - 100, np.max(np.array(y)) + 100)
        figname = self.directory[-3:-1] + '_' + self.directory[5:13]
        plt.savefig(figname + '.png')
        plt.show()

    def plot_exp(self, exptimes, constraints=None): #all the files for given telescope on given date
        master_darks = []
        if constraints is not None:
            for filename in os.listdir(self.directory):
                if fits.open(self.directory+ filename,ignore_missing_end=True)[0].header['IMAGETYP'] == 'Dark Frame':
                    for key in list(exptimes.keys()):
                        score = 0
                        for constraint in list(constraints.keys()):
                            if constraint == "SET-TEMP":
                                if float(fits.open(self.directory +filename, ignore_missing_end=True)[0].header[constraint]) > float(constraints[constraint])-2.5 and \
                                        float(fits.open(self.directory + filename, ignore_missing_end=True)[0].header[constraint]) < float(constraints[constraint])+2.5:
                                    pass
                                else:
                                    score = score + 1
                            else:
                                try:
                                    if float(fits.open(self.directory +filename, ignore_missing_end=True)[0].header[constraint]) != float(constraints[constraint]):
                                        score = score +1
                                except ValueError:
                                    if str(fits.open(self.directory +filename, ignore_missing_end=True)[0].header[constraint]) != str(constraints[constraint]):
                                        score = score +1
                        if score == 0:
                            if filename.endswith('.fit'):
                                if float(fits.open(self.directory +filename, ignore_missing_end=True)[0].header['EXPTIME']) == float(key):
                                    exptimes[key].append(filename)
        else:
            for filename in os.listdir(self.directory):
                for key in list(exptimes.keys()):
                    if filename.endswith('.fit'):
                        if float(fits.open(self.directory +filename, ignore_missing_end=True)[0].header['EXPTIME']) == float(key) and \
                                fits.open(self.directory+ filename,ignore_missing_end=True)[0].header['IMAGETYP'] == 'Dark Frame':
                            exptimes[key].append(filename)
        print(exptimes)
        for key in list(exptimes.keys()):
            print(key)
            darks = []
            for file in exptimes[key]:
                print(file, key)
                if file.endswith('.fit'):
                    darks.append(CCDData.read(self.directory + file, unit= 'adu'))
            c = Combiner(darks)
            masterdark = c.median_combine()
            masterdark = np.asarray(masterdark)
            masterdark = np.subtract(masterdark, 100) #PEDESTAL
            master_darks.append(masterdark)
        print(len(master_darks))
        y = []
        for master in master_darks:
            y.append(np.mean(master))
        x = list(exptimes.keys())
        x = [float(i) for i in x]
        plt.scatter(x,y)
        slope, intercept, r, p, std_err = stats.linregress(x, y)
        coefs = [slope, intercept]
        plt.plot(np.linspace(0, max(x)), [coefs[0]*i + coefs[1] for i in np.linspace(min(x), max(x))], color = 'r')
        plt.text(50,coefs[1]-50,'y = ' +str(coefs[0]) + '*x + ' + str(coefs[1]))
        plt.text(50, coefs[1]-75, 'r = ' + str(r) + ', r^2 = ' + str(r**2))
        title = 'Mean Pixel Count vs Exposure Time for ' + self.directory[-3:-1] + ' at '
        for constraint in list(constraints.keys()):
            title = title + str(constraint) + '=' + str(constraints[constraint]) + ' '
        plt.title(title)
        plt.xlabel('Exposure Time (s)')
        plt.ylabel('Pixel Count (ADU)')

        plt.xlim(0, np.max(np.array(x))+10)
        plt.ylim(coefs[1] - 100, np.max(np.array(y))+100)
        figname = self.directory[-3:-1] + '_' + self.directory[5:13] + '_'
        for constraint in list(constraints.keys()):
            figname = figname + str(constraint) + '_' + str(constraints[constraint]) + '_'
        plt.savefig(figname + '.png')
        plt.show()

    def tempvtime(self):
        x = []
        y = []
        for filename in os.listdir(self.directory):
            x.append(float(filename[-8:-4]))
            y.append(float(fits.open(self.directory + filename, ignore_missing_end=True)[0].header['CCD-TEMP']))
        plt.scatter(x, y)
        plt.xlabel('time')
        plt.ylabel('ccd-temp')
        plt.title('CCD Temp v Time ' + self.directory[5:-1])
        plt.savefig('CCD Temp v Time ' + self.directory[5:-1] + '.png')
        plt.show()


ccd = CCDStats('Data/20200117_p4/')
#ccd.plot_exp({'10.0':[],'30.0':[],'90.0':[],'180.0':[], '540.0':[]}, {'SET-TEMP':-30})
#ccd.create_master('Bias Frame')
ccd.plot('EXPTIME', 'pixel count', 'Dark Frame', constraint='CCD-TEMP', compile = 'mean')
