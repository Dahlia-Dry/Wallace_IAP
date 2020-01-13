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

    def create_masterbias(self):
        biases = []
        for filename in os.listdir(self.directory):
            print(filename)
            if filename.endswith('.fit'):
                biases.append(CCDData.read(self.directory + filename, unit= 'adu'))
        print(len(biases))
        c = Combiner(biases)
        masterbias = c.median_combine()
        masterbias.write(self.directory + '/masterbias.fits')

    def plot(self, exptimes, constraints=None): #all the files for given telescope on given date
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
                        if float(fits.open(self.irectory +filename, ignore_missing_end=True)[0].header['EXPTIME']) == float(key) and \
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


ccd = CCDStats('Data/20200110_p4/')
ccd.plot({'10.0':[],'30.0':[],'90.0':[],'180.0':[], '540.0':[]}, {'SET-TEMP':-30})
#create_masterbias('20200108_p3')
