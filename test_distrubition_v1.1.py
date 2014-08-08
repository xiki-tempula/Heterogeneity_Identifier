import sys
import os
import numpy as np
import prettyplotlib as ppl
import matplotlib.pyplot as plt
from scipy import optimize
import copy

def convertoriginal(amp, dwell):
		originaldwell = dwell * 9 #resolution is 111ms
		originaldwell = np.rint(originaldwell)
		originaldata = np.array([])
		for i in range(len(originaldwell)):
				originaldata = np.append(originaldata, np.ones(originaldwell[i]) * amp[i])
		return originaldata

def originalplot(savefilename, amp, dwell):
		originaldata = convertoriginal(amp, dwell)
		print 'length of the original data:', len(origindata)
		ppl.plot(originaldata, color='k')
		plt.title(savefilename + '_original')
		plt.savefig(pathfilename + '_original' + '.png', dpi = 500)
		print 'saved original figure', savefilename
		plt.close()

def convert(dwell):
	dwell *= 9 #resolution is 111ms
	dwell = np.rint(dwell)
	result = np.array([])
	for idwell in range(len(dwell)):
			if idwell % 2 == 0:
					result = np.append(result, np.ones(dwell[idwell]))
			else:
					result = np.append(result, np.zeros(dwell[idwell]))
	return result

def generaterandom(length, change, first, second):
		length = int(length)
		change = int(change)
		first = float(first)
		second = float(second)
		result = np.zeros(length)
		ranresult = np.random.random_sample((length,))
		for iran in range(length):
				if (iran< change) & (ranresult[iran] < first):
						result[iran] = 1.0
				elif (iran >= change) & (ranresult[iran] < second):
						result[iran] = 1.0
		return result

def movavg(result, window):
		#window must be an odd number
		if window % 2 == 0:
				window += 1
				#print 'Interval for moving average should be an odd number. Window += 1'
		weights = np.ones(window) / window
		ma = np.convolve(weights, result, mode='valid')
		zeros = np.zeros(np.rint((window-1)/2))
		ma = np.hstack((zeros, ma, zeros))
		#This step makes sure that the moving average sequence is as long as the original one.
		'''count = np.nonzero(ma)[0]
		print 'Moving average:', 'start', count[0], 'end', count[-1], 'vaild length', len(count), 'total length', len(ma)'''
		return ma

def movstd(ma, window):
		#window must be an odd number
		if window % 2 == 0:
				window += 1
				#print 'Interval for moving standard deviation should be an odd number. window += 1'
		std = np.zeros(len(ma))
		for index in range(window - 1, len(ma) - window + 1):
				std[index] = np.std(ma[index - (window-1)/2: index + (window-1)/2 + 1])
		#Same as the moving average. Ensuring the length of moving standard deviation is the same as the original one.
		'''count = np.nonzero(std)[0]
		print 'std:', 'start', count[0], 'end', count[-1], 'vaild length', len(count), 'total length', len(std)'''
		return std

def movavgstdtriangle(result):
		depth = np.rint(len(result)/2)+1
		#if depth >3000:
		#		depth = 3000
		for i in np.arange(3,depth,2):
				sma = movavg(result,i)
				plotcolour = np.array([])
				temp = movstd(sma,i)
				plotcolour = np.vstack((plotcolour,temp))
				print savefilename, 'processing', i, 'of', depth, 'completed'
		return plotcolour

def findmaxstd(std, window, precentage):
		#window must be an odd number
		if window % 2 == 0:
				window += 1
				#print 'Interval for finding max for moving standard deviation should be an odd number. window += 1'
		findmaxstdzeros = np.zeros(window - 1)
		findmaxstddiff = np.diff(std[np.nonzero(std)[0]])
		findmaxstddiffwhole = np.hstack((findmaxstdzeros, findmaxstddiff, findmaxstdzeros))
		fig, ax = plt.subplots(1)
		ppl.plot(np.arange(len(findmaxstddiffwhole)), findmaxstddiffwhole)
		ppl.plot(np.arange(len(findmaxstddiffwhole)), [0] * len(findmaxstddiffwhole))

		fig.savefig(pathfilename + '_window_' + str(interval) + '_diff_' + '.png',dpi=500)
		plt.close()


		#One index number is lost due to differentiation
		#Clean the zeros and add them again to prevent the huge change between
		#the last zero and the first standard deviation

		peak = np.array([])
		stdamp = np.array([])
		for check in range(window - 1, len(findmaxstddiffwhole) - (window - 1)):
				firsthalf = np.all(findmaxstddiffwhole[check - (window-1)*precentage: check] > 0)
				secondhalf = np.all(findmaxstddiffwhole[check: check + (window - 1)*precentage +1] <= 0)
				#Compensate for the loss of one number due to the differentiation
				if (firsthalf & secondhalf):
						if std[check] != 0:
								peak = np.append(peak, check)
								stdamp = np.append(stdamp, std[check])
		return peak, stdamp

def fsolvepopendiffstdamp(popendiff, interval, stdamp):
		calstdamp = np.arange(1, (interval+1)/2)
		calstdamp = (popendiff/2/((interval-1)/2) * calstdamp) ** 2
		solve = stdamp - np.sqrt(np.mean(calstdamp))
		return solve

def popendiffpeak(result, peak):
		popendiff = np.array([])
		calpeak = np.hstack((0, peak, len(result) - 1))
		for index in np.arange(1, len(calpeak) - 1):
				temp1 = np.mean(result[calpeak[index]:calpeak[index + 1]])
				temp2 = np.mean(result[calpeak[index - 1]:calpeak[index]])
				temp = abs(temp1 - temp2)
				popendiff = np.append(popendiff, temp)
		return popendiff

def popenpeak(result, peak):
		popen = np.array([])
		calpeak = np.hstack((0, peak, len(result) - 1))
		for index in range(len(calpeak)-1):
				temp = np.mean(result[calpeak[index]:calpeak[index + 1]])
				popen = np.append(popen, temp)
		return popen

def popendiffstdamp(stdamp, popendiff, interval):
		rootpopendiffstdamp = np.array([])
		for index in range(len(stdamp)):
				rootstdamp = stdamp[index]
				rootpopendiff = popendiff[index]
				temp = optimize.fsolve(fsolvepopendiffstdamp, rootpopendiff, args=(interval, rootstdamp))
				rootpopendiffstdamp = np.append(rootpopendiffstdamp, temp)
		return rootpopendiffstdamp

def funstdleastsq(x, coeffs):
		#need to specify interval before using this function
		time = coeffs[0]
		testamp = coeffs[1]
		caldatafirsthalf = [testamp/2 * -1] * abs(interval - (x - time)-1)
		caldatasecondhalf = [testamp/2] * abs(interval + (x - time)-1)
		caldata = np.hstack((caldatafirsthalf, 0 , caldatasecondhalf))
		weights = np.ones(interval) / interval
		calsma = np.convolve(weights, caldata, mode='valid')
		calstd = calsma.std()
		return calstd

def funstd(x, time, popendiff):
		#need to specify GLOBAL interval before using this function
		caldatafirsthalf = [popendiff/2 * -1] * abs(interval - (x - time)-1)
		caldatasecondhalf = [popendiff/2] * abs(interval + (x - time)-1)
		caldata = np.hstack((caldatafirsthalf, 0 , caldatasecondhalf))
		weights = np.ones(interval) / interval
		calma = np.convolve(weights, caldata, mode='valid')
		calstd = np.std(calma)
		return calstd

def cruvefit(std, window, peak, stdamp, popendiff):
		fitpeak = np.array([])
		fitpopendiff = np.array([])
		fitstdamp = np.array([])
		fitwindow = np.array([])
		failedattemp = 'interval:' + str(window)
		failedattempcount = 0
		for tpeak in range(len(peak)):
				xdata = np.arange(peak[tpeak]-(window-1)*precentage, peak[tpeak]+(window-1)*precentage+1)
				ydata = std[peak[tpeak]-(window-1)*precentage: peak[tpeak]+(window-1)*precentage+1]
				guess = [peak[tpeak], popendiff[tpeak]]
				params, params_covariance = optimize.curve_fit(funstd, xdata, ydata, guess)
				if funstd(params[0], params[0], params[1]) < (0.5*std[params[0]]):
						print 'Curve fit based on Popen solved from standard deviation peak failed.'
						print 'Peak time', peak[tpeak], 'Retry for maximum of ten times.'
						count = 1
						while count < 11:
								params, params_covariance = optimize.curve_fit(funstd,xdata,ydata,guess)
								if funstd(params[0], params[0], params[1]) < (0.5*std[params[0]]):
										count += 1
								else:
										print 'Retrying:', str(count), 'of 10','success'
										count += 100

						if count == 11:
								print peak[tpeak], 'Curve fitting failed', 'Using initial guess instead'
								failedattemp += ' ' + str(peak[tpeak])
								failedattempcount += 1
								params[0] = peak[tpeak]
								params[1] = popendiff[tpeak]
				fitpeak = np.append(fitpeak, params[0])
				fitpopendiff = np.append(fitpopendiff, params[1])
				fitstdamp = np.append(fitstdamp, funstd(params[0], params[0], params[1]))
		failedattemp += ' ' + str(failedattempcount) + ' out of ' + str(len(fitpeak)) + 'failed'
		return fitpeak, fitpopendiff, fitstdamp, failedattemp

def residuals(coeffs, y, x):
		return y - funstdleastsq(x, coeffs)

def drawpeak(peak, popendiff, step):
		xfit = np.array([])
		yfit = np.array([])
		for idraw in range(len(peak)):
				for isample in np.arange(peak[idraw]-(step-1), peak[idraw] + step):
						xfit = np.append(xfit, isample)
						yfit = np.append(yfit, funstd(isample, peak[idraw], popendiff[idraw]))
		return xfit, yfit

def drawpeaksingle(peak, popendiff, step, precentage):
		xfit = np.array([])
		yfit = np.array([])
		for isample in np.arange(peak-(step-1)*precentage, peak + (step-1)*precentage +1):
				xfit = np.append(xfit, isample)
				yfit = np.append(yfit, funstd(isample, peak, popendiff))
		return xfit, yfit

def orifitdiff(std, peak, popendiff, interval, precentage):
		xfit, yfit = drawpeaksingle(peak, popendiff, interval, precentage)
		xfit = xfit.astype(np.int64)
		calorifitdiff = np.mean((std[xfit]-yfit) ** 2)
		calorifitdiff /= std[peak]
		return calorifitdiff

def plotorifit(filename, std, peak, popendiff, interval, precentage):
		xfit, yfit = drawpeak(peak, popendiff, interval)
		xdiff = [0]
		ydiff = [0]
		for index in range(len(peak)):
				xdiff.append(peak[index]-1)
				ydiff.append(0)
				xdiff.append(peak[index])
				ydiff.append(orifitdiff(std, peak[index], popendiff[index], interval, precentage))
				xdiff.append(peak[index]+1)
				ydiff.append(0)
		xdiff.append(len(result)-1)
		ydiff.append(0)
		ydiff = np.array(ydiff) * (np.max(std)/np.max(ydiff)/2)
		if np.max(ydiff) == 0:
				print ydiff
				print np.max(ydiff)
				sys.exit(0)
		fig, ax = plt.subplots(1)
		ppl.plot(ax, xfit, yfit, label='curve_fit', linewidth=1.)
		ppl.plot(ax, np.arange(len(std)), std, label='original', linewidth=0.75)
		ppl.plot(ax, xdiff, ydiff, label='fit_diff', linewidth=1)
		ppl.legend(ax)
		ax.set_title(filename + ' window ' + str(interval) + ' original std and fit')
		fig.savefig(pathfilename + '_window_' + str(interval) + '_original_std_fit' + '.png',dpi=500)
		plt.close()
		print filename, 'Standard deviation and curve fit plot saved.'

def plotpopendiff(filename, result, peak, popendiffstdamp):
		fig, ax = plt.subplots(1)
		plotpeak = np.hstack((peak.astype(np.int64), len(result)))[::-1]
		popen = popenpeak(result, peak)[::-1]
		popenoriginal = np.empty(len(result))
		for index in range(len(plotpeak)):
				popenoriginal[:plotpeak[index]] = popen[index]

		plotpeak = plotpeak[::-1]
		popenfit = copy.copy(popenoriginal)

		for index in range(len(popendiffstdamp)-1):
				popenfit[plotpeak[index]-1] = (popenoriginal[plotpeak[index]-1] + popenoriginal[plotpeak[index]])/2 - popendiffstdamp[index]/2
				popenfit[plotpeak[index]] = (popenoriginal[plotpeak[index]-1] + popenoriginal[plotpeak[index]])/2 + popendiffstdamp[index]/2


		ppl.plot(ax, np.arange(len(popenfit)), popenfit, label='calculated', linewidth=1)
		ppl.plot(ax, np.arange(len(popenoriginal)), popenoriginal, label='original', linewidth=2)
		if np.all(popenoriginal[len(popenoriginal)/2:] > (max(popenoriginal)/2)):
				ppl.legend(ax, loc='lower right')
		elif np.all(popenoriginal[len(popenoriginal)/2:] < (max(popenoriginal)/2)):
				ppl.legend(ax, loc='upper right')
		elif np.all(popenoriginal[:len(popenoriginal)/2] < (max(popenoriginal)/2)):
				ppl.legend(ax, loc='upper left')
		elif np.all(popenoriginal[:len(popenoriginal)/2] > (max(popenoriginal)/2)):
				ppl.legend(ax, loc='lower left')

		ax.set_title(savefilename + ' window ' + str(interval) + ' original Popen and fit')
		fig.savefig(pathfilename + '_window_' + str(interval) + '_original_Popen_fit' + '.png',dpi=500)
		plt.close()
		print savefilename, 'Popen and calculated Popen plot saved.'

def determinse(result, precentage):
		attempt = np.zeros(20)
		count = 0
		interval = 101
		firstpeak = np.array([])
		lastpeak = np.array([])
		while count < 10:
				interval += 3
				print 'testing interval:', interval,
				ma = movavg(result, interval)
				std = movstd(ma, interval)
				oripeak, oristdamp = findmaxstd(std, interval, precentage)
				peaklength = np.append(oripeak, len(result)-1) - np.append(0, oripeak)

				if np.amin(peaklength) > interval:
						count += 1
						print 'Succeed.', str(count), 'out of 10.'
						attempt = np.append(attempt, 1)
				else:
						count = 0
						print 'failed.'
						attempt = np.append(attempt, 0)

				judge = np.sum(attempt[-15:])
				print 'Success rate:', str(judge), 'out of 15.'
				if judge > 11:
						print 'successful rate exceeded 0.8.'
						count = 100

				firstpeak = np.append(firstpeak, oripeak[0])
				lastpeak = np.append(lastpeak, oripeak[-1])

		if count == 10:
				start = interval - 27
				first = np.amin(firstpeak[-10:])
				last = int((len(result) - 1 - np.amax(lastpeak[-10:]))/2)
				end = np.amin(np.append(first, last))
		elif count == 100:
				startindex = -(15 - np.where(attempt[-15:] == 1)[0][0])
				start = interval + startindex * 3 +3
				first = np.amin(firstpeak[startindex:])
				last = int((len(result) - 1 - np.amax(lastpeak[startindex:]))/2)
				end = np.amin(np.append(first, last))

		if end > 800:
				end = 800
		print 'Interval: starting point', str(start), 'ending point', str(end)

		if (end - start) / 3 < 20:
				if end < 500:
						print 'WARNING: Some points may be lost because it is too close to the start or the end.'
				if start > 500:
						print 'WARNING: Model change is too frequent. The analysis may not be very accurate.'
				print 'The end of test will be manually set as', str(start + 300)
				step = np.arange(start, start + 303, 6)
		else:
				step = np.linspace(start, end, 50)
		return step.astype(np.int64)

sysinput = sys.argv

inputfilename = lambda x: x[-4:] == '.csv'
filenamelist = filter(inputfilename, sysinput)

inputfilename = lambda x: x == 'all'
if filter(inputfilename, sysinput) != []:
		filenamelist = ['04071424.csv', '030714c4.csv', '010714c5.csv', '180407c2.csv', '270307c3.csv', '210307c6.csv', '280307c6.csv', '020507c4.csv', '020507c6.csv', '230307c2.csv']

if filenamelist == []:
		print 'No file input detected. Please type in a XXX.csv file.'
		sys.exit(0)
else:
		for filename in filenamelist:
				if os.path.exists(filename) == False:
						print filename, 'does not exist!'
						sys.exit(0)
		print 'The list of files which will be processed', filenamelist

saveoriginal = lambda x: x == 'saveoriginal'
saveoriginalplot = filter(saveoriginal, sysinput)


adjustprecentage = lambda x: x[-1] == '%'
precent = filter(adjustprecentage, sysinput)

if precent != []:
		precentage = float(precent[0][:-1]) / 100
		print 'precentage input detected. precentage = ', precentage
else:
		precentage = 0.45
		print 'precentage input not detected. precentage = 0.45'

for filename in filenamelist:
		start, end, amp, dwell = np.loadtxt(filename, delimiter=',',usecols=(4,5,6,8),unpack=True)
		savefilename = filename[:-4] + '_start_' + "%.3f" % (start[0]/1000) + '_end_' + "%.3f" % (end[-1]/1000)
		print 'Begin processing', savefilename

		if os.path.isdir(os.getcwd() +'/' + savefilename) == False:
				os.mkdir(os.getcwd() +'/' + savefilename)
		pathfilename = os.getcwd() +'/' + savefilename + '/' + savefilename

		if saveoriginalplot != []:
				originalplot(savefilename, amp, dwell)

		result = convert(dwell)
		#Convert the result to a series of ones and zeros
		print 'Start finding the start point.'

		if os.path.exists(pathfilename+'step.npy') == False:
				step = determinse(result, precentage)
				np.save(pathfilename+'step.npy', step)
		else:
				print 'Existing calculation for interval found. Loading'
				step = np.load(pathfilename+'step.npy')
				step = step.astype(np.int64)
				print 'Interval: starting point', str(step[0]), 'ending point', str(step[-1])

		totalpeak = np.array([])
		totalinterval = np.array([])
		count = 1
		for interval in step:
				print 'interval:', interval
				ma = movavg(result, interval)
				#Calculate the moving average of the result
				std = movstd(ma, interval)
				#Calculate the moving standard deviation of the result
				oripeak, oristdamp = findmaxstd(std, interval, precentage)
				#Obtain the standard deviation peak and the amplitude of the peak from first order differentiation

				if oripeak == np.array([]):
						print 'No peak found'
				else:
						oripopendiffpeak = popendiffpeak(result, oripeak)
						#Calculate the difference in Popen from the raw data using the time of the peak
						oripopendiffstdamp = popendiffstdamp(oristdamp, oripopendiffpeak, interval)
						#Calculate the difference in Popen based on the value of the standard deviation peak
						#fitpeak, fitpopendiffstdamp, fitstdamp, failedattemp = cruvefit(std, interval, oripeak, oristdamp, oripopendiffstdamp)
						#failedcurvefit = np.hstack((failedcurvefit, failedattemp))
						#Using curve fit to calculate peak, Popen difference and amplitude of standard deviation
						#from the original peak, original Popen difference solved from the original amplitude of standard deviation
						#and the amplitude of standard deviation calculated from the Popen difference obtained from curve fitting
						#Peaks failed to be fitted are collected in failedcurvefit
						#fitpopendiffpeak = popendiffpeak(result, fitpeak)
						#Calculate the Popen difference form the peak obtained from curve fit
						#Just to calculate the Popen difference from a more accurate peak time
						#The data obtained seems to be the same

						while np.any(oripopendiffstdamp > (1.5 * oripopendiffpeak)):
								filterpeak = np.array([])
								filterpopendiffpeak = np.array([])
								for index in range(len(oripeak)):
										if oripopendiffstdamp[index] < (1.5 * oripopendiffpeak[index]):
												filterpeak = np.append(filterpeak, oripeak[index])
												filterpopendiffpeak = np.append(filterpopendiffpeak, oripopendiffstdamp[index])
								oripeak = filterpeak
								oripopendiffstdamp = filterpopendiffpeak
								oripopendiffpeak = popendiffpeak(result, oripeak)

						totalpeak = np.append(totalpeak, oripeak)
						totalinterval = np.append(totalinterval, [interval] * len(oripeak))


						for iprint in range(len(oripeak)):
								print 'peak', oripeak[iprint], 'Popen(avg)', oripopendiffpeak[iprint], 'Popen(std)', oripopendiffstdamp[iprint]

						if count == 10:
								#plotorifit(savefilename, std, oripeak, oripopendiffstdamp, interval, precentage)
								plotpopendiff(savefilename, result, oripeak, oripopendiffstdamp)
								count = 0
						count += 1

		fig, ax = plt.subplots(1)
		ppl.scatter(ax, totalpeak, totalinterval)
		ax.set_title(savefilename + ' Heterogeneity distrubition Start:')
		fig.savefig(pathfilename + '_Heterogeneity_distrubition' + '.png',dpi=500)
		plt.close()
		print 'Saved scatter plot:' + savefilename





'''
		if os.path.exists(savefilename + '.npy'):
				print savefilename + '.npy', 'Found. Loading the data.'
				plotcolour = np.load(savefilename + '.npy')
				print savefilename + '.npy', 'Loaded.'
		else:
				result = convert(dwell)
				plotcolour = movavgstdtriangle(result)
				np.save(savefilename, plotcolour)
				print savefilename + '.npy', 'saved'

		ppl.pcolormesh(plotcolour)
		plt.title(savefilename)
		plt.savefig(pathfilename + '.png',dpi=500)
		plt.close()
		print savefilename + '.png', 'saved'
'''
