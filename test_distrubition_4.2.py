import sys
import os
import numpy as np
import prettyplotlib as ppl
import matplotlib.pyplot as plt
from scipy import optimize, signal
import copy
import multiprocessing
import matplotlib as mpl


def originalplot(pathfilename, savefilename, amp, dwell):
	originaldata = convertoriginal(amp, dwell)
	print 'length of the original data:', len(originaldata)/9, 'ms'
	ppl.plot(originaldata, color='k')
	plt.title(savefilename + '_original')
	plt.savefig(pathfilename + '_original.png', dpi = 300)
	print 'saved original figure', savefilename
	plt.close()

def convertoriginal(amp, dwell):
	originaldwell = dwell * 9 #resolution is 111ms
	originaldwell = np.rint(originaldwell)
	originaldata = np.array([])
	for i, dwell in enumerate(originaldwell):
		originaldata = np.append(originaldata, np.ones(dwell) * amp[i])
	return originaldata

'''def convert(dwell):
	dwell *= 9 #resolution is 111ms
	dwell = np.rint(dwell)
	result = np.array([])
	for index, idwell in enumerate(dwell):
		if index % 2 == 0:
			result = np.append(result, np.ones(idwell))
		else:
			result = np.append(result, np.zeros(idwell))
	return result'''

'''def convertma(dwell, interval, jump):
	cumsum = np.cumsum(dwell)
	start = 0.0
	end = interval * 1.0
	ma = np.array([])
	while end <= cumsum[-1]:
		opentime = 0.0

		startindex = np.where(cumsum > start)[0][0]
		if startindex % 2 == 0:
			opentime += cumsum[startindex] - start

		endindex = np.where(cumsum > end)[0][0]
		if endindex % 2 == 0:
			opentime += end - cumsum[endindex - 1]

		for index in range(startindex +1, endindex):
			if index % 2 == 0:
				opentime += dwell[index]

		ma = np.append(ma, opentime / (end - start))
		start += jump
		end += jump
	ma = np.hstack((np.zeros((interval-jump)/jump/2), ma, np.zeros((interval-jump)/jump/2)))
	return ma'''

def convertma(dwell, interval, jump):
	cumsum = np.cumsum(dwell)
	start = 0.0
	end = interval * 1.0
	ma = np.array([])
	while end <= cumsum[-1]:
		opentime = 0.0

		startindex = np.where(cumsum > start)[0][0]
		if startindex % 2 == 0:
			opentime += cumsum[startindex] - start

		endindex = np.where(cumsum > end)[0][0]
		if endindex % 2 == 0:
			opentime += end - cumsum[endindex - 1]

		for index in range(startindex +1, endindex):
			if index % 2 == 0:
				opentime += dwell[index]

		ma = np.append(ma, opentime / (end-start)*1.0)

		start += jump
		end += jump
	ma = np.hstack((np.zeros((interval-jump)/jump/2), ma, np.zeros((interval-jump)/jump/2)))
	return ma

def convertma(dwell, interval, jump):
	cumsum = np.cumsum(dwell)
	start = 0.0
	end = interval * 1.0
	ma = np.array([])
	while end <= cumsum[-1]:
		popen = calculatepopen(dwell, start/jump, end/jump, jump)
		ma = np.append(ma, popen)

		start += jump
		end += jump
	ma = np.hstack((np.zeros((interval-jump)/jump/2), ma, np.zeros((interval-jump)/jump/2)))
	return ma

def determine_step_process(dwell, percentage, start, end, jump):
	attempt = np.array([])
	firstpeak = np.array([])
	lastpeak = np.array([])
	suminterval = np.array([])
	for window in np.arange(start, (end+jump*2), jump*2):
		#print 'testing interval:', interval,
		#window must be an odd number
		if window/jump % 2 == 0:
			window += jump

		ma = convertma(dwell, window, jump)
		interval = window / jump
		std = movstd(ma, interval)
		oripeak = signal.find_peaks_cwt(std, np.arange(interval * percentage, interval, interval*0.05))
		if oripeak != []:
			peaklength = np.append(oripeak, np.cumsum(dwell)[-1]/jump-1) - np.append(0, oripeak)
			if np.amin(peaklength) > interval:
				attempt = np.append(attempt, 1)
			else:
				attempt = np.append(attempt, 0)
			firstpeak = np.append(firstpeak, oripeak[0])
			lastpeak = np.append(lastpeak, oripeak[-1])
			suminterval = np.append(suminterval, interval)

	np.savez(pathfilename+'_'+str(start)+'_'+str(end), attempt = attempt, firstpeak = firstpeak, lastpeak = lastpeak, suminterval = suminterval)

'''def movavg(result, window):
	weights = np.ones(window) / window
	ma = np.convolve(weights, result, mode='valid')
	zeros = np.zeros((window-1)/2)
	ma = np.hstack((zeros, ma, zeros))
	#This step makes sure that the moving average sequence is as long as the original one.
	return ma'''

def findmax_std(std, window, percentage):
	#window must be an odd number
	if window % 2 == 0:
		window += 1
		#print 'Interval for finding max for moving standard deviation should be an odd number. window += 1'
	findmaxstd_zeros = np.zeros(window - 1)
	findmaxstd_diff = np.diff(std[np.nonzero(std)[0]])
	findmaxstd_whole = np.hstack((findmaxstd_zeros, findmaxstd_diff, findmaxstd_zeros))
	#fig, ax = plt.subplots(1)
	#ppl.plot(np.arange(len(findmaxstddiffwhole)), findmaxstddiffwhole)
	#ppl.plot(np.arange(len(findmaxstddiffwhole)), [0] * len(findmaxstddiffwhole))
	#fig.savefig(pathfilename + '_window_' + str(window) + '_diff_' + '.png',dpi=500)
	#plt.close()

	#One index number is lost due to differentiation
	#Clean the zeros and add them again to prevent the huge change between
	#the last zero and the first standard deviation

	peak = np.array([])
	stdamp = np.array([])
	for check in np.arange(window - 1, len(findmaxstd_whole) - (window - 1)):
		firsthalf = np.all(findmaxstd_whole[check - (window-1)*percentage: check] > 0)
		secondhalf = np.all(findmaxstd_whole[check: check + (window - 1)*percentage +1] <= 0)
		#Compensate for the loss of one number due to the differentiation
		if (firsthalf & secondhalf):
			if std[check] != 0:
			#I don't know why I need to do this but this can prevent some weird situations.
				peak = np.append(peak, check)
				stdamp = np.append(stdamp, std[check])
	peak = peak.astype(int)
	return peak, stdamp

def movstd(ma, interval):
	std = np.zeros(len(ma))
	for index in np.arange(interval - 1, len(ma) - (interval - 1)):
		std[index] = np.std(ma[index - (interval-1)/2: index + (interval-1)/2 + 1])
	#Same as the moving average. Ensuring the length of moving standard deviation is the same as the original one.
	return std

def calculatepeak(step, jump, dwell):
	totalstd = []
	totalpeak = np.array([])
	totalinterval = np.array([])
	totalfilterpeak = np.array([])
	totalfilterinterval = np.array([])
	#totalfilterpopen = np.array([])
	#totalfilterpopeninterval = np.array([])
	for window in step:
		if window/jump % 2 == 0:
			window += jump
		ma = convertma(dwell, window, jump)
		interval = window / jump
		#Calculate the moving average of the result
		std = movstd(ma, interval)
		#Calculate the moving standard deviation of the result
		oripeak = signal.find_peaks_cwt(std, np.arange(interval * percentage, interval, interval*0.05))
		oristdamp = std[oripeak]
		#oripeak, oristdamp = findmax_std(std, interval, percentage)
		#Obtain the standard deviation peak and the amplitude of the peak

		if len(oripeak) == 0:
			pass
		else:
			if totalstd == []:
				totalstd = std
			else:
				totalstd = np.vstack((totalstd, std))

			oripopendiffstdamp = oristdamp / np.sqrt(1.0/12)
			#Calculate the difference in Popen based on the value of the standard deviation peak
			#filterpopen, oripopendiffstdamp = filterpopendiff(oripeak, oripopendiffstdamp, popendifflimit)
			oripopendiffpeak = calpopendiffpeak(dwell, oripeak)
			#Calculate the difference in Popen from the raw data using the time of the peak
			filterpeak, oripopendiffstdamp, oripopendiffpeak = filterpeaktime(dwell, std, oripeak, oripopendiffstdamp, oripopendiffpeak)

			totalpeak = np.append(totalpeak, oripeak)
			totalinterval = np.append(totalinterval, [interval] * len(oripeak))
			totalfilterpeak = np.append(totalfilterpeak, filterpeak)
			totalfilterinterval = np.append(totalfilterinterval, [interval] * len(filterpeak))
			#totalfilterpopen = np.append(totalfilterpopen, filterpopen)
			#totalfilterpopeninterval = np.append(totalfilterpopeninterval, [interval] * len(filterpopen))
			print 'window', window
			for iprint in range(len(filterpeak)):
				print 'peak', filterpeak[iprint] * jump, 'Popen(avg)', oripopendiffpeak[iprint], 'Popen(std)', oripopendiffstdamp[iprint]

	np.savez(pathfilename+str(step[0])+'_'+str(step[-1]),
	totalpeak=totalpeak,
	totalinterval=totalinterval,
	totalfilterpeak=totalfilterpeak,
	totalfilterinterval=totalfilterinterval,
	totalstd=totalstd)

'''def calculatepopen(dwell, start, end, jump):
	cumsum = np.cumsum(dwell)/jump
	opentime = 0.0

	startindex = np.where(cumsum > start)[0][0]
	if startindex % 2 == 0:
		opentime += cumsum[startindex] - start

	endindex = np.where(cumsum >= end)[0][0]
	if endindex % 2 == 0:
		opentime += end - cumsum[endindex - 1]

	for index in range(startindex +1, endindex):
		if index % 2 == 0:
			opentime += dwell[index]/jump

	popen = opentime/(end - start)
	return popen'''

def calculatepopen(dwell, start, end, jump):
	popendwell = np.append(0, dwell)/jump
	cumsum = np.cumsum(popendwell)
	opentime = 0.0

	startindex = np.where(cumsum > start)[0][0]
	starttime = cumsum[startindex - 1] - start
	if startindex % 2 != 0:
		opentime += starttime

	endindex = np.where(cumsum >= end)[0][0]
	endtime = end - cumsum[endindex - 1]
	if endindex % 2 != 0:
		opentime += endtime

	middle = (cumsum > start) & (cumsum < end)
	opentime += np.sum(popendwell[middle & (np.arange(len(popendwell))%2 != 0)])
	time = end - start

	popen = opentime/time

	return popen


def calpopendiffpeak(dwell, peak):
	popendiff = np.array([])
	calpeak = np.hstack((0, peak, np.cumsum(dwell)[-1]/jump))
	for index in np.arange(1, len(calpeak) - 1):
		temp1 = calculatepopen(dwell, calpeak[index], calpeak[index + 1], jump)
		temp2 = calculatepopen(dwell, calpeak[index - 1], calpeak[index], jump)
		temp = abs(temp1 - temp2)
		popendiff = np.append(popendiff, temp)
	return popendiff

def filterpeaktime(dwell, std, peak, popendiffstdamp, popendiffpeak):
	for i in [1.7]:
		delete = np.where(popendiffstdamp > (i * popendiffpeak))[0]
		while np.any(delete):
			peak = np.delete(peak, delete[0])
			popendiffpeak = calpopendiffpeak(dwell, peak)
			popendiffstdamp = std[peak] / np.sqrt(1.0/12)
			delete = np.where(popendiffstdamp > (i * popendiffpeak))[0]

	for i in [0.7]:
		delete = np.where(popendiffstdamp < (i * popendiffpeak))[0]
		while np.any(delete):
			peak = np.delete(peak, delete[0])
			popendiffpeak = calpopendiffpeak(dwell, peak)
			popendiffstdamp = std[peak] / np.sqrt(1.0/12)
			delete = np.where(popendiffstdamp < (i * popendiffpeak))[0]
			
	return peak, popendiffstdamp, popendiffpeak

def popenpeak(dwell, peak):
	popen = np.array([])
	calpeak = np.hstack((0, peak, np.cumsum(dwell)[-1]/jump))
	for index in range(len(calpeak)-1):
		temp = calculatepopen(dwell, calpeak[index], calpeak[index + 1], jump)
		popen = np.append(popen, temp)
	return popen

def filterpopendiff(peak, dwell, popendifflimit):
	popen = popenpeak(dwell, peak)
	popendiff = np.diff(popen)
	while np.any(abs(popendiff) < popendifflimit):
		delete = np.nan
		for index in range(len(peak)):
			if abs(popendiff[index]) < popendifflimit:
				delete = index

		peak = np.delete(peak, delete)
		popen = popenpeak(dwell, peak)
		popendiff = np.diff(popen)

	return peak

'''def prepare_hist(X, tres):
	n=len(X)
	if n <= 300:
		nbdec = 5.0
	elif n <= 1000:
		nbdec = 8.0
	elif n <= 3000:
		nbdec = 10.0
	else:
		nbdec = 12.0
	logx = np.log(X)
	bins = np.arange(np.log(tres), np.amax(logx) + np.log(10)/nbdec, np.log(10)/nbdec)
	hist, bin_edges = np.histogram(logx, bins = bins)
	plothist = np.hstack((0, np.repeat(hist, 2), 0)).astype(float)
	plotbin_edges = np.repeat(np.exp(bin_edges), 2)
	return plotbin_edges, plothist'''

def plothist(plotdwellopen, plotdwellclose, begin, close):
	n=len(plotdwellopen)
	if n <= 300:
		nbdec = 5.0
	elif n <= 1000:
		nbdec = 8.0
	elif n <= 3000:
		nbdec = 10.0
	else:
		nbdec = 12.0
	logx = np.log(plotdwellopen)
	if (np.amax(logx) + np.log(10)/nbdec) > np.log(0.1):
		bins = np.arange(np.log(0.1), np.amax(logx) + np.log(10)/nbdec, np.log(10)/nbdec)
		plt.hist(logx, bins = bins)
		plt.title(savefilename + ' ' + str(begin) + '-' + str(close) + ' open histgram')
		plt.savefig(pathfilename + '_' + str(begin) + '-' + str(close) + '_open_histgram' + '.png',dpi=300)
		plt.close()
		print savefilename, str(begin) + '-' + str(close), 'Open Histogram saved.'

	n=len(plotdwellclose)
	if n <= 300:
		nbdec = 5.0
	elif n <= 1000:
		nbdec = 8.0
	elif n <= 3000:
		nbdec = 10.0
	else:
		nbdec = 12.0
	logx = np.log(plotdwellclose)
	if (np.amax(logx) + np.log(10)/nbdec) > np.log(0.1):
		bins = np.arange(np.log(0.1), np.amax(logx) + np.log(10)/nbdec, np.log(10)/nbdec)
		plt.hist(logx, bins = bins)
		plt.title(savefilename + ' ' + str(begin) + '-' + str(close) + ' close histgram')
		plt.savefig(pathfilename + '_' + str(begin) + '-' + str(close) + '_close_histgram' + '.png',dpi=300)
		plt.close()
		print savefilename, str(begin) + '-' + str(close), 'Close Histogram saved.'






sysinput = sys.argv

inputfilename = lambda x: x[-4:] == '.csv'
filenamelist = filter(inputfilename, sysinput)

if 'all' in sys.argv:
	filenamelist = os.listdir(os.curdir)
	filenamelist = filter(inputfilename, filenamelist)

if filenamelist == []:
	print 'No file input detected. Please type in a XXX.csv file.'
	sys.exit(0)
else:
	sysexit = 0
	for filename in filenamelist:
		if os.path.exists(filename) == False:
			print filename, 'does not exist!'
			sysexit = 1
	if sysexit == 1:
		sys.exit(0)
	else:
		print 'The list of files which will be processed:', filenamelist

adjustpercentage = lambda x: x[-1] == '%'
precent = filter(adjustpercentage, sysinput)



if precent != []:
	percentage = float(precent[0][:-1]) / 100
	print 'percentage input detected. percentage = ', percentage
else:
	percentage = 0.45
	print 'percentage input not detected. percentage = 0.45'

for filename in filenamelist:
	start, end, amp, dwell = np.loadtxt(filename, delimiter=',',usecols=(4,5,6,8),unpack=True)
	savefilename = filename[:-4] + '_s_'+ "%.3f" % (start[0]/1000) + '_e_' + "%.3f" % (end[-1]/1000)
	print 'Begin processing', savefilename

	if os.path.isdir(os.getcwd() + '/' + savefilename) == False:
		os.mkdir(os.getcwd() + '/' + savefilename)
	pathfilename = os.getcwd() + '/' + savefilename + '/' + savefilename


	#originalplot(pathfilename, savefilename, amp, dwell)
	#plot the original data

	print 'Start finding the start point.'

	jump = 0.25
	start = 10.25
	end = 50.25
	popendifflimit = 0.05
	if os.path.exists(pathfilename+'step.npy') == False:

		processlist = np.array_split(np.arange(start,end,jump*2), multiprocessing.cpu_count())
		if __name__ == '__main__':
			jobs = []
			for time in processlist[1:]:
				p = multiprocessing.Process(target=determine_step_process, args=(dwell, percentage, time[0], time[-1], jump))
				jobs.append(p)
				p.start()

		window = start - jump*2
		attempt = np.zeros(15)
		count = 0
		firstpeak = np.array([])
		lastpeak = np.array([])
		suminterval = np.array([])
		while (count < 10) & (window in np.arange(processlist[0][0]-jump*2, processlist[0][-1]+jump*2, jump*2)):
			window += jump*2
			print 'testing interval:', window, 'ms'
			ma = convertma(dwell, window, jump)

			interval = window / jump
			std = movstd(ma, interval)
			oripeak = signal.find_peaks_cwt(std, np.arange(interval * percentage, interval, interval*0.05))
			if oripeak != []:
				peaklength = np.append(oripeak, np.cumsum(dwell)[-1]/jump-1) - np.append(0, oripeak)
				if np.amin(peaklength) > interval:
					count += 1
					print 'Succeed.', str(count), 'out of 10.'
					attempt = np.append(attempt, 1)
				else:
					count = 0
					print 'Failed.'
					attempt = np.append(attempt, 0)

				judge = np.sum(attempt[-15:])
				print 'Success rate:', str(judge), 'out of 15.'
				if judge > 11:
					print 'Success rate exceeded 0.8.'
					count = 100

				firstpeak = np.append(firstpeak, oripeak[0])
				lastpeak = np.append(lastpeak, oripeak[-1])
				suminterval = np.append(suminterval, interval)
			else:
				print 'Peak not found in this interval.'


		if count > 9:
			for job in jobs:
				if job.is_alive() == True:
					job.terminate()

		if count == 10:
			start = interval - 18
			first = np.amin(firstpeak[-10:])

			last = int((np.cumsum(dwell)[-1]/jump - 1 - np.amax(lastpeak[-10:]))/2)
			end = np.amin(np.append(first, last))
		elif count == 100:
			startindex = -(15 - np.where(attempt[-15:] == 1)[0][0])
			start = interval + startindex * 2 + 2
			first = np.amin(firstpeak[startindex:])
			last = int((np.cumsum(dwell)[-1]/jump - 1 - np.amax(lastpeak[startindex:]))/2)
			end = np.amin(np.append(first, last))
		else:
			attempt = attempt[15:]
			for job in jobs:
				job.join()

			for time in processlist[1:]:
				loaddata = pathfilename+'_'+str(time[0])+'_'+str(time[-1])+'.npz'
				attempt = np.append(attempt, np.load(loaddata)['attempt'])
				suminterval = np.append(suminterval, np.load(loaddata)['suminterval'])
				firstpeak = np.append(firstpeak, np.load(loaddata)['firstpeak'])
				lastpeak = np.append(lastpeak, np.load(loaddata)['lastpeak'])
				os.remove(loaddata)
			weights = np.ones(10) / 10
			filterten = np.convolve(weights, attempt, mode='valid')
			weights = np.ones(15) / 15
			filterfifteen = np.convolve(weights, attempt, mode='valid')
			if (1 in filterten) or (0.8 in filterfifteen):
				ten = np.where(filterten == 1)[0][0]
				fifteen = np.where(filterfifteen == 0.8)[0][0]
				startindex = np.where(attempt[fifteen: fifteen+15] == 1)[0][0] + fifteen
				endindex = np.where(attempt[fifteen: fifteen+15] == 1)[0][-1] + fifteen

				if suminterval[ten] <= suminterval[startindex]:
					start = suminterval[ten]
					first = np.amin(firstpeak[ten: ten+10])
					last = int((np.cumsum(dwell)[-1]/jump - 1 - np.amax(lastpeak[ten: ten+10]))/2)
					end = np.amin(np.append(first, last))
				else:
					start = suminterval[startindex]
					first = np.amin(firstpeak[startindex: endindex])
					last = int((np.cumsum(dwell)[-1]/jump - 1 - np.amax(lastpeak[startindex: endindex]))/2)
					end = np.amin(np.append(first, last))
			else:
				print 'WARNING: The model change is too frequent. Setting starting point to', end
				start = end
				firstpeak = np.array([])
				lastpeak = np.array([])
				for window in np.arange(50.25,80.25,0.25):
					print 'Calculating last peak from window', window, 'ms'
					ma = convertma(dwell, window, jump)
					interval = window / jump
					std = movstd(ma, interval)
					oripeak = signal.find_peaks_cwt(std, np.arange(interval * percentage, interval, interval*0.05))
					firstpeak = np.append(firstpeak, oripeak[0])
					lastpeak = np.append(lastpeak, oripeak[-1])
				first = np.amin(firstpeak)
				last = int((np.cumsum(dwell)[-1]/jump - 1 - np.amax(lastpeak))/2)
				end = np.amin(np.append(first, last))
				if end*jump < 62.25:
					print 'WARNING: The first or the last model change is too close to the edge.'
					print 'Setting ending point to 801.'
					end = 80/jump

		if end*jump > 80.25:
			end = 80.25/jump

		if (end - start)/2 < 50:
			end = start + 300
			print 'WARNING: The first or the last model change is too close to the edge.'
			print 'Setting ending point to', str(end*jump)

		print 'Interval: starting point', str(start*jump), 'ending point', str(end*jump)
		start *= jump
		end *= jump
		step = np.linspace(start, end, 50)
		step = np.rint(step*4) / 4
		np.save(pathfilename+'step.npy', step)
	else:
		print 'Existing calculation for interval found. Loading'
		step = np.load(pathfilename+'step.npy')
		print 'Interval: starting point', str(step[0]), 'ending point', str(step[-1])

	processlist = np.array_split(step, multiprocessing.cpu_count())

	if __name__ == '__main__':
		jobs = []
		for time in processlist:
			time = time.astype(int)
			p = multiprocessing.Process(target=calculatepeak, args=(time, jump, dwell))

			jobs.append(p)
			p.start()

		for job in jobs:
			job.join()

	sumpeak = np.array([])
	suminterval = np.array([])
	sumfilterpeak = np.array([])
	sumfilterinterval = np.array([])
	totalstd = []
	#sumfilterpopen = np.array([])
	#sumfilterpopeninterval = np.array([])
	for time in processlist:
		time = time.astype(int)
		loaddata = pathfilename+str(time[0])+'_'+str(time[-1])+'.npz'
		sumpeak = np.append(sumpeak, np.load(loaddata)['totalpeak'])
		suminterval = np.append(suminterval, np.load(loaddata)['totalinterval'])
		sumfilterpeak = np.append(sumfilterpeak, np.load(loaddata)['totalfilterpeak'])
		sumfilterinterval = np.append(sumfilterinterval, np.load(loaddata)['totalfilterinterval'])
		#sumfilterpopen = np.append(sumfilterpopen, np.load(loaddata)['totalfilterpopen'])
		#sumfilterpopeninterval = np.append(sumfilterpopeninterval, np.load(loaddata)['totalfilterpopeninterval'])
		if totalstd == []:
			totalstd = np.load(loaddata)['totalstd']
		else:
			totalstd = np.vstack((totalstd, np.load(loaddata)['totalstd']))
		os.remove(loaddata)

	'''fig, ax = plt.subplots(1)
	ppl.scatter(ax, sumpeak, suminterval*jump)
	ppl.scatter(ax, sumfilterpeak, sumfilterinterval*jump)
	ax.set_xlim(0, len(totalstd[0]))
	ax.set_ylim(min(suminterval*jump), max(suminterval*jump))
	ax.set_title(savefilename + ' Heterogeneity distrubition')
	fig.savefig(pathfilename + '_Heterogeneity_distrubition' + '.png',dpi=300, transparent = True)
	plt.close()
	print 'Saved scatter plot:', savefilename'''


	plt.pcolormesh(np.arange(len(totalstd[0]))*jump, step, totalstd, cmap=mpl.cm.Reds)
	plt.scatter(sumpeak*jump, suminterval*jump, c= 'y')
	#plt.scatter(sumfilterpopen*jump, sumfilterpopeninterval*jump, c='b')
	plt.scatter(sumfilterpeak*jump, sumfilterinterval*jump, c ='g')

	plt.xlim(0, len(totalstd[0])*jump)
	plt.ylim(step[0], step[-1])
	plt.title(savefilename + ' Heterogeneity distrubition')
	plt.savefig(pathfilename + '_Heterogeneity_distrubition' + '.png',dpi=300)
	plt.close()
	print 'Saved scatter plot:', savefilename

	sortsumpeak = np.sort(sumfilterpeak)
	sortsumpeakdiff = np.diff(sortsumpeak)
	diffdistrubition = np.sort(sortsumpeakdiff.copy())

	limit = np.mean(diffdistrubition)

	fig, ax = plt.subplots(1)
	ppl.plot(np.arange(len(sortsumpeakdiff)), sortsumpeakdiff)
	ppl.plot(np.arange(len(diffdistrubition)), diffdistrubition)
	ppl.plot(np.arange(len(diffdistrubition)), np.ones(len(diffdistrubition))*limit)
	ax.set_title(savefilename + ' difference distrubition')
	fig.savefig(pathfilename + '_difference_distrubition' + '.png',dpi=300)
	plt.close()
	print 'Saved peak time difference plot:', savefilename

	csvpeak = np.array([])
	split = np.hstack((0, np.where(sortsumpeakdiff>limit)[0] + 1, len(sortsumpeak)))

	for index in np.arange(len(split)-1):
		if (split[index + 1] - split[index]) > 15:
			csvpeak = np.append(csvpeak, np.mean(sortsumpeak[split[index]: split[index + 1]]))
	csvpeak = csvpeak.astype(np.int64)

	csvpeak = filterpopendiff(csvpeak, dwell, popendifflimit)

	peak = csvpeak.copy()
	start, end, amp, dwell = np.loadtxt(filename, delimiter=',',usecols=(4,5,6,8),unpack=True)
	csvpeak = np.hstack((csvpeak.astype(np.int64), np.cumsum(dwell)[-1]/jump - 1))
	popen = popenpeak(dwell, peak)
	csvstart = np.hstack((0, peak.astype(float)))*jump + start[0]
	csvend = np.hstack((peak.astype(float), np.cumsum(dwell)[-1]/jump-1))*jump + start[0]
	meanopentime = np.array([])
	meanclosetime = np.array([])
	plotcsvpeak = np.append(0, csvpeak)*jump
	cumsum = np.cumsum(dwell)/jump
	if len(plotcsvpeak) > 1:
		for index in np.arange(1, len(plotcsvpeak)):
			rangehist = (cumsum>plotcsvpeak[index-1]) & (cumsum<plotcsvpeak[index])
			plothist(dwell[rangehist & (np.arange(len(dwell))%2==0)], dwell[rangehist & (np.arange(len(dwell))%2!=0)], plotcsvpeak[index-1], plotcsvpeak[index])
			meanopentime = np.append(meanopentime, np.mean(dwell[rangehist & (np.arange(len(dwell))%2==0)]))
			meanclosetime = np.append(meanclosetime, np.mean(dwell[rangehist & (np.arange(len(dwell))%2!=0)]))
	csvtext = np.transpose(np.vstack((csvpeak, csvstart, csvend, popen, meanopentime, meanclosetime)))
	np.savetxt(pathfilename+'.csv', csvtext, delimiter=',')

	plotcsvpeak = np.append(0, csvpeak)*jump
	cumsum = np.cumsum(dwell)/jump
	if len(plotcsvpeak) > 1:
		for index in np.arange(1, len(plotcsvpeak)):
			rangehist = (cumsum>plotcsvpeak[index-1]) & (cumsum<plotcsvpeak[index])
			plothist(dwell[rangehist & (np.arange(len(dwell))%2==0)], dwell[rangehist & (np.arange(len(dwell))%2!=0)], plotcsvpeak[index-1], plotcsvpeak[index])

	plotpeak = np.hstack((0, peak, peak + 1, np.cumsum(dwell)[-1]/jump-1))
	plotpeak = np.sort(plotpeak)*jump
	#plotpeak += start[0]
	plotpopen = np.repeat(popen, 2)

	plotoriamp = (np.repeat(amp, 2) - np.amin(amp)) / np.average(amp)/2
	plotorix = np.hstack((0, np.repeat(np.cumsum(dwell)[:-1], 2), np.cumsum(dwell)[-1]))
	plt.plot(plotorix, plotoriamp, c='black', linewidth = 0.5)
	plt.plot(plotpeak, plotpopen, c='red')
	plt.ylim(0, 1.2)
	plt.xlim(0, plotorix[-1])
	plt.title(savefilename + ' Popen')
	plt.savefig(pathfilename + '_Popen.png',dpi=300)
	plt.close()
	print savefilename, 'Popen plot saved.'
