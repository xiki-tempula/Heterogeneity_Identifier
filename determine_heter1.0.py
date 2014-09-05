import sys
import os
import numpy as np
import prettyplotlib as ppl
import matplotlib.pyplot as plt
from scipy import optimize, signal
import copy
import multiprocessing
import matplotlib as mpl

def calculatepopen(dwell, start, end, JUMP):
	popendwell = np.append(0, dwell)/JUMP
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

def fastcalpopen(arbdwell, popendwell, cumsum, start, end): #start and open need to be in arbitry unit
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

def genpopendiff(arbdwell, step):
	length = np.floor(np.sum(arbdwell))
	popendiff = np.zeros(length)
	popendwell = np.append(0, arbdwell)
	cumsum = np.cumsum(popendwell)
	for i in np.arange(step, length - step + 1):
		popendiff[i] = abs(
		fastcalpopen(arbdwell, popendwell, cumsum, i- step, i) -
		fastcalpopen(arbdwell, popendwell, cumsum, i, i + step))
	return popendiff

def filterpeaktime(dwell, peak, oripopendiff, JUMP):
	popendiffpeak = calpopendiffpeak(dwell, peak, JUMP)

	for i in np.arange(2.0, 1.2, -0.1):
		delete = np.where(oripopendiff[peak] > (i * popendiffpeak))[0]
		while np.any(delete):
			peak = np.delete(peak, delete[0])
			popendiffpeak = calpopendiffpeak(dwell, peak, JUMP)
			delete = np.where(oripopendiff[peak] > (i * popendiffpeak))[0]

	for i in np.arange(0.1, 0.8, 0.1):
		delete = np.where(oripopendiff[peak] < (i * popendiffpeak))[0]
		while np.any(delete):
			peak = np.delete(peak, delete[0])
			popendiffpeak = calpopendiffpeak(dwell, peak, JUMP)
			delete = np.where(oripopendiff[peak] < (i * popendiffpeak))[0]

	return peak

def calpopendiffpeak(dwell, peak, JUMP):
	calpopendiff = np.array([])
	calpeak = np.hstack((0, peak, np.sum(dwell)/JUMP))
	for index in np.arange(1, len(calpeak) - 1):
		temp = abs(
		calculatepopen(dwell, calpeak[index], calpeak[index + 1], JUMP) -
		calculatepopen(dwell, calpeak[index - 1], calpeak[index], JUMP))
		calpopendiff = np.append(calpopendiff, temp)
	return calpopendiff

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

	nbdec = 5.0 + (n > 300) * 3 + (n > 1000) * 2 + (n > 3000) * 2

	logx = np.log(plotdwellclose)
	if (np.amax(logx) + np.log(10)/nbdec) > np.log(0.1):
		bins = np.arange(np.log(0.1), np.amax(logx) + np.log(10)/nbdec, np.log(10)/nbdec)
		plt.hist(logx, bins = bins)
		plt.title(savefilename + ' ' + str(begin) + '-' + str(close) + ' close histgram')
		plt.savefig(pathfilename + '_' + str(begin) + '-' + str(close) + '_close_histgram' + '.png',dpi=300)
		plt.close()
		print savefilename, str(begin) + '-' + str(close), 'Close Histogram saved.'

def popenpeak(dwell, peak):
	popen = np.array([])
	calpeak = np.hstack((0, peak, np.cumsum(dwell)[-1]/JUMP))
	for index in range(len(calpeak)-1):
		temp = calculatepopen(dwell, calpeak[index], calpeak[index + 1], JUMP)
		popen = np.append(popen, temp)
	return popen

def filterpopendiff(peak, dwell, POPENDIFFLIMIT):
	popen = popenpeak(dwell, peak)
	popendiff = np.diff(popen)
	while np.any(abs(popendiff) < POPENDIFFLIMIT):
		delete = np.nan
		for index in range(len(peak)):
			if abs(popendiff[index]) < POPENDIFFLIMIT:
				delete = index

		peak = np.delete(peak, delete)
		popen = popenpeak(dwell, peak)
		popendiff = np.diff(popen)

	return peak

def determine_step_process(dwell, PERCENTAGE, start, end, JUMP):
	attempt = np.array([])
	firstpeak = np.array([])
	lastpeak = np.array([])
	suminterval = np.array([])
	arbdwell = dwell / JUMP
	for interval in np.arange(start, (end+JUMP), JUMP)/JUMP:
		#print 'testing interval:', interval,
		std = genpopendiff(arbdwell, interval)
		oripeak = signal.find_peaks_cwt(std, np.arange(interval * PERCENTAGE, interval, interval*0.05))
		if oripeak != []:
			peaklength = np.append(oripeak, np.sum(dwell)/JUMP-1) - np.append(0, oripeak)
			if np.amin(peaklength) > interval:
				attempt = np.append(attempt, 1)
			else:
				attempt = np.append(attempt, 0)
			firstpeak = np.append(firstpeak, oripeak[0])
			lastpeak = np.append(lastpeak, oripeak[-1])
			suminterval = np.append(suminterval, interval)

	np.savez(pathfilename+'_'+str(start)+'_'+str(end), attempt = attempt, firstpeak = firstpeak, lastpeak = lastpeak, suminterval = suminterval)


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
	PERCENTAGE = float(precent[0][:-1]) / 100
	print 'percentage input detected. percentage = ', PERCENTAGE
else:
	PERCENTAGE = 0.5
	print 'percentage input not detected. percentage = 0.5'

for filename in filenamelist:
	start, end, amp, dwell = np.loadtxt(filename, delimiter=',',usecols=(4,5,6,8),unpack=True)
	savefilename = filename[:-4] + '_s_'+ "%.3f" % (start[0]/1000) + '_e_' + "%.3f" % (end[-1]/1000)
	print 'Begin processing', savefilename

	makedirectory = os.path.join(os.getcwd(), savefilename)
	if os.path.isdir(makedirectory) == False:
		os.mkdir(makedirectory)
	pathfilename = os.path.join(makedirectory, savefilename)

	print 'Start finding the start point.'

	JUMP = 0.25
	start = 10.25
	end = 50.25
	POPENDIFFLIMIT = 0.05
	LOWLIMIT = 15
	arbdwell = dwell/JUMP

	if os.path.exists(pathfilename+'_step.npy') == False:
		processlist = np.array_split(np.arange(start,end,JUMP), multiprocessing.cpu_count())
		if __name__ == '__main__':
			jobs = []
			for time in processlist[1:]:
				p = multiprocessing.Process(target=determine_step_process, args=(dwell, PERCENTAGE, time[0], time[-1], JUMP))
				jobs.append(p)
				p.start()

		window = start - JUMP
		attempt = np.zeros(15)
		count = 0
		firstpeak = np.array([])
		lastpeak = np.array([])
		suminterval = np.array([])
		while (count < 10) & (window in np.arange(processlist[0][0]-JUMP, processlist[0][-1]+JUMP, JUMP)):
			window += JUMP
			print 'testing interval:', window, 'ms'
			interval = window / JUMP
			std = genpopendiff(arbdwell, interval)
			oripeak = signal.find_peaks_cwt(std, np.arange(interval * PERCENTAGE, interval, interval*0.05))
			if oripeak != []:
				peaklength = np.append(oripeak, np.cumsum(dwell)[-1]/JUMP-1) - np.append(0, oripeak)
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
				job.terminate()

		if count == 10:
			start = interval - 9
			first = np.amin(firstpeak[-10:])

			last = np.sum(dwell)/JUMP - 1 - np.amax(lastpeak[-10:])
			end = np.amin(np.append(first, last))

		elif count == 100:
			startindex = -(15 - np.where(attempt[-15:] == 1)[0][0])
			start = interval + startindex + 1
			first = np.amin(firstpeak[startindex:])
			last = np.sum(dwell)/JUMP - 1 - np.amax(lastpeak[startindex:])
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
					last = np.sum(dwell)/JUMP - 1 - np.amax(lastpeak[ten: ten+10])
					end = np.amin(np.append(first, last))
				else:
					start = suminterval[startindex]
					first = np.amin(firstpeak[startindex: endindex])
					last = sum(dwell)/JUMP - 1 - np.amax(lastpeak[startindex: endindex])
					end = np.amin(np.append(first, last))
			else:
				print 'WARNING: The model change is too frequent. Setting starting point to', end
				start = end
				firstpeak = np.array([])
				lastpeak = np.array([])
				for window in np.arange(50.25,80.25,0.25):
					print 'Calculating last peak from window', window, 'ms'
					ma = convertma(dwell, window, JUMP)
					interval = window / JUMP
					std = movstd(ma, interval)
					oripeak = signal.find_peaks_cwt(std, np.arange(interval * PERCENTAGE, interval, interval*0.05))
					firstpeak = np.append(firstpeak, oripeak[0])
					lastpeak = np.append(lastpeak, oripeak[-1])
				first = np.amin(firstpeak)
				last = np.sum(dwell)/JUMP - 1 - np.amax(lastpeak)
				end = np.amin(np.append(first, last))


		print 'Calculated starting point:', start*JUMP, 'ms'
		print 'Calculated ending point:', end*JUMP, 'ms'
		if end*JUMP > 80.25:
			end = 80.25/JUMP

		if (end - start) < 50:
			end = start + 50
			print 'WARNING: The first or the last model change is too close to the edge.'
			print 'Setting ending point to', str(end*JUMP)

		print 'Interval: starting point', str(start*JUMP), 'ending point', str(end*JUMP)
		start *= JUMP
		end *= JUMP
		step = np.linspace(start, end, 50)
		step = np.rint(step*4) / 4
		np.save(pathfilename+'_step.npy', step)
	else:
		print 'Existing calculation for interval found. Loading'
		step = np.load(pathfilename+'_step.npy')
		print 'Interval: starting point', str(step[0]), 'ending point', str(step[-1])

	pcolor = np.empty([len(step), np.floor(np.sum(arbdwell))])
	oripeak, oriinterval, filterpeak, filterinterval = ([] for i in range(4))

	for index, interval in enumerate(step/JUMP):
		print interval
		oripopendiff = genpopendiff(arbdwell, interval)
		pcolor[index,:] = oripopendiff
		peak = signal.find_peaks_cwt(oripopendiff, np.arange(interval*0.3, interval, interval*0.05))
		oripeak = np.append(oripeak, peak)
		oriinterval = np.append(oriinterval, [interval] * len(peak))
		peak = filterpeaktime(dwell, peak, oripopendiff, JUMP)
		filterpeak = np.append(filterpeak, peak)
		filterinterval = np.append(filterinterval, [interval] * len(peak))




	plt.pcolormesh(np.arange(len(pcolor[0]))*JUMP, step, pcolor, cmap=mpl.cm.Reds)
	plt.scatter(oripeak*JUMP, oriinterval*JUMP, c= 'y')
	plt.scatter(filterpeak*JUMP, filterinterval*JUMP, c ='g')

	plt.xlim(0, len(pcolor[0])*JUMP)
	plt.ylim(step[0], step[-1])
	plt.title(savefilename + ' Heterogeneity distrubition')
	plt.savefig(pathfilename + '_Heterogeneity_distrubition' + '.png',dpi=300)
	plt.close()
	print 'Saved scatter plot:', savefilename

	sumfilterpeak = filterpeak




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
		if (split[index + 1] - split[index]) > LOWLIMIT:
			csvpeak = np.append(csvpeak, np.mean(sortsumpeak[split[index]: split[index + 1]]))
	csvpeak = csvpeak.astype(np.int64)

	csvpeak = filterpopendiff(csvpeak, dwell, POPENDIFFLIMIT)

	peak = csvpeak.copy()
	start, end, amp, dwell = np.loadtxt(filename, delimiter=',',usecols=(4,5,6,8),unpack=True)
	csvpeak = np.hstack((csvpeak.astype(np.int64), np.cumsum(dwell)[-1]/JUMP - 1))
	popen = popenpeak(dwell, peak)
	csvstart = np.hstack((0, peak.astype(float)))*JUMP + start[0]
	csvend = np.hstack((peak.astype(float), np.cumsum(dwell)[-1]/JUMP-1))*JUMP + start[0]
	meanopentime = np.array([])
	meanclosetime = np.array([])

	plotcsvpeak = np.append(0, csvpeak)*JUMP
	cumsum = np.cumsum(dwell)/JUMP
	if len(plotcsvpeak) > 1:
		for index in np.arange(1, len(plotcsvpeak)):
			rangehist = (cumsum>plotcsvpeak[index-1]) & (cumsum<plotcsvpeak[index])
			plothist(dwell[rangehist & (np.arange(len(dwell))%2==0)], dwell[rangehist & (np.arange(len(dwell))%2!=0)], plotcsvpeak[index-1], plotcsvpeak[index])
			meanopentime = np.append(meanopentime, np.mean(dwell[rangehist & (np.arange(len(dwell))%2==0)]))
			meanclosetime = np.append(meanclosetime, np.mean(dwell[rangehist & (np.arange(len(dwell))%2!=0)]))
	csvtext = np.transpose(np.vstack((csvpeak, csvstart, csvend, popen, meanopentime, meanclosetime)))
	np.savetxt(pathfilename+'.csv', csvtext, delimiter=',')

	plotpeak = np.hstack((0, peak, peak + 1, np.cumsum(dwell)[-1]/JUMP-1))
	plotpeak = np.sort(plotpeak)*JUMP
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
