#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
# for take-home add your functions

rf_Fs = 2.4e6
if_Fs = 240e3

rf_Fc = 100e3
if_Fc = 16e3

rf_taps = 151
if_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_decim = 5
# add other settings for audio, like filter taps, ...
def myCoeff(Fs, Fc, Ntaps):
	Norm = Fc/Fs*2
	arr = []

	for i in range(Ntaps):
		if i == (Ntaps-1)/2:
			hTemp = Norm
		else:
			hTemp = Norm * (math.sin(math.pi*Norm*(i-(Ntaps-1)/2)))/(math.pi*Norm*(i-(Ntaps-1)/2))
		hTemp = hTemp * math.pow(math.sin(i*math.pi/Ntaps),2)
		arr.append(hTemp)
	return arr

def mylfilter(coeff, data):    #My lfilter for single pass
	coeff_length = len(coeff)
	data_length = len(data)
	sum = 0
	arr = np.zeros(data_length)

	for i in range(data_length):
		sum=0
		for j in range(coeff_length):
			if i < (j+1):
				break
			sum = sum + coeff[j]*data[i-j-1]
		arr[i] = sum
		#print(str(i) + "    " +str(data_length))
	return arr

def mylfilter_w_block(coeff, data, size, buffer):    # My lfilter for block processing
	coeff_length = len(coeff)
	data_length = len(data)
	buffer_length = len(buffer)
	num_special_treat = coeff_length - 1
	buffer_out = []
	filtered_data = []
	sum = 0
	buffer_counter = 0

	for i in range(0, size):
		sum=0
		buffer_counter = i+1
		for j in range(coeff_length-1, -1, -1):
			if i > num_special_treat:
				if i < j:
					break
				sum = sum + coeff[j]*data[i-j]
				if i == (size - 1):
					buffer_out.append(data[i-j])
			elif ((i <= num_special_treat) and (buffer_counter < coeff_length)):
				sum = sum + buffer[buffer_counter]*coeff[j]
				buffer_counter+=1
			elif ((i <= num_special_treat) and (buffer_counter >= coeff_length)):
				sum = sum + coeff[j]*data[i-j]
		filtered_data.append(sum)
		#print(str(i) + "    " +str(data_length))
	return buffer_out, filtered_data

def myDemod(I, Q, prev_i = 0.0, prev_q = 0.0):    #My demod function

	fm_demod = np.empty(len(I))
	numerator = 0
	denominator = 0

	for i in range(len(I)):
		numerator = (I[i] * (Q[i]-prev_q)) - (Q[i] * (I[i]-prev_i))
		denominator = math.pow(I[i], 2) + math.pow(Q[i], 2)
		fm_demod[i] = numerator/denominator
		prev_i = I[i]
		prev_q = Q[i]

	return fm_demod, prev_i, prev_q
# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
il_vs_th = 1

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/iq_samples.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)
		audio_coeff = signal.firwin(if_taps, if_Fc/(if_Fs/2), window=('hann'))
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = myCoeff(if_Fs, if_Fc, if_taps)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_phase = 0
	prev_i = 0
	prev_q = 0
	# add state as needed for the mono channel filter

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)
	#audio_data = []
	coeff_length = len(audio_coeff)
	buffer = np.ones(coeff_length)
	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		if il_vs_th == 0:
			# already given to you for the in-lab
			# take particular notice of the "special" state-saving
			fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

		else:
			# you will need to implement your own FM demodulation based on:
			# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
			# see more comments on fmSupportLib.py - take particular notice that
			# you MUST have also "custom" state-saving for your own FM demodulator
			fm_demod, prev_i, prev_q = myDemod(i_ds, q_ds, prev_i, prev_q)

		# extract the mono audio data through filtering
		if il_vs_th == 0:
		# 	# to be updated by you during the in-lab session based on lfilter
		# 	# same principle as for i_filt or q_filt (but different arguments)
		 	audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod) #change as needed
		else:
			#audio_filt = mylfilter(audio_coeff, fm_demod)
			buffer, audio_filt = mylfilter_w_block(audio_coeff, fm_demod, len(fm_demod), buffer)
		# 	# to be updated by you for the takehome exercise
		# 	# with your own code for BLOCK convolution
		# 	audio_filt = ... change as needed
		# downsample audio data
		# to be updated by you during in-lab (same code for takehome)
		audio_block = audio_filt[::audio_decim] #change as needed

		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		#
		audio_data = np.concatenate((audio_data, audio_block))
		#

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			fmPlotPSD(ax1, audio_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Mono')
			# ... change as needed

			# plot PSD of selected block after downsampling mono audio
			fmPlotPSD(ax2, audio_data, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')
			# ... change as needed

			# save figure to file
			fig.savefig("../data/fmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	# plt.show()
