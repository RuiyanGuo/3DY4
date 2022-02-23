#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import sys
import math

def filter_block_processing(audio_data, \
       block_size, \
       audio_Fc, \
       audio_Fs, \
       N_taps):

 # derive filter coefficients
 firwin_coeff = signal.firwin(N_taps,
        audio_Fc/(audio_Fs/2),
        window=('hann'))

 # we assume the data is stereo as in the audio test file
 filtered_data = np.empty(shape = audio_data.shape)

 # start at the first block (with relative position zero)
 position = 0

 # intiial filter state - state is the size of the impulse response minus 1
 # we need to channels for the state (one for each audio channel)
 filter_state = np.zeros(shape = (len(firwin_coeff)-1, 2))

 while True:

  # filter both left and right channels
  filtered_data[position:position+block_size, 0], filter_state[:,0] = \
   signal.lfilter(firwin_coeff, 1.0, \
   audio_data[position:position+block_size, 0], zi = filter_state[:,0])
  # the filter state has been saved only for the first channel above
  # you will need to adjust the code for the second channel below
  filtered_data[position:position+block_size, 1], filter_state[:,1] = \
   signal.lfilter(firwin_coeff, 1.0, \
   audio_data[position:position+block_size, 1], zi = filter_state[:,1])

  position += block_size

  # the last incomplete block is ignored
  if position > len(audio_data):
   break

 return filtered_data

def filter_single_pass(audio_data, audio_Fc, audio_Fs, N_taps):

	# derive filter coefficients
	firwin_coeff = signal.firwin(N_taps,
								audio_Fc/(audio_Fs/2),
								window=('hann'))

	# we assume the data is stereo as in the audio test file
	filtered_data = np.empty(shape = audio_data.shape)

	# filter left channel
	filtered_data[:,0] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,0])
	# filter stereo channel
	filtered_data[:,1] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,1])
'''
	# filter left channel
	filtered_data[:,0] = mylfilter(myCoeff(audio_Fs, audio_Fc, N_taps), audio_data[:,0])
	# filter stereo channel
	filtered_data[:,1] = mylfilter(myCoeff(audio_Fs, audio_Fc, N_taps), audio_data[:,1])

	# filter left channel
	filtered_data[:,0] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,0])
	# filter stereo channel
	filtered_data[:,1] = signal.lfilter(firwin_coeff, 1.0, audio_data[:,1])
'''


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

def mylfilter(coeff, data):   # My lfilter for singlepass
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
		print(str(i) + "    " +str(data_length))
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
		print(str(i) + "    " +str(data_length))
	return buffer_out, filtered_data

def myBlockfilter_implementation(audio_data, block_size, audio_Fc, audio_Fs, N_taps):  # My implementation for block filter
	# derive filter coefficients
	coeff = myCoeff(audio_Fs, audio_Fc, N_taps)
	coeff_length = len(coeff)
	buffer1 = np.ones(coeff_length)
	buffer2 = np.ones(coeff_length)
	data_length = len(audio_data[:,0])
	filtered_data = np.empty(shape = audio_data.shape)
	position1 = 0
	position2 = 0

	while True:
		buffer1, filtered_data[position1:position1+block_size, 0] = mylfilter_w_block(coeff, audio_data[position1:position1+block_size, 0], block_size, buffer1)
		buffer2, filtered_data[position2:position2+block_size, 1] = mylfilter_w_block(coeff, audio_data[position2:position2+block_size, 1], block_size, buffer2)

		position1 += block_size
		position2 += block_size
		if (position1 + block_size) > (data_length):
			break

	return filtered_data

def cli_error_msg():

	# error message to provide the correct command line interface (CLI) arguments
	print('Valid arguments:')
	print('\trc:  reference code')
	print('\til1: in-lab 1')
	print('\tth:  take-home')
	sys.exit()

# audio test file from: https://www.videvo.net/royalty-free-music/
if __name__ == "__main__":

	if len(sys.argv[0:]) != 2:
		cli_error_msg()

	# use use wavfile from scipy.io for handling .wav files
	print('Opening audio file (.wav format)')
	audio_Fs, audio_data = wavfile.read("../data/audio_test.wav")
	print(' Audio sample rate = {0:d} \
		\n Number of channels = {1:d} \
		\n Numbef of samples = {2:d}' \
		.format(int(audio_Fs), audio_data.ndim, len(audio_data)))

	if (sys.argv[1] == 'rc'): # runs the reference code (rc)

		print('Reference code for processing streams divided in blocks')

		# you can control the cutoff frequency and number of taps
		single_pass_data = filter_single_pass(audio_data, \
							audio_Fc = 10e3, \
							audio_Fs = audio_Fs, \
							N_taps = 51)

		# write filtered data back to a .wav file
		wavfile.write("../data/single_pass_filtered.wav", \
					audio_Fs, \
					single_pass_data.astype(np.int16))

	elif (sys.argv[1] == 'il1'):

		print('In-lab experiment 1 for processing streams divided in blocks')

		# you can control also the block size
		block_processing_data = filter_block_processing(audio_data, \
							block_size = 1000, \
							audio_Fc = 10e3, \
							audio_Fs = audio_Fs, \
							N_taps = 51)

		wavfile.write("../data/block_processing_filtered.wav", \
					audio_Fs, \
					block_processing_data.astype(np.int16))

	elif (sys.argv[1] == 'th'):

		print('Take-home exercise for processing streams divided in blocks')

		# for specific details check the lab document
		block_processing_data = myBlockfilter_implementation(audio_data, block_size = 1000, audio_Fc = 10e3, audio_Fs = audio_Fs, N_taps = 51)

		wavfile.write("../data/block_processing_filtered.wav", \
					audio_Fs, \
					block_processing_data.astype(np.int16))

		# it is suggested that you add plotting while troubleshooting
		# if you plot in the time domain, select a subset of samples,
		# from a particular channel (or both channels) e.g.,
		# audio_data[start:start+number_of_samples, 0]

	else:

		cli_error_msg()

	# plt.show()
