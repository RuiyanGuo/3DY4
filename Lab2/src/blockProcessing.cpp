/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <cmath>

#define PI 3.14159265358979323846
// function for computing the impulse response (reuse from previous experiment)
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// allocate memory for the impulse response
	float Norm = Fc/Fs*2;
	float hTemp;
	h.clear(); h.resize(num_taps, 0.0);

	for (int i = 0; i < num_taps; i++)
    {
		  if (i == (num_taps-1)/2)
			 {
				 hTemp = Norm;
			 }
			else
			 {
				 hTemp = Norm * (sin(PI*Norm*(i-(num_taps-1)/2)))/(PI*Norm*(i-(num_taps-1)/2));
			 }
			 hTemp = hTemp * pow(sin(i*PI/num_taps),2);
			 h[i] = hTemp;
	  }
}

// function for computing the impulse response (reuse from previous experiment)
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)  //my singlepass filter
{
	// allocate memory for the output (filtered) data
	float coeff_length = h.size();
	float data_length = x.size();
	float sum = 0;
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	for (int i = 0; i < data_length; i++)
	  {
		  sum = 0;
			for (int j = 0; j < coeff_length; j++)
			  {
				  if(i < j)
					 {
						 break;
					 }
					 sum = sum + h[j]*x[i-j];
			  }
			y[i] = sum;
	  }
	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
}

std::vector<float> slicing(std::vector<float>& arr, int a, int b)
{

    // Starting and Ending iterators
    auto start = arr.begin() + a;
    auto end = arr.begin() + b;

    // To store the sliced vector
    std::vector<float> result(b - a);

    // Copy vector using copy function()
    std::copy(start, end, result.begin());

    // Return the final sliced vector
    return result;
}

std::vector<float> mylfilter_w_block(std::vector<float> &coeff, std::vector<float> &data, int size, std::vector<float> &buffer)   // My lfilter for block processing
{
	int coeff_length = coeff.size();
	int data_length = data.size();
	int buffer_length = buffer.size();
	int num_special_treat = (coeff_length - 1);
	std::vector<float> buffer_out;
	std::vector<float> filtered_data;
	float sum = 0;
	int buffer_counter = 0;

	for (int i = 0; i < size; i++)
	{
		sum=0;
		buffer_counter = i+1;
		for (int j = (coeff_length-1); j >= 0; j--)
		{
			if (i > num_special_treat)
			 {
				if (i < j)
					{break;}
				sum = sum + coeff[j]*data[i-j];
				if (i == (size - 1))
					{buffer_out.push_back(data[i-j]);}
			 }
			else if ((i <= num_special_treat) && (buffer_counter < coeff_length))
			 {
				sum = sum + buffer[buffer_counter]*coeff[j];
				buffer_counter+=1;
			 }
			else if ((i <= num_special_treat) && (buffer_counter >= coeff_length))
			 {
				sum = sum + coeff[j]*data[i-j];
			 }
		}//inner for
		filtered_data.push_back(sum);
		//print(str(i) + "    " +str(data_length))
	}//outer for
  buffer = buffer_out;
	return filtered_data;
}

std::vector<float> myBlockfilter_implementation(std::vector<float> &audio_data, int block_size, float audio_Fc, float audio_Fs, unsigned short int N_taps)  // My implementation for block filter
{
	// derive filter coefficients
	std::vector<float> h;
	std::vector<float> sliced_data;
	impulseResponseLPF(audio_Fs, audio_Fc, N_taps, h);
	//coeff = myCoeff(audio_Fs, audio_Fc, N_taps)
	int coeff_length = h.size();
	std::vector<float> buffer;
	buffer.clear();
	buffer.resize(coeff_length, 1);

  std::vector<float> temp;

	int data_length = audio_data.size();
	std::vector<float> filtered_data;
	filtered_data.clear();
	filtered_data.resize(data_length, 0);
	int position = 0;

	while (true)
	{
		sliced_data = slicing(audio_data, position, position+block_size);
		temp = mylfilter_w_block(h, sliced_data, block_size, buffer);
		for (int i = position; i < position+block_size; i++)
		  {
			  filtered_data[i] = temp[i - position];
		  }
		position += block_size;
		if ((position + block_size) > data_length)
		 {
			break;
		 }
	}

	return filtered_data;
}
// function to read audio data from a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the Python script that can prepare this type of files
// directly from .wav files
void read_audio_data(const std::string in_fname, std::vector<float> &audio_data)
{
	// file descriptor for the input to be read
	std::ifstream fdin(in_fname, std::ios::binary);
	if(!fdin) {
		std::cout << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cout << "Reading raw audio from \"" << in_fname << "\"\n";
	}
	// search for end of file to count the number of samples to be read
	fdin.seekg(0, std::ios::end);
	// we assume the Python script has written data in 32-bit floats
	const unsigned int num_samples = fdin.tellg() / sizeof(float);

	// allocate memory space to store all the samples
	audio_data.clear(); audio_data.resize(num_samples);
	// back to the beginning of the file to read all samples at once
	fdin.seekg(0, std::ios::beg);
	// do a single read for audio data from the input file stream
	fdin.read(reinterpret_cast<char*>(&audio_data[0]), \
						num_samples*sizeof(float));
	// close the input file
	fdin.close();
}

// function to split an audio data where the left channel is in even samples
// and the right channel is in odd samples
void split_audio_into_channels(const std::vector<float> &audio_data, std::vector<float> &audio_left, std::vector<float> &audio_right)
{
	for (unsigned int i=0; i<audio_data.size(); i++) {
		if (i%2 == 0)
			audio_left.push_back(audio_data[i]);
		else
			audio_right.push_back(audio_data[i]);
	}
}

// function to write audio data to a binary file that contains raw samples
// represented as 32-bit floats; we also assume two audio channels
// note: check the python script that can read this type of files
// and then reformat them to .wav files to be run on third-party players
void write_audio_data(const std::string out_fname, const std::vector<float> &audio_left, const std::vector<float> &audio_right)
{
	// file descriptor for the output to be written
	if (audio_left.size() != audio_right.size()) {
		std::cout << "Something got messed up with audio channels\n";
		std::cout << "They must have the same size ... exiting\n";
		exit(1);
	} else {
		std::cout << "Writing raw audio to \"" << out_fname << "\"\n";
	}
	std::ofstream fdout(out_fname, std::ios::binary);
	for (unsigned int i=0; i<audio_left.size(); i++) {
		// we assume we have handled a stereo audio file
		// hence, we must interleave the two channels
		// (change as needed if testing with mono files)
		fdout.write(reinterpret_cast<const char*>(&audio_left[i]),\
								sizeof(audio_left[i]));
		fdout.write(reinterpret_cast<const char*>(&audio_right[i]),\
								sizeof(audio_right[i]));
	}
	fdout.close();
}

int main()
{
	// assume the wavio.py script was run beforehand to produce a binary file
	const std::string in_fname = "../data/float32samples.bin";
	// declare vector where the audio data will be stored
	std::vector<float> audio_data;
	// note: we allocate memory for audio_data from within this read function
	read_audio_data(in_fname, audio_data);

	// set up the filtering flow
	float Fs = 44100.0;	// sample rate for our "assumed" audio (change as needed for 48 ksamples/sec audio files)
	float Fc = 10000.0;	// cutoff frequency (explore ... but up-to Nyquist only!)
	int block_size = 1000;
	// number of FIR filter taps (feel free to explore ...)
	unsigned short int num_taps = 51;

	// impulse response (reuse code from the previous experiment)
	std::vector<float> h;
	impulseResponseLPF(Fs, Fc, num_taps, h);
	// note: memory for the impulse response vector and output data vectors
	// should be allocated from within the corresponding functions
	// (as for the previous experiment, from where you should reuse your code)

	// there is one more point before filtering is done:
	// recall we assume there are two channels in the audio data
	// the channels must be handled separately by your DSP functions, hence
	// split the audio_data into two channels (audio_left and audio_right)

	// declare vectors where the audio left/right channels will be stored
	std::vector<float> audio_left, audio_right;
	// note: we allocate the memory for the left/right channels
	// from within the split function that is called in the code below
	split_audio_into_channels(audio_data, audio_left, audio_right);

	// convolution code for filtering (reuse from the previous experiment)

 /*	Code for single pass
	std::vector<float> single_pass_left, single_pass_right;
	convolveFIR(single_pass_left, audio_left, h);
	convolveFIR(single_pass_right, audio_right, h);
 */
  //Code for block processing
	std::vector<float> block_left, block_right;
	block_left = myBlockfilter_implementation(audio_left, block_size, Fc, Fs, num_taps);
	block_right = myBlockfilter_implementation(audio_right, block_size, Fc, Fs, num_taps);

	// note: by default the above convolution produces zero on the output stream
	// YOU will need to update the convolveFIR and impulseResponseLPF functions
	// create a binary file to be read by wavio.py script to produce a .wav file
	// note: small adjustments will need to be made to wavio.py, i.e., you should
	// match the filenames, no need for self-checks as default Python code, ...
	const std::string out_fname = "../data/float32filtered.bin";
	//write_audio_data(out_fname, single_pass_left,	single_pass_right);    //For singlepass
	write_audio_data(out_fname, block_left,	block_right);                  //For Block processing

	return 0;
}
