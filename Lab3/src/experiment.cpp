/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#define PI 3.14159265358979323846

std::vector<float> arange(float start, float step, float end)
{
	std::vector<float> result;
	for (float i = start; i < end; i+=step)
	  {
		  result.push_back(i);
	  }
	return result;
}

std::vector<std::complex<float>> slicing(std::vector<std::complex<float>> &arr, int a, int b)
{

    // Starting and Ending iterators
    auto start = arr.begin() + a;
    auto end = arr.begin() + b;

    // To store the sliced vector
    std::vector<std::complex<float>> result(b - a);

    // Copy vector using copy function()
    std::copy(start, end, result.begin());

    // Return the final sliced vector
    return result;
}

std::vector<float> estimatePSD(std::vector<float> &samples, std::vector<float> &freq, int nfft, float Fs)
{
	// rename the NFFT argument (notation consistent with matplotlib.psd)
	// to freq_bins (i.e., frequency bins for which we compute the spectrum)
	int freq_bins = nfft;
	// frequency increment (or resolution of the frequency bins)
	float df = Fs/freq_bins;
	// create the frequency vector to be used on the X axis
	// for plotting the PSD on the Y axis (only positive freq)
	freq = arange(0, Fs/2, df);
	// design the Hann window used to smoothen the discrete data in order
	// to reduce the spectral leakage after the Fourier transform
	std::vector<float> hann;
	hann.resize(freq_bins);

  for(int i = 0; i < freq_bins; i++)
	  {
		  hann[i] = pow(sin(i*PI/freq_bins),2);
	  }

	// create an empty list where the PSD for each segment is computed
	std::vector<float> psd_list;
	psd_list.clear();

	// samples should be a multiple of frequency bins, so
	// the number of segments used for estimation is an integer
	// note: for this to work you must provide an argument for the
  // number of frequency bins not greater than the number of samples!
	int no_segments = int(floor((samples.size())/float(freq_bins)));

  std::vector<float> windowed_samples;
	// iterate through all the segments
  for (int k = 0; k < no_segments; k++)
	{
		// apply the hann window (using pointwise multiplication)
		// before computing the Fourier transform on a segment
		windowed_samples.clear();

   for (int n = 0; n < int(hann.size()); n++)
    {
	    windowed_samples.push_back(samples[k*freq_bins+n]*hann[n]);
    }

		//compute the Fourier transform using the built-in FFT from numpy
		std::vector<std::complex<float>> Xf;

    DFT(windowed_samples, Xf);
		// note, you can check how MUCH slower is DFT vs FFT by replacing the
		// above function call with the one that is commented below
		//
		// Xf = DFT(windowed_samples)
		//
		// note: the slow impelementation of the Fourier transform is not as
		// critical when computing a static power spectra when troubleshooting
		//
		// note also: time permitting a custom FFT can be implemented

		// since input is real, we keep only the positive half of the spectrum
		// however, we will also add the signal energy of negative frequencies
		// to have a better a more accurate PSD estimate when plotting
		Xf = slicing(Xf, 0, int(freq_bins/2));  //keep only positive freq bins

		std::vector<float> psd_seg;
		float temp = 0;
		for (int j = 0; j < int(Xf.size()); j++)
		  {
				//temp = pow(abs(Xf[j]), 2) * (1/(Fs*freq_bins/2)) *2; // compute signal power, add the energy from the negative freq bins
				temp = pow(abs(Xf[j])/(freq_bins/2), 2); // compute signal power, add the energy from the negative freq bins
        psd_seg.push_back(temp);
		  }
		// translate to the decibel (dB) scale
		for (int q = 0; q < int(psd_seg.size()); q++)
		  {
			  psd_seg[q] = 10*log10(psd_seg[q]);
		  }
		// append to the list where PSD for each segment is stored
		// in sequential order (first segment, followed by the second one, ...)
		for(int w = 0; w < int(psd_seg.size()); w++)
		  {
			  psd_list.push_back(psd_seg[w]);
		  }
  }

	// compute the estimate to be returned by the function through averaging
	std::vector<float> psd_est;
	psd_est.resize(int(freq_bins/2), 0);
	// iterate through all the frequency bins (positive freq only)
	// from all segments and average them (one bin at a time ...)
  for (int e = 0; e < int(freq_bins/2); e++)
	  {
		  for (int r = 0; r < no_segments; r++)
			  {
				  psd_est[e] += psd_list[e + r*int(freq_bins/2)];
			  }
			psd_est[e] = psd_est[e] / no_segments;
	  }
	// the frequency vector and PSD estimate
	return psd_est;
}

int main()
{
	// binary files can be generated through the
	// Python models from the "../model/" sub-folder
	const std::string in_fname = "../data/fm_demod_10.bin";
	std::vector<float> bin_data;
	readBinData(in_fname, bin_data);

	// generate an index vector to be used by logVector on the X axis
	std::vector<float> vector_index;
	genIndexVector(vector_index, bin_data.size());
	// log time data in the "../data/" subfolder in a file with the following name
	// note: .dat suffix will be added to the log file in the logVector function
	logVector("demod_time", vector_index, bin_data);

	// take a slice of data with a limited number of samples for the Fourier transform
	// note: NFFT constant is actually just the number of points for the
	// Fourier transform - there is no FFT implementation ... yet
	// unless you wish to wait for a very long time, keep NFFT at 1024 or below
	std::vector<float> slice_data = \
		std::vector<float>(bin_data.begin(), bin_data.begin() + NFFT);
	// note: make sure that binary data vector is big enough to take the slice

	// declare a vector of complex values for DFT
  std::vector<std::complex<float>> Xf;
	// ... in-lab ...
	// compute the Fourier transform
	DFT(slice_data, Xf);
	// the function is already provided in fourier.cpp

	// compute the magnitude of each frequency bin
	// note: we are concerned only with the magnitude of the frequency bin
	// (there is NO logging of the phase response)
	std::vector<float> Xmag;
	// ... in-lab ...
	// compute the magnitude of each frequency bin
	// the function is already provided in fourier.cpp
	computeVectorMagnitude(Xf, Xmag);
	// log the frequency magnitude vector
	vector_index.clear();
	genIndexVector(vector_index, Xmag.size());
	logVector("demod_freq", vector_index, Xmag); // log only positive freq

	// for your take-home exercise - repeat the above after implementing
	// your OWN function for PSD based on the Python code that has been provided
	// note the estimate PSD function should use the entire block of "bin_data"
	//
	// ... complete as part of the take-home ...
	//
	float Fs = 240e3;
	std::vector<float> freq;
  std::vector<float> est_psd = estimatePSD(bin_data, freq, NFFT, Fs);
  std::cout << "Writing raw audio to \"" << bin_data.size() << "\"\n";
	vector_index.clear();
	genIndexVector(vector_index, est_psd.size());
	logVector("demod_psd", vector_index, est_psd);
	// if you wish to write some binary files, see below example
	//
	// const std::string out_fname = "../data/outdata.bin";
	// writeBinData(out_fname, bin_data);
	//
	// output files can be imported, for example, in Python
	// for additional analysis or alternative forms of visualization

	// naturally, you can comment the line below once you are comfortable to run GNU plot
	std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";

	return 0;
}
