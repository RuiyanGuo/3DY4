## Software Defined Radio Project
Real-time SDR for mono/stereo FM and RDS

The project realizes a real-time implementation of a computing system that operates in a form factor-constrained environment.

Project Description:

Designed a Software Defined Radio system to transform real time signal into audio file

• Developed the project in Python for simplification and debug, translated to C++ for better
  performance
  
• Designed FIR filters to filter out signals with unwanted frequencies

• Designed a RDS demodulator with combined implementation of Rational Resampler and
  Root-raised cosine filter to demodulate the data for further processing
  
• Implemented Manchester and differential decoding technique for data processing

• Developed frame synchronization error detection algorithms to make sure that the decoded
  data is consistent with the input data
