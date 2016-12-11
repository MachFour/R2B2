
R2B2 - Real 2ime Beat 2racker.
Written by Max Fisher

*********
To run:
*********

The main function is contained in r2b2.m. It will read in an audio file,
calculate beat times, and output another audio file with the calculated
beats superimposed as beeps on top of the original audio. Note that the
algorithm only outputs predicted beat times after about 6 seconds of audio,
so listen for beeps after then.

It has only been tested with .wav files at 44.1kHz sample rate.

To run on an audio file called 'my_audio.wav', located in 'my_directory',
call the MATLAB function

>>> r2b2('my_audio', 'my_directory');

The script will then generate a file called 'annotated_my_audio.wav', and
put it into 'my_directory'. This file can then be listened to using any
normal music player application.

***************
Class structure
***************

Abstract classes:

feature_extractor -
for something that takes audio and spits out feature data, which may consist of
one or several dimensions of data over time

tempo_alignment_estimator - for something that takes in (single dimensional)
feature data and creates several tempo and beat alignment estimates at regular intervals

Concrete classes:

odf_klapuri -
Implements Klapuri's onset detection function

odf_mfcc -
Computes a feature made of differenced MFCCs

tae_autocorrelation -
Performs autocorrelation-based tempo detection, then for each tempo estimate
in a frame of features, finds beat alignment by sliding set of equally spaced impulses
along the feature frame, and choosing the amount of shift that maximises their
dot product

bp_viterbi - Uses the Viterbi algorithm to model probabilities of things, which can then
be used to takes in tempo and beat alignment estimates over time, aggregate them,
and predict most likely future beat times


Helper functions:

autocorrelation.m - Calculates autocorrelation in an improved way
mel_filter.m - for klapuri's feature, and also MFCC calculation
beat_alignment_function.m - for beat alignment calculation
