
R2B2 aka 'Meat tracker'

To run: use the function R2B2 which looks really nice now

In its current form, the function computes tempo and phase estimates
and plots a sample of frames showing likely beat locations

***************
Class structure
***************

Abstract classes:

feature_extractor -
for something that takes audio and spits out feature data, which may consist of
one or several dimensions of data over time

tempo_phase_estimator - for something that takes in (single dimensional)
feature data and creates several tempo and phase estimates at regular intervals

beat_predictor - for something that takes in tempo and phase estimates
over time, aggregates them, and predicts a single most likely future beat time

Concrete classes:

odf_klapuri -
Implements Klapuri's onset detection function

tpe_autocorrelation -
Performs autocorrelation-based tempo detection, then for each tempo estimate
in a frame of features, finds phase by sliding set of equally spaced impulses
along the feature frame, and choosing the amount of shift that maximises their
dot product

tpe_autocorrelation2 -
Also performs autocorrelation-based tempo detection, but using a different method
of autocorrelation

Helper functions:

autocorrelation.m - Calculates autocorrelation
autocorrelation2.m - Calculates autocorrelation in a better way
mel_filter.m - for klapuri's feature, and also MFCC calculation (later)
beat_alignment_function.m - for phase calculation
