% Class implementing the Viterbi Algorithm,
% used to do prediction of beat times.

classdef bp_viterbi < beat_predictor

properties
	num_features;

	% how many samples do we have for each observation frame, from each feature?
	feature_frame_length;

	% conversion factor from observation frame to time
	feature_sample_rate;

	% how many features do we have?
	num_feature_frames;

	% just counts how many frames have been added
	frame_index;

end

methods
	function initialise(this, num_features, feature_sample_rate, ...
			feature_frame_length, num_feature_frames)
		this.num_features = num_features;
		this.feature_frame_length = feature_frame_length;
		this.feature_sample_rate = feature_sample_rate;
		this.num_feature_frames = num_feature_frames;
		this.frame_index = 1;

	end

	% gives the Viterbi algorithm a frame of observations (feature data), so that it
	% can compute the most likely state for the current frame. The set of tempo
	% estimates is used to narrow down the search of possible tempos to only a small
	% subset identified by the tempo/phase estimator as likely.

	% PARAMETERS:
	% tempo_estimates = cell(1, this.num_features)
	%	  is a cell matrix containing a set of tempo estimates (no confidences)
	%	  for the current frame. There is one set for each feature, which means the cell
	%	  matrix is 1 row (since it's only one frame's worth of data) by n columns, where
	%	  n is the number of features. The tempo estimates should be in SECONDS.

	% feature_data = cell(1, this.num_features)
	%	  is a 1xn cell matrix containing the raw features, or observations.
	%	  There is a column in the cell matrix for each feature, and each cell is a
	%	  (single dimensional) vector of length equal to the feature frame length.

	function add_observations(this, feature_data, tempo_estimates)
		% decide which tempos to use group similar tempo estimates from different
		% features, adding more weight to each of them from tpe_autocorrelation tempo
		% peak picking: if estimates are closer than min_lag/2, treat them as the
		% same tempo.

		% values from tempo_phase_estimator
		max_bpm = 320;
		min_lag = 60/max_bpm;
		% value from tpe_autocorrelation
		MAX_TEMPO_PEAKS = 8;

		% stores the set of tempos to search over in this iteration of the
		% Viterbi algorithm
		candidate_tempos = zeros(this.num_features*MAX_TEMPO_PEAKS, 1)
		candidate_tempo_idx = 1;

		for n = 1:this.num_features
			tempo_estimates_n = tempo_estimates{n};
			for estimate_idx = 1:length(tempo_estimates_n)
				candidate_tempos(candidate_tempo_idx) = tempo_estimates_n(estimate_idx, 1);
				candidate_tempo_idx = candidate_tempo_idx + 1;
			end
		end
		% trim list to its actual length
		candidate_tempos = candidate_tempos(1:end)

		% now group together similar tempo estimates; use mean shift clustering.
		clustered_candidate_tempos = mean_shift_cluster(candidate_tempos, min_lag/2);


	end

end % methods


% MODEL SUMMARY:
% Each time the tempo and beat location needs to be estimated, we use the Viterbi
% algorithm to choose the most probable tempo period and beat location, according to
% the following model:

% for the kth set of feature frames (one for each of the n features), we have a state
% vector: X[k] = (T[k], B[k], M[k], W_i[k]). W_i are a set of n weights, one for each
% feature.

% T[k] is the dominant tempo period, measured in seconds, for that frame

% B[k] is the location of the beat (in seconds, measured from the beginning of the
% audio) closest to the time at which the current feature frame ends (but not after)

% M[k] is an indicator variable which represents the rhythmic meter that the current
% slice of audio is in.  I'm not sure whether to make this for simple or compound
% time, or duple or triple time.  This may be used to determine the shape of the beat
% alignment function, or to look at harmonics in the autocorrelation

% W_i[k] are n weight variables, one for each feature, which are used to determine
% how important each feature's onsets are in predicting beat locations. We model
% these as hidden variables, but there may be a better way of doing it. This is
% because the relevance of a feature may change in different sections of the music,
% and so the weights are allowed to vary over time, according to some observed
% confidence measure which we have yet to determine exactly.

% Transition probabilities essentially maintain continuity between estimates, and are
% all reasonably intuitive. Progression of tempo, meter and weights are assumed
% independent, and next beat location is dependent only on the previous tempo as well
% as its own previous value.  Each probability is more or less a Gaussian
% distribution centred around a projection of the previous state one time step
% forward (in the case of the beat time) or just centred around the previous value
% itself.

% The rest of this description focuses on the observation probabilities, which encode
% musical knowledge.  For the latter, we define F1, F2, F3, ... to be each of the n
% features, and write P(Observations | State) := P(F1, F2, F3, ... | T[n] = t, B[n] =
% b, M[n] = m, W_i[n] = w_i)           (1)

% We assume that, given the state, the observed features are independent (since they
% are meant to measure different aspects of the audio).  So then, we can split up the
% joint probability into a product for each feature F_i

% = (product over i of) P(F_i | T[n] = t, B[n] = b, M[n] = m, W_i[n] = w_i))     (2)

% Now the probability of observing the feature frame, given the tempo, beat location,
% and meter, is assumed to be proportional to the strength of the autocorrelation
% function at lag t, multiplied by a sequence of impulses spaced t seconds a part and
% shifted so that the latest one is at b seconds, multiplied by some metric of
% musical measure, for example the energy of the harmonics of the autocorrelation at
% harmonics 2 and 4 vs 3 and 6 (this comes from Davies and Plumbley,
% "Context-Dependent Beat tracking of Musical Audio") in equation form, this gives
% (up to a normalising constant for each term)

% = (product over i of) ACF_i(t)*BAF_i(t, b)*metric_measure_i(t, m)              (3)

% where ACF_i is the autocorrelation of feature i's frame, BAF is the 'beat alignment
% function' i.e. the inner product of the t-spaced impulses with feature i and
% metric_measure_i is some metric on the feature used to establish the meter

% the last problem is to integrate the weights. If we simply multiply each term by a
% weight, it doesn't really make sense. What's the probability of a observing feature
% given its weight?  I took a liberal view, and used the weight as a measure of
% 'how much discriminatory power' each feature has. If the weight is low, then we
% shouldn't be 'listening' to that feature so much if it's high, that feature has
% 'an important measure' of the observation probability.  One way to disregard one
% particular feature's term was to make the probability artificially higher, so a
% mismatch between the model and that feature's data isn't penalised so much.
% Therefore, what could be done is to raise the previous probabilities in (3) to the
% power of the weight.  When the weight is 0, the probability term becomes 1, thus
% effectively removing it from the product.  When the weight is 1, then the
% probability is unchanged. % This gives,

% P(Observations | State)
% := P(F1, F2, F3, ... | T[n] = t, B[n] = b, M[n] = m, W_i[n] = w_i)
%  = (product over i of) (ACF_i(t)*BAF_i(t, b)*metric_measure_i(t, m))^w_i       (4)

% A restriction could be imposed such that the weights must multiply to give 1, or
% 1/e.  The latter has the nice property that, taking the logarithm of (4) gives a
% weighed sum of log-probabilities such that the magnitudes of the weights sum to 1

% One problem with the model as it stands is its large search space. Other than
% discretising some of the variables (e.g weights are multiples of 0.1), one way to
% make the search tractable is to pick the peaks of the autocorrelation function for
% each feature - as is currently done - and evaluate state candidates only at these
% likely tempo values. If S tempo values are allowed as candidates at each iteration,
% and W is the number of weight values allowed, then the number of calculations is of
% the order of S*S*W*n, where n is the number of features.  Does restricting S to
% maybe 16 mean that the use of the Viterbi Algorithm is unwarranted?


% prior probability distribution on tempos
% lognormal distribution, with parameters taken
% from Klapuri's paper (Equation 22).
% He takes it from Parncutt (1994) "A perceptual model
% of pulse salience"
methods (Static)

	function p = tempo_prior_prob(t1)
		% parameters for the distribution.
		% These are genre dependent! But Klapuri uses these anyway.
		scale = 0.55;
		shape = 0.28;
		% p = lognpdf(t1, log(scale), shape);
		p = exp(-1*(log(t1/scale))^2/(2*shape^2)) / (t1*shape*sqrt(2*pi));
	end


	% probability of a tempo change between any two frames
	% in a song. t1 and t2 are the previous and next tempo periods respectively,
	% measured in seconds
	% This is taken from Klapuri's paper (equation 21).
	function p = tempo_transition_prob(t1, t2)
		% first term is a normal distribution as a function
		% of the ratio of t1 and t2 - doubling tempo is as likely
		% as halving it.

		ratio_sigma = 0.2; % variance for normal distribution
		ratio_prob = exp(-0.5*(ln(t2/t1)/ratio_sigma)^2) / sqrt(2*pi*ratio_sigma^2);
		p = ratio_prob * tempo_prior_prob(t2);
	end

	% transition probability of beat position (B[k-1] = b1, T[k-1] = t1) -> B[k] = b2
	% meant to reflect the fact that the expected beat location is
	% some multiple of t1 seconds after the previous one.
	% we make the transition probability a normal distribution
	% centred around the expected beat time.
	% Unlike the previous one, we don't use a prior, since the prior was
	% uniform over times
	function p = beat_position_transition_prob(b1, t1, b2)
		% choose k to minimise b2 - (b1 + k*t1)
		k = round((b2 - b1)/t1);

		% heuristic: 95% of samples lie within 2 standard deviations
		% of the mean. maybe 95% of beats lie within a quaver (tempo period/8)
		% of the expected time. then 2*sigma = t1/8
		beat_sigma = t1/16;
		% p = normcdf(b2, b1 + k*t1, beat_sigma);
		p = exp(-0.5*((b2 - (b1 + k*t1))/beat_sigma)^2)/sqrt(2*pi*beat_sigma^2);
	end

	% beats can start anywhere in the audio, make this uniform
	% but over what period of time? I suppose any uniform distribution is
	% proportional to 1, and we're going to be use proportional probabilities anyway.
	function p = beat_position_prior_prob(~)
		warning('Using unnormalised beat_position prior');
		p = 1;
	end

	% Transition probability for the meter variable M[k-1] = m1 -> M[k] = m2
	% for the moment this is just a small, constant probability of changing,
	% since the transition probability does not include current observations
	% if we know the meter to be one thing, we don't usually is would just
	% change randomly, but making the probability zero doesn't seem like a
	% good idea either. How do we decide what it is in the first place?
	function p = meter_variable_transition_prob(m1, m2)
		% probability of changing meter in any given frame?
		% quite small. but it's obviously dependent on the piece
		change_prob = 0.01;
		% check for valid values
		if (strcmp(m2, 'simple') || strcmp(m2, 'compound')) ...
				&& (strcmp(m1, 'simple') || strcmp(m1, 'compound'))
			if strcmp(m2, m1)
				p = 1 - change_prob;
			else
				p = change_prob;
			end
		else
			p = 0;
		end
	end

	% more pieces are probably simple time than complex, but should
	% that be modelled here? For now, we use a uniform distribution.
	function p = meter_variable_prior_prob(m1)
		if strcmp(m1, 'simple')
			p = 0.5;
		elseif strcmp(m1, 'compound')
			p = 0.5;
		else
			p = 0;
		end
	end

	% all weights equal initially. Should the product be 1/e?
	% In which case the initial weights are e^(1/n), where n is the number
	% of features
	function v = feature_weight_initial_value
		% have to pass in number of features from somewhere
		warning('Number of features undefined');
		v = exp(-1/n);
	end



	% need an Observation and a State class

	% Compute the forward message of the hidden Markov model, given the most recent
	% observations. This is a distribution over all possible states.
	% in other words, do one step of filtering.
	% the initial forward message is the initial prior.

	% THIS WILL NEED OPTIMISING!!!

	% PARAMETERS
	% current_message = matrix(num_tempos, num_beat_alignments)
	% 	is the probability of each state given observations up to this point
	% current_feature_frame = matrix(feature_frame_length, num_features);
	%   is the matrix of feature data from each feature, for the current 
	%   feature frame
	function new_message = update_forward_message(current_message, current_feature_frame)
		% new_forward_message(X_t)
		% = Prob(X_t | e_{1...t})
		% = K*(sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)
		%		*P(current_observations | X_t)
		% where K is a normalising constant, so that the sum of
		% probabilities for each X_t is 1.

		% this will store the feature frame when it's updated each time.
		% rows index time, from earliest to most recent feature samples,
		% while columns index different features.

		% stores autocorrelation data, indexed in a similar way to above.
		% currently, autocorrelation returns an array of half the length of the
		% feature frame due to the way it works.

		autocorrelation_data = zeros(feature_frame_length/2, num_features);
		beat_alignment_data  = zeros(feature_frame_length/2, ...
			feature_frame_length/2, num_features);

		% state: (tempo, beat alignment) ... for now
		new_message = zeros(num_tempos, num_beat_alignments);
		observation_probability = zeros(num_tempos, num_beat_alignments);
		transition_probability = zeros(num_tempos, num_beat_alignments);

		% for all features
			% calculate acf
			% (for all meter variables - later)
				% for all tempos
					% calculate beat alignment
					% for all beat alignments

		% EXAMPLE VALUES: FIX
		all_tempos = 1:1024;
		all_beat_alignments = 1:1024;

		for n = 1:num_features % -> currently we're just adding up all features equally!
			feature_frame_n = current_feature_frame(:, n);
			autocorrelation_data(:, n) = autocorrelation(feature_frame_n);
			for t = all_tempos
				beat_alignment_data(:, t, n) = beat_alignment_function(feature_frame_n, t);
				for b = all_beat_alignments
					%state = [t, b];
					% unnormalised
					observation_probability(t, b) = observation_probability(t, b) + ...
						autocorrelation_data(t, n)*beat_alignment_function(b, t, n);
				end
			end
		end

		% need to calculate transition probability from all states to all
		% other states?????? 
		for t2 = all_tempos
			for b2 = all_beat_alignments
				% new_state = (t2, b2)
				% this sums the probabilities of going from each of the old states
				% to the (fixed) new state (t2, b2), ignoring new information
				trans_prob_t2b2 = 0;
				for t1 = all_tempos
					for b1 = all_beat_alignments % oh god
						% old_state = (t1, b1)
						P_t1_t2 = tempo_transition_prob(t1, t2);
						P_b1_t1_t2 = beat_position_transition_prob(b1, t1, b2);
						trans_prob_t2b2 = trans_prob_t2b2 + ...
							P_t1_t2*P_b1_t1_t2*current_message(t1, b1);
					end
				end
				transition_probability(t2, b2) = trans_prob_t2b2;
			end
		end

		% now multiply them and make it sum to 1

		sum_updated_probs = 0;
		for t = all_tempos
			for b = all_beat_alignments
				% probability of each state given new observations
				updated_prob_tb = observation_probability(t, b)*transition_probability(t, b);
				sum_updated_probs = sum_updated_probs + updated_prob_tb;
				new_message(t, b) = updated_prob_tb;
			end
		end

		new_message = new_message/sum_updated_probs;
	end

	% uses the prior distribution to generate an initial forward message
	% P(X_0 | null), i.e. the probability of any state given no information
	function new_message = generate_initial_forward_message
		new_message = zeros(num_tempos, num_beat_alignments);

		% EXAMPLE VALUES: FIX
		all_tempos = 1:1024;
		all_beat_alignments = 1:1024;

		% for all states
		for t = all_tempos
			for b = all_beat_alignments
				% state = (t, b);
				new_message(t, b) = tempo_prior_prob(t)*beat_position_prior_prob(b);
			end
		end

		% make the whole thing sum to 1 - as some priors are only proportional
		new_message_sum = sum(sum(new_message));
		new_message = new_message/new_message_sum;
	end

end


end

