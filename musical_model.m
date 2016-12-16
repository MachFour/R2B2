classdef musical_model
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
	% how to normalise this? If density values are above 1, is this a problem?
	function p = lognormal_density(x, mu, sigma)
		p = exp(-0.5*((log(x) - mu)./sigma).^2)./(x.*sigma*sqrt(2*pi));
	end

	function p = normal_density(x, mu, sigma)
		p = exp(-0.5*((   (x) - mu)./sigma).^2)./(1.*sigma*sqrt(2*pi));
	end

	function p = tempo_prior_prob(t1)
		% parameters for the distribution.
		% These are genre dependent! But Klapuri uses these anyway.
		%scale = 0.55;
		%shape = 0.28;
		%p = musical_model.lognormal_density(t1, log(scale), shape);

		% From: PREFERRED TEMPO RECONSIDERED - Dirk Moelants
		f0 = 125/60; % 125 BPM in Hz
		fext = 1./t1; % tempo period in Hz
		beta = 8;
		fext_2 = fext.^2;
		f0_2 = f0^2;
		fext_4 = fext.^4;
		f0_4 = f0^4;

		p = ((f0_2 - fext_2).^2 + beta*fext_2).^-0.5 - (fext_4 + f0_4).^-0.5;
	end


	% probability of a tempo change between any two frames
	% in a song. t1 and t2 are the previous and next tempo periods respectively,
	% measured in seconds
	% This is taken from Klapuri's paper (equation 21).
	function p = tempo_transition_prob(old_tempo, new_tempo)
		if old_tempo <= 0 || new_tempo <= 0
			error('tempo values must be positive');
		end
		% first term is a normal distribution as a function
		% of the ratio of t1 and t2 - doubling tempo is as likely
		% as halving it.

		% in fact, what's actually going on (Klapuri doesn't explain), is that we
		% say that log(t2/t1) is normally distributed with mean 0 and some variance.
		% This means that t2, given the value of t1, should be lognormally
		% distributed with mean equal to log(t1), and that variance.

		% ALTHOUGH, the mode of the lognormal distribution is not the same as
		% the mean, so there's a slight difference between the expected value
		% and the most likely value.
		% The mode of a lognormal distribution with parameters mu and sigma is
		% exp(mu - sigma^2) = exp(mu)exp(-sigma^2).
		% So in order to get the mode equal to t1, we set the mu parameter to be
		% log(t1)*exp(sigma^2)

		% Question: what's the variance? Should it be proportional to tempo?
		% Or the probability of the last state?
		%sigma = 0.1;
		%mu = log(t1)*exp(sigma^2);

		%p = musical_model.lognormal_density(t2, mu, sigma) * ...
		%	musical_model.tempo_prior_prob(t2);
		% tempo percentage change
		p = 1 - abs(new_tempo - old_tempo)/old_tempo;
	end

	% transition probability is meant to reflect the fact that the expected new beat
	% location is some multiple of the tempo after the previous one.
	% This must allow for an arbitrary frame hope size and so we must be able to
	% calculate the exact number of tempo multiples that are expected.

	% The end of the frame will advance by frame_hop_size samples, so the alignment
	% of the previous absolute beat time with respect to the new frame end time,
	% is b1 (which is a negative beat alignment) minus the frame hop size.
	% We then need to add multiples of t1 to get to a number in the range
	% [-t1 + 1, 0]. This is the new expected beat time. This cam be expressed as
	% b2' = mod(b1 - feature_hop_size, -1*t1)
	% (where the return value has the same sign as the divisor)

	function p = beat_location_transition_prob(old_state, new_state)
		hop_size = new_state.params.feature_hop_size;
		feature_sample_rate = new_state.params.feature_sample_rate;

		old_tempo = old_state.tempo_samples;
		old_alignment = old_state.beat_alignment;

		new_tempo = new_state.tempo_samples;


		% variance in beat location transition prob:
		% If the tempo hasn't changed, the difference in beat locations should be
		% almost exactly a multiple of the tempo. Assume that maybe 95% of intervals
		% between beats at the same tempo lie within a quaver (tempo period/8) of the
		% ideally expected time.

		% for samples from a normal distribution, it is  expected that 95% of them
		% lie within 2 standard deviations of its mean. Then 2*sigma = tempo_period/8

		% If the old and new tempos are different, it is not so clear what to
		% do. Which tempo period should we take?
		% For now, just take the average. This may be a terrible hack - but it has the
		% interpretation of a gradually changing tempo, maybe.

		tempo = new_tempo;
		tempo_period = tempo/feature_sample_rate;
		tempo_percent_change = abs(old_tempo - new_tempo)/old_tempo;

		if tempo <= 0
			error('tempo value must be positive');
		end

		% Rather than 1/16 of the tempo, use the percentage change in tempo,
		% times the tempo. (This becomes the absolute difference in tempo
		% times).

		sigma = tempo_period/16;
		if tempo_percent_change > 1/16
			sigma = sigma*(1+tempo_percent_change);
		end

		expected_alignment = mod(old_alignment - hop_size, -tempo);

		% make the probability a normal distrubution with mean equal to the
		% absolute beat location corresponding to the expected alignment
		expected_beat_loc = new_state.frame_end_time + ...
			expected_alignment/feature_sample_rate;

		new_beat_loc = new_state.beat_location;

		p = musical_model.normal_density(new_beat_loc, expected_beat_loc, sigma);

	end

	% beats can start anywhere in the audio, make this uniform
	% but over what period of time? I suppose any uniform distribution is
	% proportional to 1, and we're going to be use proportional probabilities anyway.
	function p = beat_position_prior_prob(~)
		%warning('Using unnormalised beat_position prior');
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

	function prob = transition_prob(old_state, new_state)
		t1 = old_state.tempo_period;
		t2 = new_state.tempo_period;

		% by independence
		prob = musical_model.tempo_transition_prob(t1, t2);
		%beat_location_prob = ...
		%	musical_model.beat_location_transition_prob(old_state, new_state);
		%prob = tempo_prob * beat_location_prob;
	end

	function prob = prior_prob(state)
		tempo_prior = musical_model.tempo_prior_prob(state.tempo_period);
		beat_pos_prior = musical_model.beat_position_prior_prob(state.beat_location);
		prob = tempo_prior*beat_pos_prior;
	end

end % methods (Static)


end % classdef
