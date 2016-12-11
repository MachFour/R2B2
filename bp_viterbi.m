% bp_viterbi.m

% Class implementing the Viterbi Algorithm,
% used to do prediction of beat times.

% Author: Max Fisher

classdef bp_viterbi < handle

properties
	% counts the frame number
	frame_number;

	% each frame, we calculate the probability that the 'state' (i.e. tempo, beat
	% location, meter, etc ...) of the music in this frame is one of a set of
	% estimated 'likely' states. (Since calculating the probability for all possible
	% states is too much). These variables store a cell array of the states, and a
	% vector of probabilties to match.
	% The probabilities measure the probability that the 'true' (and unobservable)
	% state of the frame is equal to the corresponding entry in the current_states
	% vector, given the 'observations' made of the music so far.
	current_states;
	current_probabilities;
	current_tempos;

	winning_states;
	winning_probabilities;

	data_output_suffix = '-beat-times.txt';

	name;
	num_features;

	% processing_params object, used for sample to time conversions.
	params;

end

methods
	function viterbi = bp_viterbi(params, name, num_features)
		if nargin == 0
			warning('bp_viterbi initialised with default values');
			params = {};
			name = 'bp_viterbi';
			num_features = 0;
		end
		viterbi.params = params;
		viterbi.name = name;
		viterbi.num_features = num_features;

		viterbi.frame_number = 0;

	end

	% compute initial states and probabilities, given no tempos
	function initialise(this)
		% restrict to 50-180 BPM initially - low bpm makes it really
		% inefficient!!! (there are a huge number of states)
		min_initial_bpm = 45;
		max_initial_bpm = 180;

		sample_to_bpm_factor = this.params.feature_sample_rate*60;
		min_tempo_samples = round(sample_to_bpm_factor/max_initial_bpm);
		max_tempo_samples = round(sample_to_bpm_factor/min_initial_bpm);

		% choose a sparser set of tempos to start with
		initial_tempos = (min_tempo_samples:2:max_tempo_samples)';
		beat_coarseness = 6;
		initial_states = this.generate_all_states(initial_tempos, beat_coarseness);
		this.current_states = initial_states;
		this.current_probabilities = this.initial_forward_message(initial_states);

		% allocate some space for winning estimates and probabilities;
		this.winning_states = cell(50, 1);
		this.winning_probabilities = zeros(50, 1);

	end


	% converts tempo and beat alignment in samples to a model state which
	% measures. Note that beat alignment is negative
	function s = create_state(this, tempo_samples, beat_alignment)
		s = model_state(this.params, this.frame_number, ...
			tempo_samples, beat_alignment);
	end

	% gives the Viterbi algorithm a frame of observations (feature data), so that it
	% can compute the most likely state for the current frame.

	% PARAMETERS:
	% feature_data = matrix(feature_frame_length, this.num_features)
	%	is a 1xn cell matrix containing the raw features, or observations.
	%	There is a column in the cell matrix for each feature, and each cell is a
	%	(single dimensional) vector of length equal to the feature frame length.
	% tempo_alignment_estimates = cell(1, this.num_features)
	%	is a cell matrix containing a set of tempo and alignment estimates (no
	% 	confidences) for the current frame. There is one set for each feature,
	%	which means the cell matrix is 1 row (since it's only one frame's worth of
	%	data) by n columns, where n is the number of features. The tempo and alignment
	%	estimates should be in SAMPLES. These estimates will be used to generate the
	%	states that are searched for this frame


	function step_frame(this, feature_data, tempo_alignment_estimates)
		this.frame_number = this.frame_number + 1;

		% [OLD NOTES ON CLUSTERING]
		% these are the tempos that will be searched over
		% the second column represents the size of the cluster. We can
		% normalise this to give each cluster a weight.
		% -> extension: track which features had estimates in the cluster
		% containing the tempo that is finally picked as the most likely
		% one for this frame, and score those features higher.

		% for now, just use all the tempo/alignment estimates given as states,
		% even if they're close together.

		% find which states to search over, using the tempo/alignment estimates.

		max_tp_estimates_per_feature = ...
			this.params.max_tempo_peaks * this.params.max_alignment_peaks;

		% initialise states cell array

		% frame updates are O(n^2) in this number!!
		% Should be well under 1000 for speed, I think. 4000 takes minutes.
		% This number would be multiplied by other stuff too, if the model becomes
		% more complex.
		max_states_to_search = max_tp_estimates_per_feature*this.num_features;

		new_states = cell(max_states_to_search, 1);
		num_new_states = 0;

		% add all of the estimates
		for n = 1:this.num_features
			feature_n_estimates = tempo_alignment_estimates{n};
			for estimate_idx = 1:size(feature_n_estimates, 1)
				curr_tp_estimate = feature_n_estimates(estimate_idx, :);
				tempo_estimate = curr_tp_estimate(1);
				alignment_estimate = curr_tp_estimate(3);

				num_new_states = num_new_states + 1;
				new_states{num_new_states} = ...
					this.create_state(tempo_estimate, alignment_estimate);
			end
		end

		% Trim cell array down to size
		new_states = new_states(1:num_new_states);

		past_states = this.current_states;
		past_probabilities = this.current_probabilities;

		% in exceptional circumstances...
		if isempty(new_states)
			% reuse our past states again
			new_states = past_states;
		end

		new_probs = this.update_forward_message(past_states, past_probabilities, ...
			new_states, feature_data);

		% overwrite past information now we don't need them any more
		this.current_states = new_states;
		this.current_probabilities = new_probs;

		% Now find argmax of the probabilities
		% what if there's no most likely state?
		[winning_state, winning_prob] = this.most_likely_state(new_states, new_probs);
		this.winning_states{this.frame_number} = winning_state;
		this.winning_probabilities(this.frame_number) = winning_prob;

		% plot(this.current_probabilities);

	end


	% Calculate the beat times predicted by this algorithm, after frames have been
	% added and winning states determined.
	% Note that this method is inherently non-realtime, but it serves to demonstrate
	% the functionality. Predictions for a given time is still be made only
	% using past estimate data

	% Each state contains variables of tempo period and beat location. These
	% respectively represent the most likely tempo (period length, in seconds), and
	% the absolute time at which a beat most recently occurred (most probably), in
	% the frame for which the state was created.
	% Since predictions are causal, the beat location will always be a time before
	% its corresponding feature frame ends.
	% Predictions of future beat times are done by adding multiples of the tempo
	% period to the beat location, until the end of the next frame is reached.
	function beat_times = compute_beat_times(this)
		num_predictions = this.frame_number;

		% assume 3 beats predicted per row of estimate data,
		% just to be able to preallocate something
		beat_times = zeros(3*num_predictions, 1);
		num_beats_predicted = 0;

		for k = 1:num_predictions
			% output beat predictions using the kth winning tempo/alignment
			% estimate, up to the estimate time of the next frame
			estimate_time = this.params.estimate_time(k);
			% note the estimates for the final frame may extend slightly over the
			% end of the audio
			next_estimate_time = this.params.estimate_time(k+1);

			% these should be in seconds.
			kth_winning_state = this.winning_states{k};
			if isempty(kth_winning_state)
				warning('no winning state in frame %d (%.1f s)', k, estimate_time);
				% no beats for this frame?
				continue;
			end
			kth_tempo = kth_winning_state.tempo_period;
			kth_beat_time = kth_winning_state.beat_location;

			predicted_beat_time = kth_beat_time;

			% so we think there's a beat here, but this was some time in the past.
			% it may even be before the end of the previous feature frame.
			% we need to add multiples of the tempo period until we get to after the
			% end of the current feature frame time.

			while predicted_beat_time <= estimate_time
				predicted_beat_time = predicted_beat_time + kth_tempo;
			end

			% now output predictions as frame times
			% do we need to allow for latency?
			overlap_allowed = 0.01;
			while predicted_beat_time <= next_estimate_time + overlap_allowed
				num_beats_predicted = num_beats_predicted + 1;
				beat_times(num_beats_predicted) = predicted_beat_time;

				% try to predict another beat time
				predicted_beat_time = predicted_beat_time + kth_tempo;
			end
		end
		% trim beat times array to remove zeros?
		beat_times = beat_times(1:num_beats_predicted);
	end

	function output_beat_times(this, data_directory)
		filename = strcat(data_directory, '/', this.predictor_name, ...
			this.data_output_suffix);
		outfile = fopen(filename, 'w+');

		for beat_index = 1:length(this.beat_times)
			beat_time = this.beat_times(beat_index);
			if beat_index > 1 && beat_time ~= 0
				fprintf(outfile, '%.3f\n', this.beat_times(beat_index));
			end
		end
	end



	% generates all possible model states for the given range of tempos, in samples
	% PARAMETERS
	%	tempos: range of tempos to consider. Should be a column vector.
	%	beat_alignment_coarseness: generate a state for every n_th possible
	%	beat alignment. For all possible states, set this to 1. Can set to
	%	higher to prevent generating a large number of initial states.
	function states = generate_all_states(this, tempos, beat_alignment_coarseness)
		num_tempos = length(tempos);
		% vague estimate of number of states needed. Will increase if
		% more parameters are added to the model.
		num_states = round(mean(tempos)^2);

		states = cell(num_states, 1);
		num_states = 0;

		for tempo_idx = 1:num_tempos
			tempo = tempos(tempo_idx);
			% this is with reference to the end of the frame
			for beat_alignment = 0:beat_alignment_coarseness:tempo-1
				state = this.create_state(tempo, -1*beat_alignment);
				num_states = num_states + 1;
				states{num_states} = state;
			end
		end

		% trim returned lists to have length equal to the actual number of states
		states = states(1:num_states);
	end
end

methods (Static)

	% uses the prior distribution to generate an initial forward message
	% P(X_0 | null), i.e. the probability of any state given no information
	function probs = initial_forward_message(initial_states)
		probs = zeros(size(initial_states));

		for state_idx = 1:length(initial_states)
			state = initial_states{state_idx};
			probs(state_idx) = musical_model.prior_prob(state);
		end

		% make the whole thing sum to 1 - as some priors are only proportional
		probs = probs/sum(probs);
	end

	% Compute the forward message of the hidden Markov model, given the most recent
	% observations. This is a distribution over all possible states.
	% in other words, do one step of filtering.
	% the initial forward message is the initial prior.

	% PARAMETERS
	% old_states, old_probs = cell(num_old_states, 1), zeros(num_old_states, 1)
	% 	are vectors of states and a probability for each state respectively, given the
	% 	past observations. This need not contain all possible states, as due to
	% 	computational constraints, we calculate probabilities only for a subset of
	% 	all possible tempos. By definition, all other states have a probability of
	% 	zero, given the previous observations.
	% new_states = cell(num_new_states, 1)
	% 	is a list of the new states to calculate the probabilities for, which is a
	% 	subset of all possible states that is heuristically determined by the
	% 	tempo_alignment_estimator. Note that the new states don't have to be a subset of
	% 	the old states.
	% observations = matrix(feature_frame_length, num_features);
	% 	is the matrix of feature data from each feature, for the current
	% 	feature frame
	% tempos = cell(num_new_states, 1)
	% 	is a vector of states that will be part of the new forward message.
	% 	the new message
	function new_probs = update_forward_message(past_states, past_probabilities, ...
			new_states, new_observations)
		% new_forward_message(X_t)
		% = Prob(X_t | e_{1...t})
		% = K*(sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)
		%		*P(current_observations | X_t)
		% where K is a normalising constant, so that the sum of
		% probabilities for each X_t is 1.

		% calculate P(current_observations | X_t) for each X_t that we are
		% considering. Note that these are proportional: the probability over all
		% possible observations given some X_t may not sum to 1.
		observation_probs = bp_viterbi.compute_observation_probs(new_observations, ...
			new_states);

		% Calculate, for each new state X_t,
		% (sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)
		transition_probs = bp_viterbi.compute_transition_probs(past_states, ...
			past_probabilities, new_states);

		% now multiply them and make it sum to 1
		new_probs = observation_probs.*transition_probs;
		new_probs = new_probs/sum(new_probs);
	end

	function new_probs = compute_transition_probs(old_states, old_probs, new_states)
		num_new_states = size(new_states, 1);
		num_old_states = size(old_states, 1);

		% Calculate, for each new state X_t,
		% (sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)

		new_probs = zeros(num_new_states, 1);

		for new_state_idx = 1:num_new_states
			new_state = new_states{new_state_idx};
			new_state_prob = 0;

			% we sum the probabilities of going from each of the old states
			% to the (fixed) new state (t2, b2), ignoring new information
			for old_state_idx = 1:num_old_states
				old_state = old_states{old_state_idx};
				old_state_prob = old_probs(old_state_idx);

				% prob from old_state to new_state
				transition_prob = ...
						musical_model.transition_prob(old_state, new_state);
				new_state_prob = new_state_prob + transition_prob*old_state_prob;
			end
			new_probs(new_state_idx) = new_state_prob;
		end
		% normalise? -> no, as transition probabilities are already normalised.
	end



	% calculates observations for a set of states under a certain restriction of
	% which tempos to use.

	% PARAMETERS:
	% tempos = matrix(num_tempos, 1);
	% 	is a list of find which tempos in SAMPLES we need to search over.
	% observations = matrix(feature_window_length, num_features)
	% 	is a matrix of feature data, windowed for the current frame
	function probs = compute_observation_probs(observations, states)
		frame_length = size(observations, 1);
		num_features = size(observations, 2);
		num_states = size(states, 1);

		% first, calculate the autocorrelation function for each feature.
		% enhancement -> pass in the values for each tempo, from the tempo_alignment
		% estimator.

		acf_data = zeros(frame_length, num_features);

		for n = 1:num_features
			acf_data(:, n) = autocorrelation(observations(:, n));
		end

		% Now, for each of the tempos found in the list of state, calculate the
		% beat alignment function. To save recalculating the BAF, we create a map of
		% tempos to index of the beat alignment function table, so we can do lookups
		% only knowing the tempo value (as found in the state)
		num_tempos = 0;
		idx_for_tempo = containers.Map('KeyType', 'int32', 'ValueType', 'uint32');

		% record longest tempo just to preallocate size of beat alignment function
		% array. (since the beat alignment function is as long as the given tempo,
		% in samples. So clearly, they won't all be this long, just the slowest is.)
		longest_tempo_period = 0;

		% create the map
		for state_idx = 1:num_states
			state = states{state_idx};
			tempo_of_state = state.tempo_samples;
			if tempo_of_state > longest_tempo_period
				longest_tempo_period = tempo_of_state;
			end
			if ~isKey(idx_for_tempo, tempo_of_state)
				% add to map
				num_tempos = num_tempos + 1;
				idx_for_tempo(tempo_of_state) = num_tempos;
			end
		end

		baf_data = zeros(longest_tempo_period, num_tempos, num_features);

		% calculate beat alignment data.
		% Yes I know this is clumsy. If I implement clustering of tempo and beat 
		% alignment estimates, then this might not be needed any more, since we can
		% just pass in the value.
		all_tempos = idx_for_tempo.keys;
		for i = 1:num_tempos
			ith_tempo = all_tempos{i};
			feature_frame_n = observations(:, n);
			for n = 1:num_features
				baf_data(1:ith_tempo, idx_for_tempo(ith_tempo), n) = ...
					beat_alignment_function(feature_frame_n, ith_tempo);
			end
		end

		% initialise result vector
		probs = zeros(num_states, 1);

		for state_idx = 1:num_states
			curr_state = states{state_idx};
			curr_tempo = curr_state.tempo_samples;
			% to index into ACF
			acf_tempo_idx = curr_tempo + 1;
			% to index into BAF column
			baf_tempo_idx = idx_for_tempo(curr_state.tempo_samples);
			% correct for 1-indexing, and that beat alignment is negative
			beat_align_idx = -1*curr_state.beat_alignment + 1;

			% -> currently we just add up all features equally!
			observation_prob = 0;
			for n = 1:num_features
				% this is our model assumption. 'Probability' here really means a
				% proportional measure; the value may not be normalised.
				acf_at_tempo = acf_data(acf_tempo_idx, n);

				% ideas
% 				acf_at_double_tempo = 0.5*acf_data(floor(curr_tempo/2) + 1, n) + ...
% 					0.5*acf_data(ceil(curr_tempo/2) + 1, n);
% 				% if the tempo is fast enough to have half of its value in the
% 				% acf, increase it by this much
% 				if 2*curr_tempo <= frame_length/2
% 					acf_at_half_tempo = acf_data(2*curr_tempo + 1, n);
% 				else
% 					acf_at_half_tempo = 1;
% 				end

				baf_at_tempo_and_alignment = baf_data(beat_align_idx, baf_tempo_idx, n);
				observation_prob_for_feature = acf_at_tempo*baf_at_tempo_and_alignment;
				% times feature_weight(n)?

				observation_prob = observation_prob + observation_prob_for_feature;
			end

			probs(state_idx) = observation_prob;
		end
	end


	function [state, prob] = most_likely_state(states, probs)
		if isempty(states)
			state = {};
			return
		end
		if ~isequal(size(states), size(probs))
			error('probability and state vectors have different sizes!');
		end

		highest_state_probability = 0;
		% arbitrary
		most_likely_state = states{1};
		num_states = size(states, 1);

		for state_idx = 1:num_states
			% this is, in theory, the probability of this state being the actual one,
			% given all observations so far.
			state = states{state_idx};
			state_probability = probs(state_idx);

			if state_probability > highest_state_probability
				most_likely_state = state;
				highest_state_probability = state_probability;
			end
		end

		state = most_likely_state;
		prob = highest_state_probability;
	end

end % methods (Static)

end

