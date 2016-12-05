% bp_viterbi.m

% Class implementing the Viterbi Algorithm,
% used to do prediction of beat times.

% Author: Max Fisher

classdef bp_viterbi < handle

properties
	% counts the frame number
	frame_idx;

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

	winning_states;

	data_output_suffix = '-beat-times.txt';

	name;
	num_features;

	% processing_params object, used for sample to time conversions.
	params;

end

methods
	function initialise(this, params, predictor_name, num_features)
		this.params = params;
		this.predictor_name = predictor_name;
		this.num_features = num_features;

		this.frame_idx = 1;


		% restrict to 50-180 BPM initially - low bpm makes it really
		% inefficient!!! (there are a huge number of states)
		min_initial_bpm = 45;
		max_initial_bpm = 180;

		sample_to_bpm_factor = this.params.feature_sample_rate*60;
		min_tempo_samples = round(sample_to_bpm_factor/max_initial_bpm);
		max_tempo_samples = round(sample_to_bpm_factor/min_initial_bpm);

		% choose a sparser set of tempos to start with
		initial_tempos = (min_tempo_samples:2:max_tempo_samples)';
		this.generate_initial_forward_message(initial_tempos);

		% allocate some space for winning estimates;
		this.winning_states = cell(100, 1);

	end

	% uses the prior distribution to generate an initial forward message
	% P(X_0 | null), i.e. the probability of any state given no information
	function generate_initial_forward_message(this, initial_tempos)
		states = this.generate_all_states(initial_tempos);
		probs = zeros(size(states));

		for state_idx = 1:length(states)
			state = states{state_idx};
			probs(state_idx) = musical_model.prior_prob(state);
		end

		% make the whole thing sum to 1 - as some priors are only proportional
		probs = probs/sum(probs);
		this.current_states = states;
		this.current_probabilities = probs;
	end

	% converts tempo and phase in samples to a model state which measures
	function s = state_from_tempo_and_alignment(this, tempo_samples, beat_alignment)
		feature_sample_rate = this.params.feature_sample_rate;
		frame_end_time = this.params.prediction_time(this.frame_idx);

		% calculate tempo and absolute beat location in seconds
		tempo_period = tempo_samples/feature_sample_rate;
		beat_location = frame_end_time - beat_alignment/feature_sample_rate;
		s = model_state(tempo_period, tempo_samples, beat_location, beat_alignment);
	end

	% gives the Viterbi algorithm a frame of observations (feature data), so that it
	% can compute the most likely state for the current frame.

	% PARAMETERS:
	% tempo_estimates = cell(1, this.num_features)
	%	  is a cell matrix containing a set of tempo estimates (no confidences)
	%	  for the current frame. There is one set for each feature, which means the cell
	%	  matrix is 1 row (since it's only one frame's worth of data) by n columns, where
	%	  n is the number of features. The tempo estimates should be in SAMPLES.

	% feature_data = cell(1, this.num_features)
	%	  is a 1xn cell matrix containing the raw features, or observations.
	%	  There is a column in the cell matrix for each feature, and each cell is a
	%	  (single dimensional) vector of length equal to the feature frame length.

	function step_frame(this, feature_data, tempo_estimates)
		% these are the tempos that will be searched over
		% the second column represents the size of the cluster. We can
		% normalise this to give each cluster a weight.

		% -> extension: track which features had estimates in the cluster
		% containing the tempo that is finally picked as the most likely
		% one for this frame, and score those features higher.
		tempo_clusters = this.cluster_tempo_estimates(tempo_estimates);

		% for the moment we won't use which features had estimates in which cluster;
		% just extract the actual clusters that were identified.
		num_clusters = size(tempo_clusters, 1);
		tempos_to_search = zeros(num_clusters, 1);
		for i = 1:num_clusters
			% clusters are in the first column of the cell array
			tempos_to_search(i) = tempo_clusters{i, 1};
		end

		this.update_forward_message(feature_data, tempos_to_search);

		% find the most likely state
		most_likely_state = {};
		highest_state_probability = 0;
		num_states = size(this.current_states, 1);

		for state_idx = 1:num_states
			% this is, in theory, the probability of this state being the actual one,
			% given all observations so far.
			state = this.current_states(state_idx);
			state_probability = this.current_probabilities(state_idx);

			if state_probability > highest_state_probability
				most_likely_state = state;
				highest_state_probability = state_probability;
			end
		end

		this.winning_states(state_idx) = most_likely_state;
	end

	% Narrows down the search of possible tempos to only a small
	% subset identified by the tempo/phase estimator as likely.
	% Different features will have picked different peaks for their autocorrelation
	% functions, but if the different estimates are close enough, we treat them as
	% both estimating the same tempo.  The idea is that if a lot of different
	% features agree on the tempo, then it's more likely to be that tempo.

	% PARAMETERS:
	% tempo_estimates = cell(1, this.num_features)
	%	  is a cell matrix containing a set of tempo estimates (no confidences)
	%	  for the current frame. There is one set for each feature, which means the cell
	%	  matrix is 1 row (since it's only one frame's worth of data) by n columns, where
	%	  n is the number of features. The tempo estimates should be in SAMPLES.
	function clustered_tempos = cluster_tempo_estimates(this, tempo_estimates)

		% stores the set of tempos to search over in this iteration of the
		% Viterbi algorithm
		% each tempo is tagged with the number of the feature that it comes
		% from.
		candidate_tempos = zeros(this.num_features*this.params.max_tempo_peaks, 2);
		candidate_tempo_idx = 1;

		for feature_idx = 1:this.num_features
			tempo_estimates_list = tempo_estimates{feature_idx};
			for estimate_idx = 1:length(tempo_estimates_list)
				tempo_estimate = tempo_estimates_list(tempo_estimate);
				candidate_tempos(candidate_tempo_idx, 1) = tempo_estimate;
				candidate_tempos(candidate_tempo_idx, 2) = feature_idx;

				candidate_tempo_idx = candidate_tempo_idx + 1;
			end
		end
		% trim list to its actual length
		candidate_tempos = candidate_tempos(1:candidate_tempo_idx, :);

		% now group together similar tempo estimates;
		% use mean shift clustering *BY BPM*.
		% 2BPM seems like a reasonable
		% cluster width

		sample_to_bpm_factor = this.params.feature_sample_rate*60;

		bpm_distance = @(x, y) sqrt(sum(sample_to_bpm_factor.*(1./x - 1./y)).^2);
		cluster_width = 2; % BPM
		clustered_tempos = mean_shift_cluster(candidate_tempos, ...
			bpm_distance, cluster_width);
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
		num_predictions = size(this.winning_states, 1);

		% assume 3 beats predicted per row of estimate data,
		% just to be able to preallocate something
		beat_times = zeros(3*num_predictions, 1);
		beat_index = 1;

		for k = 1:num_predictions
			% output beat predictions using the kth winning tempo/phase
			% estimate, up to the estimate time of the next frame
			estimate_time = this.params.prediction_time(k);
			% note the estimates for the final frame may extend slightly over the
			% end of the audio
			next_estimate_time = this.params.prediction_time(k+1);
			% this doesn't mean much if k = 1 but it won't cause any problems
			prev_estimate_time = this.params.prediction_time(k-1);

			% these should be in seconds.
			kth_winning_state = this.winning_states{k};
			kth_tempo = kth_winning_state.tempo_period;
			kth_beat_time = kth_winning_state.beat_location;

			if prev_estimate_time > kth_beat_time || kth_beat_time > estimate_time
				warning(strcat('Winning beat time estimate for frame %d', ...
					'occurs before the end of frame %d!'), k, k-1);
			end

			predicted_beat_time = kth_beat_time + kth_tempo;
			% do we need to allow for latency?
			overlap_allowed = 0.01;
			while predicted_beat_time <= next_estimate_time + overlap_allowed
				beat_times(beat_index) = predicted_beat_time;
				beat_index = beat_index + 1;
				% try to predict another beat time
				predicted_beat_time = predicted_beat_time + kth_tempo;
			end
		end
		% trim beat times array to remove zeros?
		beat_times = beat_times(1:beat_index);
	end

	function output_beat_times(this, data_directory)
		filename = strcat(data_directory, '/', this.predictor_name, ...
			this.DATA_OUTPUT_SUFFIX);
		outfile = fopen(filename, 'w+');

		for beat_index = 1:length(this.beat_times)
			beat_time = this.beat_times(beat_index);
			if beat_index > 1 && beat_time ~= 0
				fprintf(outfile, '%.3f\n', this.beat_times(beat_index));
			end
		end
	end


	% need an Observation and a State class?

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
	% 	tempo_phase_estimator. Note that the new states don't have to be a subset of
	% 	the old states.
	% observations = matrix(feature_frame_length, num_features);
	% 	is the matrix of feature data from each feature, for the current
	% 	feature frame
	% tempos = matrix(num_tempos, 1)
	% 	is the set of tempos (in samples) to consider when generating the states for
	% 	the new message
	function update_forward_message(this, observations, tempos)
		% new_forward_message(X_t)
		% = Prob(X_t | e_{1...t})
		% = K*(sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)
		%		*P(current_observations | X_t)
		% where K is a normalising constant, so that the sum of
		% probabilities for each X_t is 1.

		this.current_tempos = tempos;
		new_states = this.generate_all_states(tempos);

		% calculate P(current_observations | X_t) for each X_t that we are
		% considering. Note that these are proportional: the probability over all
		% possible observations given some X_t may not sum to 1.
		observation_probs = this.compute_observation_probs(observations, ...
			new_states, tempos);

		% Calculate, for each new state X_t,
		% (sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)
		transition_probs = this.compute_transition_probs(this.current_states, ...
			this.current_probabilities, new_states);

		% now multiply them and make it sum to 1
		new_probs = observation_probs.*transition_probs;
		new_probs = new_probs/sum(new_probs);

		this.current_states = new_states;
		this.current_probabilites = new_probs;

	end


	% generates all possible model states for the given range of tempos, in samples
	% XXX OBSERVATION PROBABILITY CALCULATIONS DEPEND ON CORRECT ORDERING OF THE STATES
	% DO NOT CHANGE THIS FUNCTION WITHOUT FIXING THIS
	function states = generate_all_states(this, tempos)
		num_tempos = length(tempos);
		% conservative estimate of number of states needed. Will increase if
		% more parameters are added to the model.
		num_states = num_tempos^2;

		states = cell(num_states, 1);
		num_states = 0;

		for tempo_idx = 1:length(tempos)
			tempo = tempos(tempo_idx);
			% this is with reference to the end of the frame
			for beat_alignment = 0:tempo-1
				state = this.state_from_tempo_and_alignment(tempo, beat_alignment);
				num_states = num_states + 1;
				states{num_states} = state;
			end
		end

		% trim returned lists to have length equal to the actual number of states
		states = states(1:num_states);
	end
end

% transition and observation probability calculation
methods (Static)
	function new_probs = compute_transition_probs(old_states, old_probs, new_states)
		num_new_states = size(new_states, 1);
		num_old_states = size(old_states, 1);

		% Calculate, for each new state X_t,
		% (sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)

		new_probs = zeros(num_new_states, 1);

		for new_state_idx = 1:num_new_states
			new_state = new_states{new_state_idx};

			% we sum the probabilities of going from each of the old states
			% to the (fixed) new state (t2, b2), ignoring new information
			for old_state_idx = 1:num_old_states
				old_state = old_states{old_state_idx};
				old_state_prob = old_probs(old_state_idx);

				% prob from old_state to new_state
				old_new_prob = musical_model.transition_prob(old_state, new_state);

				new_probs(new_state_idx) = new_probs(new_state_idx) + ...
					old_state_prob*old_new_prob;
			end
		end
		% normalise???
	end



	% calculates observations for all possible states under a certain restriction of
	% which tempos to use.

	% PARAMETERS:
	% tempos = matrix(num_tempos, 1);
	% 	is a list of find which tempos in SAMPLES we need to search over.
	% observations = matrix(feature_window_length, num_features)
	% 	is a matrix of feature data, windowed for the current frame
	function probs = compute_observation_probs(observations, states, tempos)
		num_tempos = size(tempos, 1);
		num_states = size(states, 1);

		probs = zeros(num_states, 1);

		frame_length = size(observations, 1);
		num_features = size(observations, 2);

		% first, calculate the autocorrelation function for each feature, and the
		% beat alignment function for each of the valid tempos.

		% currently the autocorrelation function returns a vector of half the input
		% length.
		acf_length = frame_length/2;
		acf_data = zeros(acf_length, num_features);

		% the beat alignment function is as long as the given tempo (in samples)
		baf_length = max(tempos);
		baf_data = zeros(baf_length, num_tempos, num_features);

		% calculate autocorrelation and beat alignment data;
		for n = 1:num_features
			feature_frame_n = observations(:, n);
			frame_n_acf = autocorrelation(feature_frame_n);
			acf_data(:, n) = frame_n_acf;
			for tempo_idx = 1:num_tempos
				tempo_for_baf = tempos(tempo_idx);
				baf_data(1:tempo_for_baf, tempo_idx, n) = ...
					beat_alignment_function(feature_frame_n, tempo_for_baf);
			end
		end

		% So that we can iterate over the generated states, we create a map of the
		% set of tempos to their index in the input array
		tempo_idx_map = containers.Map(tempos', 1:num_tempos);

		for state_idx = 1:num_states
			curr_state = states{state_idx};
			state_tempo = curr_state.tempo_samples;
			% this is with reference to the end of the frame
			state_beat_alignment = state.beat_alignment;

			tempo_idx = tempo_idx_map(state_tempo);
			% correct for 1-indexing
			beat_align_idx = state_beat_alignment + 1;

			% -> currently we just add up all features equally!
			state_observation_prob = 0;
			for n = 1:num_features
				% this is our model assumption. 'Probability' here really means a
				% proportional measure; the value may not be normalised.
				feature_observation_prob = acf_data(tempo_idx, n)*...
					baf_data(beat_align_idx, tempo_idx, n);
					% times feature_weight(n)?
				state_observation_prob = state_observation_prob + ...
					feature_observation_prob;
			end

			probs(state_idx) = state_observation_prob;
		end
	end

end % methods (Static)

end

