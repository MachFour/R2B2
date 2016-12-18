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
		initial_states = this.generate_initial_states;
		this.current_states = initial_states;
		this.current_probabilities = viterbi_helper.initial_forward_message(initial_states);

		% allocate some space for winning estimates and probabilities;
		this.winning_states = cell(50, 1);
		this.winning_probabilities = zeros(50, 1);

	end

	function states = generate_initial_states(this)
		% restrict to 50-180 BPM initially - low bpm makes it really
		% inefficient!!! (there are a huge number of states)
		min_initial_bpm = 45;
		max_initial_bpm = 160;

		sample_to_bpm_factor = this.params.feature_sample_rate*60;
		min_tempo_samples = round(sample_to_bpm_factor/max_initial_bpm);
		max_tempo_samples = round(sample_to_bpm_factor/min_initial_bpm);

		% choose a sparser set of tempos to start with
		initial_tempos = (min_tempo_samples:3:max_tempo_samples)';
		beat_coarseness = 10;
		states = this.generate_all_states(initial_tempos, beat_coarseness);
	end



	% converts tempo and beat alignment in samples to a model state which
	% measures. Note that beat alignment is negative
	function s = create_state(this, tempo_samples, beat_alignment)
		s = model_state(this.params, this.frame_number, tempo_samples, beat_alignment);
	end

	function s = copy_past_state(this, past_state)
		if past_state.frame_number > this.frame_number
			error('past state has a larger frame number than current frame number');
		end
		tempo = past_state.tempo_samples;
		old_beat_alignment = past_state.beat_alignment;
		frame_difference = this.frame_number - past_state.frame_number;
		frame_hop_size = this.params.feature_hop_size;

		beat_alignment = viterbi_helper.project_beat_alignment(...
			old_beat_alignment, tempo, frame_difference, frame_hop_size);
		s = model_state(this.params, this.frame_number, tempo, beat_alignment);
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

		%testing: add initial states
		%initial_states = this.generate_initial_states;

		% initialise states cell array

		% frame updates are O(n^2) in this number!!
		% Should be well under 1000 for speed, I think. 4000 takes minutes.
		% This number would be multiplied by other stuff too, if the model becomes
		% more complex.
		max_states_to_search = max_tp_estimates_per_feature*this.num_features;

		new_states = cell(max_states_to_search, 1);
		num_new_states = 0;

		%for state_idx = 1:length(initial_states)
		%	num_new_states = num_new_states + 1;
		%	new_states{num_new_states} = initial_states{state_idx};
		%end

		% add all of the estimates
		for n = 1:size(tempo_alignment_estimates, 2)
			feature_n_estimates = tempo_alignment_estimates{n};
			for estimate_idx = 1:size(feature_n_estimates, 1)
				curr_tp_estimate = feature_n_estimates(estimate_idx, :);
				tempo_estimate = curr_tp_estimate(1);
				alignment_estimate = curr_tp_estimate(3);

				num_new_states = num_new_states + 1;
				new_states{num_new_states} = ...
					this.create_state(tempo_estimate, alignment_estimate);
				% add tempo harmonics as estimates?
			end
		end

		% Trim cell array down to size
		new_states = new_states(1:num_new_states);

		% add the previously winning state (or top two>?)
		if this.frame_number > 1
			prev_winning_state = this.winning_states{this.frame_number - 1};
			new_states{end} = this.copy_past_state(prev_winning_state);
		end

		past_states = this.current_states;
		past_probabilities = this.current_probabilities;

		% in exceptional circumstances...
		if isempty(new_states)
			% reset to initial states again?
			new_states = this.generate_initial_states;
		end

		new_probs = viterbi_helper.update_forward_message(past_states, ...
			past_probabilities, new_states, feature_data);

		% overwrite past information now we don't need them any more
		this.current_states = new_states;
		this.current_probabilities = new_probs;

		% Now find argmax of the probabilities
		% what if there's no most likely state?
		[winning_state, winning_prob] = viterbi_helper.most_likely_state(new_states, new_probs);
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

end

