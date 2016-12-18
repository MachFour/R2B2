
% viterbi_helper.m

% Class of static helper functions for the Viterbi Algorithm
% Author: Max Fisher

classdef viterbi_helper < handle

methods (Static)
	% plots a cell array of states, showing for each one, the observed beat location,
	% and the next predicted beat time given the state's tempo
	function plot_states(states)
		num_states = length(states);
		if num_states == 0
			return
		end

		frame_end_time = states{1}.frame_end_time;
		next_frame_end_time = frame_end_time + states{1}.params.feature_hop_size/...
			states{1}.params.feature_sample_rate;
		prev_frame_end_time = frame_end_time + (frame_end_time - next_frame_end_time);

		figure;
		hold on;
		% show line at end time of frame
		stem(frame_end_time, num_states, 'gx');
		stem(next_frame_end_time, num_states, 'rx');
		stem(prev_frame_end_time, num_states, 'rx');
		% do this so that the points with the same index come out in the right colour
		for state_idx = 1:length(states)
			state = states{state_idx};

			predicted_locs = state.beat_location:state.tempo_period:next_frame_end_time;
			scatter(predicted_locs, state_idx*ones(size(predicted_locs)), 'kx');

			observed_locs = state.beat_location:-state.tempo_period:prev_frame_end_time;
			scatter(observed_locs, state_idx*ones(size(observed_locs)), 'ko');
			% add a line to join the dots
			plot([min(observed_locs), max(predicted_locs)], state_idx*[1, 1], 'b');
		end
		xlabel(sprintf('Plot of states for frame ending at %.1f s', frame_end_time));
		xlabel('Time (s)');
		ylabel('State index');
	end



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
		observation_probs = viterbi_helper.compute_observation_probs(new_observations, ...
			new_states);

		% Calculate, for each new state X_t,
		% (sum over all states x_t-1) P(X_t | x_t-1)*current_forward_message(x_t-1)
		transition_probs = viterbi_helper.compute_transition_probs(past_states, ...
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
		% NO -> This should be a maximum over x_t, not a sum

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
				transition_prob = old_state_prob * ...
						musical_model.transition_prob(old_state, new_state);

				if transition_prob > new_state_prob
					new_state_prob = transition_prob;
				end
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
		probs = zeros(num_states, 1);

		% calculate autocorrelation function of the current frame for each feature.
		% autocorrelation function returns a vector of half the frame_length
		acf_length = frame_length/2;
		acf_data = zeros(acf_length, num_features);

		for n = 1:num_features
			feature_frame_n = observations(:, n);
			acf = autocorrelation(feature_frame_n);
			acf_data(:, n) = acf;
		end

		% calculate the beat alignment data for each unique tempo among the states

		max_tempo_lag = frame_length/2 - 1;
		baf_data = zeros(max_tempo_lag + 1, max_tempo_lag + 1, num_features);
		baf_calculated = zeros(max_tempo_lag + 1, 1);

		for state_idx = 1:num_states
			state = states{state_idx};
			tempo = state.tempo_samples;
			% convert to 1-indexing
			tempo_idx = tempo + 1;
			if ~baf_calculated(tempo_idx)
				for n = 1:num_features
					feature_frame_n = observations(:, n);
					baf = beat_alignment_function(feature_frame_n, tempo);
					baf_data(1:tempo, tempo_idx, n) = baf;

				end
				baf_calculated(tempo_idx) = 1;
			end
		end

		% now calculate the observation probabilities

		for state_idx = 1:num_states
			state = states{state_idx};
			tempo = state.tempo_samples;
			tempo_idx = tempo + 1;
			% also correct for beat alignment being negative
			baf_idx = -1*state.beat_alignment + 1;

			% -> currently we just add up all features equally!
			observation_prob = 0;

			for n = 1:num_features
				% this is our model assumption. 'Probability' here really means a
				% proportional measure; the value may not be normalised.
				acf_height = acf_data(tempo_idx, n);

				% ideas
% 				acf_at_double_tempo = 0.5*acf_data(floor(tempo/2) + 1, n) + ...
% 					0.5*acf_data(ceil(tempo/2) + 1, n);
% 				% if the tempo is fast enough to have half of its value in the
% 				% acf, increase it by this much
% 				if 2*tempo <= frame_length/2
% 					acf_at_half_tempo = acf_data(2*tempo + 1, n);
% 				else
% 					acf_at_half_tempo = 1;
% 				end

				baf_height = baf_data(baf_idx, tempo_idx, n);

				observation_prob_for_feature = acf_height*baf_height;
				% + baf at half_tempo_and_alignment?
				% (if generated by that state)

				if 2*tempo + 1 < frame_length/2;
					baf_height_half_tempo = beat_alignment_function(...
						feature_frame_n, tempo*2);
				end

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

	% how to project the beat alignment of a past state into an equivalent beat
	% alignment several frames later, with beats occuring at the given tempo.
	% Everything is measured in samples.
	function projected_alignment = project_beat_alignment(old_alignment, tempo, ...
			frame_hops, hop_size)
		if frame_hops < 0
			error('number of frame hops has to be positive');
		end

		projected_alignment = mod(old_alignment - frame_hops*hop_size, -tempo);
	end

end % methods (Static)

end

