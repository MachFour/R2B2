% class used to contain the state of a frame (in terms of tempo, phase, etc.)
% for use with the bp_viterbi class
% see musical_model.m for a description of the following model

classdef model_state

properties
	% tempo period (in seconds): i.e. time between beats
	tempo_period;
	% time when the latest beat was in the current frame.
	% measured in absolute time (i.e. with respect to 0 seconds of audio)
	beat_location;

end % properties

properties (Dependent)
	tempo_bpm
end

methods
	% constructor
	function state = model_state(tempo_period, beat_location)
		if nargin == 0
			warning('model_state class constructed with default values');
			state.tempo_period = 0;
			state.beat_location = 0;
		else
			state.tempo_period = tempo_period;
			state.beat_location = beat_location;
		end
	end

	function t = get.tempo_bpm(this)
		t = 60/this.tempo_period;
	end

end % methods

end
