bp = bp_viterbi;
p = processing_params;
bp.initialise(p, 'bp_viterbi', 8);
new_tempos = (100:10:200)';
new_states = bp.generate_all_states(new_tempos, 8);
bp.update_forward_message(bp.current_states, bp.current_probabilities, ...
	new_states, abs(randn(1024, 8)));
