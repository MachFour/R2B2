bp = bp_viterbi;
bp.initialise(p, 'bp_viterbi', 4);
new_tempos = (100:10:200)';
bp.update_forward_message(abs(randn(1024, 8)), new_tempos);
