import numpy as np
from autodeer.hardware import Interface,Waveform,Sequence,SequenceTable
from autodeer.ResPro import resonatorProfile
import autodeer.hardware.pulses as pulses
import autodeer.tools as tools


def sequence_nutation(awg:Interface,p_start:int,p_step:int,nx:int,h:int,IF:float=2):

    awg.Abort()
    awg.deleteAllWaveforms(ch=3)
    awg.resetSequenceTable()

    seq_table = SequenceTable(awg)
    seq_table.table = []
    
    ps_length_table = np.linspace(p_start,p_start + nx * p_step,nx)
    
    for iD,length in enumerate(ps_length_table):
        iD = iD + 1
        wave = Waveform(awg,
                        pulses.rectPulse(12e9,length*1e-9,IF*1e9,1,0),
                        None)
        wave.pulse_blanking()
        wave.enforce_gradualrity()
        wave.set_amp('max')

        length = wave.num_points

        wave.define_new_waveform(iD,ch=3)
        wave.send_waveform(iD,ch=3)
        
        # Build the sequences

        seq = Sequence(awg)
        if iD == 1: # The first element of a sequence is different to the subsequent ones
            seq.encode(
                        {"data_cmd":0,
                    "init_marker_seq":1,
                    "marker_enable":1,
                    "seg_adv_mode":3,
                    "seq_loop_cnt":1024,
                    "seg_loop_cnt":h,
                    "seg_id":iD,
                    "seg_start_offset":0,
                    "seg_end_offset":length-1}
                    )
        elif iD == nx:
            seq.encode(
                        {"data_cmd":0,
                    "init_marker_seq":0,
                    "end_marker_seq":1,
                    "end_marker_sce":1,
                    "marker_enable":1,
                    "seg_adv_mode":3,
                    "seq_loop_cnt":1024,
                    "seg_loop_cnt":h,
                    "seg_id":iD,
                    "seg_start_offset":0,
                    "seg_end_offset":length-1}
                    )            
        else:
            seq.encode(
                        {"data_cmd":0,
                    "marker_enable":1,
                    "seg_adv_mode":3,
                    "seq_loop_cnt":1024,
                    "seg_loop_cnt":h,
                    "seg_id":iD,
                    "seg_start_offset":0,
                    "seg_end_offset":length-1}
                    )


        # Append Sequences to the sequnece table
        seq_table.table.append(seq)
        tools.progress_bar_frac(iD,nx)
    seq_table.write_to_str()
    awg.setSequenceTable(seq_table,0,3)

    pass



def deer_pulse_5p(awg:Interface,tpump:int,nu,res_prof:resonatorProfile=None):

    awg.Abort()
    awg.deleteAllWaveforms(ch=3)
    awg.resetSequenceTable()

    seq_table = SequenceTable(awg)
    seq_table.table = []
    
    # Positive Hyperbolic Secant Pulse
    [p_real,p_imag] = pulses.HSorder1(12,
                                tpump,
                                nu,
                                resonator=res_prof,
                                HSbeta=10,
                                HSorder=[1,6],
                                phase=0,
                                scale=1)

    wave1 = Waveform(awg,np.vstack([p_real,p_imag]))
    wave1.pulse_blanking()
    wave1.enforce_gradualrity()
    wave1.set_amp('max')
    wave1._awg_check()
    wave1.define_new_waveform(1,ch=3)
    wave1.send_waveform(1,ch=3)


    # Negative Hyperbolic Secant Pulse

    [n_real,n_imag] = pulses.HSorder1(12,
                                tpump,
                                np.flip(nu),
                                resonator=res_prof,
                                HSbeta=10,
                                HSorder=np.flip([1,6]),
                                phase=0,
                                scale=1)
    wave2 = Waveform(awg,np.vstack([n_real,n_imag]))
    wave2.pulse_blanking()
    wave2.enforce_gradualrity()
    wave2.set_amp('max')
    wave2._awg_check()
    wave2.define_new_waveform(2,ch=3)
    wave2.send_waveform(2,ch=3)

    # Build the sequence table
    seq1 = Sequence(awg)
    seq1.encode(
        {"data_cmd":0,
    "init_marker_seq":1,
    "marker_enable":1,
    "seg_adv_mode":3,
    "seq_loop_cnt":1,
    "seg_loop_cnt":1,
    "seg_id":1,
    "seg_start_offset":0,
    "seg_end_offset":wave1.num_points - 1}
    )
    seq2 = Sequence(awg)
    seq2.encode(
        {"data_cmd":0,
    "init_marker_seq":0,
    "end_marker_seq":1,
    "end_marker_sce":1,
    "marker_enable":1,
    "seg_adv_mode":3,
    "seq_loop_cnt":1,
    "seg_loop_cnt":1,
    "seg_id":2,
    "seg_start_offset":0,
    "seg_end_offset":wave2.num_points - 1}
    )

    seq_table = SequenceTable(awg)
    seq_table.table = [seq1,seq2]
    seq_table.write_to_str()

    seq_table.write_to_str()
    awg.setSequenceTable(seq_table,0,3)

    return [p_real+1j*p_imag,n_real+1j*n_imag]
