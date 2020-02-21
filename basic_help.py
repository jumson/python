#!/opt/conda/bin python
# 
##################################################
# Some functions I use frequently
# 
# Author: Jon Munson
# Description: Helper functions and classes
# 
##################################################
 
import numpy as np
from gnuradio import blocks
from gnuradio import digital
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
 
import threading
class AsyncWrite(threading.Thread): 
 
    def __init__(self, out):
 
        # calling superclass init
        threading.Thread.__init__(self) 
        self.text = ''
        self.out = out
 
    def log(self, text):
        self.text = text
        f = open(self.out, "a+")
        f.write(self.text)
        f.close()
 
        #print("Finished background file write to", self.out)
        # usage:     background = AsyncWrite('stuff to write\n', filename)
                #    background.start()
        # usage:     l = AsyncWrite(filename)
                #    l.log('log this \n')
 
class Slices:
    def __init__(self, 
                 symbols=[], 
                 ones_long=0, 
                 ones_short=0, 
                 zeros_long=0, 
                 zeros_short=0,
                ones_threshold=0,
                zeros_threshold=0,
                transitions=[],
                burst=[],
                 burst_ends = [],
                samp_rate=400000,
                min_width=20):
        #this will be an array of 0's and 1's
        self.symbols = symbols
        #these are the lengths of each kind of pulse, in microseconds
        self.ones_long = ones_long
        self.ones_short = ones_short
        self.zeros_long = zeros_long
        self.zeros_short = zeros_short
        # these are the values used t odiscriminate between longs and shorts
        self.ones_threshold = ones_threshold
        self.zeros_threshold = zeros_threshold
        # these are the positions in the slices where there are transitions
        self.transitions = transitions
        # this is the boolean burst itself, starts and ends with Trues
        # usually created with a sliced file: dat_sliced = np.fromfile(infile, dtype="float32")
                                            # my_slices.burst = np.array(dat_sliced,np.bool)
        self.burst = burst
        # The sample rate at which the slices/samples were processed
        self.samp_rate = samp_rate
        # the amount of samples in a row to suppose a legitimate pulse
        self.min_width = min_width
        self.burst_ends = burst_ends
 
class Burst:
    def __init__(self, 
                 name='',
                 modulation='',
                 encoding='',
                 slices=Slices(),
                 samp_rate=8e6,
                 center_freq=0,
                 working_samp_rate=4e5,
                 offset=0,
                 threshold=.5,
                 filter_cutoff=100e3,
                 filter_transition=10e3,
                 raw_input_file='',
                 demodulated_output='',
                 data_arr=[]):
        self.name = name
        self.modulation = modulation
        self.encoding = encoding
        self.slices = slices
        self.samp_rate = samp_rate
        self.working_samp_rate = working_samp_rate
        self.offset = offset
        self.threshold = threshold
        self.filter_cutoff = filter_cutoff
        self.filter_transition = filter_transition
        self.raw_input_file = raw_input_file
        self.demodulated_output = demodulated_output
        self.center_freq = center_freq
        self.data_arr = data_arr
 
    def load_data(self):
        try:
            self.data_arr = np.fromfile(self.raw_input_file, dtype="complex64")
        except:
            print("ensure Burst.raw_input_file has a legit sample path/file")
    
    def duration(self):
        if len(self.data_arr) == 0:
            self.load_data()
        print(len(self.data_arr)/self.samp_rate,'seconds')
 
class Colors:
    # Foreground:
    MAGENTA = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE = '\033[97m'
    BLACK = '\033[90m'
    # Formatting
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'    
    # End colored text
    END = '\033[0m'
    NC ='\x1b[0m' # No Color
    DEFAULT = '\033[99m'
 
c = Colors()
 
def filter_raw(burst_object):
    ### OOK_Demod  script
    ook_slicer = gr.top_block() # Define the container
 
            ##################################################
            # Variables
            ##################################################
    ook_slicer.working_samp_rate = working_samp_rate = burst_object.working_samp_rate
    ook_slicer.samp_rate = samp_rate = burst_object.samp_rate
    ook_slicer.threshold = threshold = burst_object.threshold
    ook_slicer.offset = offset = burst_object.offset
    ook_slicer.filter_transition = filter_transition = burst_object.filter_transition
    ook_slicer.filter_decimation = filter_decimation = int(samp_rate/working_samp_rate)
    ook_slicer.filter_cutoff = filter_cutoff = burst_object.filter_cutoff
    ook_slicer.raw_data_file = raw_data_file = burst_object.raw_input_file
    ook_slicer.filename_filtered = filename_filtered = raw_data_file+'.filtered'
    
 
            ##################################################
            # Blocks
            ##################################################
 
    #ook_slicer.throttle = blocks.throttle(gr.sizeof_gr_complex*1, working_samp_rate,True)
    #ook_slicer.multiply_constant = blocks.multiply_const_vcc((4, ))
    ook_slicer.freq_xlating_fir_filter = filter.freq_xlating_fir_filter_ccc(filter_decimation, (firdes.low_pass(1, samp_rate, filter_cutoff, filter_transition)), offset, samp_rate)
    ook_slicer.file_source = blocks.file_source(gr.sizeof_gr_complex*1, raw_data_file, False)
    ook_slicer.file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, filename_filtered, False)
    ook_slicer.file_sink_0.set_unbuffered(False)
    #ook_slicer.digital_binary_slicer = digital.binary_slicer_fb()
    #ook_slicer.complex_to_mag_squared = blocks.complex_to_mag_squared(1)
    #ook_slicer.blocks_uchar_to_float_0 = blocks.uchar_to_float()
    #ook_slicer.add_const = blocks.add_const_vff((-1*threshold, ))
 
            ##################################################
            # Connections
            ##################################################
    ook_slicer.connect((ook_slicer.file_source, 0), (ook_slicer.freq_xlating_fir_filter, 0))
    ook_slicer.connect((ook_slicer.freq_xlating_fir_filter, 0), (ook_slicer.file_sink_0, 0))
 
    # This thing now just pumps it out!
    ook_slicer.start() # Start the flow graph
    ook_slicer.wait() # wait for it to finish and return
 
    ## the returned object is mutated, and some basic info changes
    burst_object.raw_input_file = filename_filtered
    burst_object.samp_rate = burst_object.working_samp_rate
    return burst_object
 
def slice_ook(burst_object):
    ## do an auto threshold function here -- mutate burst_object.threshold
    ### OOK_Demod  script
    ook_slicer = gr.top_block() # Define the container
 
            ##################################################
            # Variables
            ##################################################
    ook_slicer.threshold = threshold = burst_object.threshold
    ook_slicer.filename_raw_rx = filename_raw_rx = burst_object.raw_input_file 
    ook_slicer.filename_sliced = filename_sliced = filename_raw_rx+'.sliced'
 
            ##################################################
            # Blocks
            ##################################################
 
    ook_slicer.file_source = blocks.file_source(gr.sizeof_gr_complex*1, filename_raw_rx, False)
    ook_slicer.complex_to_mag_squared = blocks.complex_to_mag_squared(1)
    ook_slicer.add_const = blocks.add_const_vff((-1*threshold, ))
    ook_slicer.digital_binary_slicer = digital.binary_slicer_fb()
    ook_slicer.blocks_uchar_to_float_0 = blocks.uchar_to_float()
    ook_slicer.file_sink_0 = blocks.file_sink(gr.sizeof_float*1, filename_sliced, False)
    ook_slicer.file_sink_0.set_unbuffered(False)
 
            ##################################################
            # Connections
            ##################################################
    ook_slicer.connect((ook_slicer.file_source, 0), (ook_slicer.complex_to_mag_squared, 0))
    ook_slicer.connect((ook_slicer.complex_to_mag_squared, 0), (ook_slicer.add_const, 0))
    ook_slicer.connect((ook_slicer.add_const, 0), (ook_slicer.digital_binary_slicer, 0))
    ook_slicer.connect((ook_slicer.digital_binary_slicer, 0), (ook_slicer.blocks_uchar_to_float_0, 0))
    ook_slicer.connect((ook_slicer.blocks_uchar_to_float_0, 0), (ook_slicer.file_sink_0, 0))
 
    # This thing now just pumps it out!
 
    ook_slicer.start() # Start the flow graph
    ook_slicer.wait() # wait for it to finish and return
    burst_object.demodulated_output = filename_sliced
    dat_sliced = np.fromfile(burst_object.demodulated_output, dtype="float32")
    sliced_bool = np.array(dat_sliced,np.bool)
    burst_object.slices = Slices(burst=sliced_bool,samp_rate=burst_object.samp_rate)
    return burst_object



def shift_dat(dat,shift,samp_rate):
    ### Shift it
    Fs=samp_rate
    shifting = shift # for pos this shifts "up" or moves everything "left" to center on a higher freq
    # To mix the data down, generate a complex exponential 
    # with phase -f_shift/Fs
    fc = np.exp(-1.0j*2.0*np.pi* shifting/Fs*np.arange(len(dat)))
    shifted = dat * fc
    return shifted
 
    
def get_transitions(sliced_data):
    #create boolean values from the floats -- zero is false, everything else is True
    numpy_bools = np.array(sliced_data,np.bool)
    
    # first get value of where transitions occure
    diff_pos = np.diff(numpy_bools)
    
    # these values are a list of positions of transitions
    diff_pos_loc = np.where(diff_pos)[0]
    
    return diff_pos_loc
 
def get_bad_bursts(transitions, bad_burst_width=20):
    # roll it left so I can ieterate over both and compare current to next values
    next_transitions = np.roll(transitions,-1)
    
    # this stores the location of 'illegitimate' bursts, one-offs
    abberent_bursts = []
    
    # now we look for the illegitimate bursts and get a list of them:
    counter = 0
    for current_loc, next_loc in zip(transitions, next_transitions):
        if next_loc < current_loc:
            burst_slice_range = [next_loc, current_loc]
            continue
        if next_loc - current_loc < bad_burst_width:
            abberent_bursts.append(counter)
        counter = counter + 1
    return abberent_bursts
 
def get_burst_widths(sliced_data,transistion_list):
    true_width_list = []
    false_width_list = []
    counter = 0
    numpy_bools = np.array(sliced_data,np.bool)
    transitions_next = np.roll(transistion_list,-1)
    
    for current_loc, next_loc in zip(transistion_list, transitions_next):
        if next_loc < current_loc:
            burst_slice_range = [next_loc, current_loc]
            continue
        width = next_loc-(current_loc+1)
        value = sum(numpy_bools[current_loc+1:next_loc])
        if value > 10:
            true_width_list.append(width)
        else:
            false_width_list.append(width)
    return [true_width_list,false_width_list,burst_slice_range]
 
def get_discriminators(my_slices):
    
    my_slices.transitions = get_clean_trans(my_slices.burst, my_slices.min_width)
    
    # [true_width_list,false_width_list,burst_slice_range]
    burst_info = get_burst_widths(my_slices.burst,my_slices.transitions)
    # print('burst widths, info:'+str(burst_info))
    
    my_slices.burst_ends = burst_info[2]
    #my_slices.burst = my_slices.burst[(burst_ends[0]+1):(burst_ends[1]+1)]
    
    sorted_true = sorted(burst_info[0])
    true_len = len(sorted_true)
    sorted_false = sorted(burst_info[1])
    false_len = len(sorted_false)
    
    my_slices.ones_threshold = sorted_true[0]+(sorted_true[true_len-1] - sorted_true[0])/2
    my_slices.zeros_threshold = sorted_false[0]+(sorted_false[false_len-1] - sorted_false[0])/2
    
    return my_slices
 
def get_clean_trans(sliced_data, min_width = 20):
    # get initial transistion locations
    i_transitions = get_transitions(sliced_data)
    # print('initial transition locations:'+str(i_transitions))
    # get list of spurious, bad, 'illegitimate' bursts
    bad_bursts = get_bad_bursts(transitions=i_transitions, bad_burst_width=min_width)
    # print('bad bursts:'+str(bad_bursts))
    # this deletes abberations, making the authoritative list of transitions
    transitions = np.delete(i_transitions,bad_bursts)
    # print('final transition locations:'+str(transitions))
    return transitions
 
def get_symbols(my_slices):
    high = False
    last = 0
    high_length_long = []
    high_length_short = []
    low_length_long = []
    low_length_short = []
    
    for val in my_slices.transitions:
        pulse = 0
        pulse = val - last
        last = val
        if high:
            if pulse > my_slices.ones_threshold:
                my_slices.symbols.append(1)
                my_slices.symbols.append(1)
                high_length_long.append(pulse)
            else:
                my_slices.symbols.append(1)
                high_length_short.append(pulse)
            high = False
            continue
        if not high:
            if pulse > my_slices.zeros_threshold:
                my_slices.symbols.append(0)
                my_slices.symbols.append(0)
                low_length_long.append(pulse)
            else:
                my_slices.symbols.append(0)
                low_length_short.append(pulse)
            high = True
    # These convert the samples into microseconds --- based on workinf samp_rate
 
    my_slices.ones_long = np.average(high_length_long)*(1000000/my_slices.samp_rate)
    my_slices.ones_short = np.average(high_length_short)*(1000000/my_slices.samp_rate)
    my_slices.zeros_long = np.average(low_length_long)*(1000000/my_slices.samp_rate)
    my_slices.zeros_short = np.average(low_length_short)*(1000000/my_slices.samp_rate)
    
    return my_slices
    
def get_codes(filename):
    # opens the json, looks at the burst base64 codes #
    # returns a list of all codes in the capture file
    import json
    import decimal  #this gets floats back out from json
 
    rec_file = filename
    #rec_file = 'lgdo2_CF_315.018M_400k1540760759961675.json'
    with open(rec_file) as infile:
        data = json.load(infile,parse_float=decimal.Decimal)
    rec_data = json.loads(data)
 
    base64_list = []
    for thing in range(rec_data['sample_bursts_detected']):
        thing = thing + 1
        ## symbol list is all the zeros and ones, in bitarrays
        symbol_list.append(bitarray(rec_data['bursts'][str(thing)]['symbols']))
        ## the base_64 list converts the bitarrays to base64 encoding
        # those look like: b'qqqoASabbaabaTbbbaTSaaabSbaSbSaabaababA=\n'
        base64_list.append(binascii.b2a_base64(symbol_list[thing-1].tobytes()))
    
    return base64_list