#!/usr/bin/env python
# Purpose: This script loads a raw iq samples recorded by a 
#          Software Defined Radio.
# Script command line args: sample_rate, offset, center_freq, 
#                           filename
# Output: The symbols represented by the signal. 
#         A symbol is a series detectable data points that will 
#         map out to a bit (1 or 0)
 
# helper functions that I keep in basic_help.py
import basic_help as bh
 
# standard library for doing data science
import numpy as np 
 
# to parse command line options
from optparse import OptionParser
 
# for good logs and output
import time
import json
import decimal  #this gets floats back out from json
 
# doing some deconding of symbols
from bitarray import bitarray 
import binascii

## Setup the option parser
parser = OptionParser()
parser.add_option("-i", "--infile", dest="infile",
                  help="what raw iq data file to use?", metavar="INFILE")
parser.add_option("-o", "--outfile", dest="outfile",
                  help="what data file to output?", metavar="OUTFILE")       
parser.add_option("-s", "--samp_rate", dest="samp_rate", type="float", default=8e6,
                  help="what is the input file sample rate?", metavar="SAMP_RATE")    
parser.add_option("-t", "--threshold", dest="threshold", type="float",
                  help="what is the threshold for slicing?", metavar="THRESHOLD")                  
parser.add_option("-e", "--offset", dest="offset", type="float", default=0.0,
                  help="what is the offset from center?", metavar="OFFSET")      
parser.add_option("-r", "--transition", dest="transition", type="float", default=10e3,
                  help="what is the transition width for lowpass filter?", metavar="TRANSITION")      
parser.add_option("-w", "--working_samp", dest="working_samp", type="float", default=400e3,
                  help="what is the working sample rate?", metavar="WORKING_SAMPLE_RATE")      
parser.add_option("-c", "--cutoff", dest="cutoff", type="float", default=100e3,
                  help="what is the filter cuttoff, stop?", metavar="CUTOFF")   
parser.add_option("-f", "--center_freq", dest="center_freq", type="float", default=315.018e6,
                  help="what is the center freq?", metavar="CENTER_FREQ")    
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()
# this wil be a script that takes several arguments and slices a raw ook modulated sample


# We are inside some application


if not options.infile:
    print("\t You must specify an infile.  \n\t Sample rate defaults to 8.0Ms \n\t threshold defaults to auto-generated")
    exit(-1)

#print(options)

my_burst = bh.Burst()
my_burst.name = options.infile
#raw_file = './gdo_2M_200k_315.218M_400k'
my_burst.center_freq = options.center_freq
my_burst.samp_rate = options.samp_rate
my_burst.offset = options.offset
my_burst.threshold = options.threshold
my_burst.raw_input_file = options.infile



# in seconds, for the gdo it is .03
burst_seperation_factor = .03
# in seconds, used for gdo 
preamble_gap_factor = .003
# minimum symbol width, in seconds -- about half what is normal
#for gdo, period is .0003000
min_symbol_factor = .00002
#min_symbol_factor = .00008125

logfile = './raw_to_symbols.log'
l = bh.AsyncWrite(logfile)
start_t = time.localtime()
start_time = time.time()
script_start = time.asctime(start_t)
gathered_log =                ('\n********************** Log START ***********************\n')
gathered_log = gathered_log + ('\n\tStart  Time: '+ str(script_start))
gathered_log = gathered_log + ('\n\tInput  File: '+ str(my_burst.raw_input_file))
gathered_log = gathered_log + ('\n\tSample Rate: {0}'.format(my_burst.samp_rate))
gathered_log = gathered_log + ('\n\tCenter Freq: {0}'.format(my_burst.center_freq))
gathered_log = gathered_log + ('\n\tBurst Seperation Detection width (s): {0:1.4f}'.format(burst_seperation_factor))
gathered_log = gathered_log + ('\n\tPreamble Gap detection width     (s): {0:1.6f}'.format(preamble_gap_factor))
gathered_log = gathered_log + ('\n\tMinimal Symbol detection width   (s): {0:1.9f}'.format(min_symbol_factor))
gathered_log = gathered_log + ('\n\tRunning ****\n')
l.log(gathered_log)
gathered_log = ''

'''
Showing relationship between samples (no unit) , sample rates (samples per second) , and duration of pulses (seconds)
samp_rate * duration = samples
samples / samp_rate = duration

Example:
samp_rate = 400e3
samples = 200
duration = .0005
400e3 * .0005 = 200
200 / 400e3 = .0005

'''
sample_uid = "{0:10d}".format(int(start_time*1000000))
sample_dict = {'name':my_burst.raw_input_file,
               'uid':sample_uid,
               'date':script_start,
               'samp_rate':my_burst.samp_rate,
               'center_freq':my_burst.center_freq,
               'burst_seperation_factor':burst_seperation_factor,
               'preamble_gap_factor':preamble_gap_factor,
               'min_symbol_factor':min_symbol_factor
              }
if not options.outfile:
    outfile = my_burst.raw_input_file+'.json'
else: 
    outfile = options.output

dist_between_bursts = burst_seperation_factor*my_burst.samp_rate

#load the data into numpy array (and into Burst.data_arr)
my_burst.load_data()
print(my_burst.data_arr)
# the complex data is converted into magnitudes
mags = np.absolute(my_burst.data_arr)

# the threshold is automatically calculated by finding the maximum magnitude, then dividing it by two
# this assumes no "noise" was included directly in this filtered data
max_val = np.amax(mags)
mythresh = max_val/2
gathered_log = gathered_log + ('\n\tWhole Sample max amplitude: '+ str(max_val))
gathered_log = gathered_log + ('\n\tWhole Sample threshold: '+ str(mythresh))
sample_dict['sample_max_amp'] = max_val.tolist()
sample_dict['sample_initial_threshold'] = mythresh.tolist()
if options.threshold:
    mythresh = options.threshold
    sample_dict['sample_initial_threshold'] = mythresh
    gathered_log = gathered_log + ('\n\tUsing threshold: '+ str(mythresh))

# numpy creates an array of positions where the magnitudes are higher than the threshold value
# this accomplishes nearly the same thing as gnuradio binary slicer, a kind of metadata though
sliced_data = np.array(np.where(mags > mythresh )[0])
            
# the differences between each consecutive number is calculated and placed in this new array
new = np.array(np.ediff1d(sliced_data))
# gathered_log = gathered_log + ('\n \t initial sliced positions of past threshold values (new) \n') + str(new.size)

# if the differences between the sliced data postions are greater than the average distance between the normal bursts,
# then the number indicates a new burst, and this array holds the starting location of each noew burst.
burst_segments = np.where(new > dist_between_bursts )[0]
# print burst_segments
# chunk through and slice out each burst to analyze seperately
cur = 0
the_bursts = []
for nex in burst_segments:
    the_bursts.append(np.absolute(my_burst.data_arr[sliced_data[cur]:sliced_data[nex]]))
    cur = nex+1

# detect the gap between the preamble and data, using same method as above, but on single burst
# based on specgram, looks like near .001 seconds, anything smaller should be treated as gap between data symbols
gap_after_preamble = preamble_gap_factor*my_burst.samp_rate

gathered_log = gathered_log + ('\n\tBursts Detected: ' + str(len(the_bursts)))
l.log(gathered_log)
sample_dict['sample_bursts_detected'] = len(the_bursts)    
gathered_log = ''

burst_number = 0
sample_dict['bursts'] = {}
# now generate symbols for all bursts in this sample
for burst in the_bursts:
    
    burst_number = burst_number + 1
    sample_dict['bursts'][str(burst_number)] = {'number':burst_number}
    # the complex data is converted into magnitudes
    mags = np.absolute(burst)

    # numpy creates an array of positions where the magnitudes are higher than the threshold value
    # this accomplishes nearly the same thing as gnuradio binary slicer, a kind of metadata though
    sliced_data = np.array(np.where(mags > mythresh )[0])

    # the differences between each consecutive number is calculated and placed in this new array
    new = np.array(np.ediff1d(sliced_data))
    gathered_log = gathered_log + ('\n\n\t******************* Burst '+str(burst_number)+' Data ***********************')
### Work on detecting the preamble --
    # if the differences between the sliced data postions are greater than the average distance between the normal bursts,
    # then the number indicates a new burst, and this array holds the starting location of each noew burst.
    preamble_loc = np.where(new > gap_after_preamble )[0]
    #print("preamble size and thing",preamble_loc.size,preamble_loc)
    preamble = False
    ## basically just stopped checking for premables -- this needs further work for non preamble bursts
    if preamble_loc.size < 0 and preamble_loc[0] is not 0:
        #print("preamble loc & size",preamble_loc.size,preamble_loc)
        gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' preamble gap (s): '+ str( new[preamble_loc][0] / my_burst.samp_rate ))
        sample_dict['bursts'][str(burst_number)]['preamble_gap'] = (new[preamble_loc][0] / my_burst.samp_rate)
        #preamble = True
        if preamble_loc[0] is 0: preamble = False
    else:
        sample_dict['bursts'][str(burst_number)]['preamble_gap'] = 'none detected'
        preamble = False

    ## now to find widths of minimum symbols --- need to have an idea of minimum sumbol length in samples
    # zooming in on plot of preamble might help -- choose something like half the preamble bits
    # they look about 75 samples wide so...for samp_rate of 4e5, that is 187 microseconds (.0001875 seconds)
    # half that is 
    min_symbol_length = int(min_symbol_factor * my_burst.samp_rate)
    symbol_vals = []
    # this slice is the preamble: (sliced_data[0],sliced_data[preamble_loc[0]])
    # and this is the base  (sliced_data[preamble_loc[0]+1],the_bursts[0].size)
    if preamble:
        print("preamble loc 0",preamble_loc[0])
        the_burst_pramble = np.absolute(burst[sliced_data[0]:sliced_data[preamble_loc[0]]])
        preamble_threshold  = np.amax(the_burst_pramble)/2
        gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Recalculated preamble threshold (s): '+ str(preamble_threshold))
        sample_dict['bursts'][str(burst_number)]['preamble_threshold'] = preamble_threshold.tolist()
        preamble_bit_locs = np.array(np.where(the_burst_pramble > preamble_threshold )[0])  
        preamble_gap_locs = np.array(np.where(the_burst_pramble < preamble_threshold )[0])  
        preamb_thresh_bit_diffs = np.diff(preamble_bit_locs)
        preamb_thresh_gap_diffs = np.diff(preamble_gap_locs)
        preamb_bit_locations = np.where(preamb_thresh_bit_diffs > min_symbol_length)
        preamb_gap_locations = np.where(preamb_thresh_gap_diffs > min_symbol_length)
        preamb_bit_lengths = np.diff(preamb_bit_locations)
        preamb_gap_lengths = np.diff(preamb_gap_locations)
        # now combine the lists of locations of gaps and bits, and put them in order
        symbol_locs = np.array([])
        bit_list = []
        gap_list = []
        for locs in preamb_bit_locations:
            val = preamble_bit_locs[locs]
            symbol_locs = np.append(symbol_locs,val)
            bit_list.append(val)
        
        for locs in preamb_gap_locations:
            val = preamble_gap_locs[locs]
            symbol_locs = np.append(symbol_locs,val)
            gap_list.append(val)


        symbol_list = np.sort(symbol_locs)
        # presumt the symbols start with a 1?
        
        # iterate through the_burst_base, testing the lengths/values of the symbols and assigning symbols
        for pos in symbol_list:
            test = np.in1d(bit_list,pos)    #test = np.isin(bit_list,pos)
            if len(np.where(test == True)[0]):
                symbol_vals.append(1)
            test = np.in1d(gap_list,pos)    #test = np.isin(gap_list,pos)
            if len(np.where(test == True)[0]):
                symbol_vals.append(0)           
        if len(symbol_vals) > 0:
            sample_dict['bursts'][str(burst_number)]['preamble'] = symbol_vals
            gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' preamble detected: '+ str(symbol_vals))
            # putting gap in between preamble and rest of data
            symbol_vals = symbol_vals + [0,0,0,0,0,0,0,0,0,0]
        else:
            gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' preamble not detected :(')
            sample_dict['bursts'][str(burst_number)]['preamble'] = 'none detected'


    if preamble:
        the_burst_base = np.absolute(burst[sliced_data[preamble_loc[0]+1]:]) #burst.size
    else:
        the_burst_base = np.absolute(burst[:]) #burst.size

    # this takes half the maximum value found in the magnitueds     
    base_threshold  = np.amax(the_burst_base)/2
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Recalculated base data threshold(s): '+ str(base_threshold))
    sample_dict['bursts'][str(burst_number)]['base_threshold'] = base_threshold.tolist()

    #these are the locations where the magnitudes reach the threshold
    base_bit_locs = np.array(np.where(the_burst_base > base_threshold )[0])

    # opposite list for gaps
    base_gap_locs = np.array(np.where(the_burst_base < base_threshold )[0]) 

    # now where the difference between locations is only one, that is part of a consecutive bit, where it is more than one, that is a gap
    # this returns an array of th edifferences in consecutive values
    base_thresh_bit_diffs = np.diff(base_bit_locs)
    # same for gaps
    base_thresh_gap_diffs = np.diff(base_gap_locs)

    #   base_bit_locs[base_thresh_gap_diffs[base_bit_locations]] == location in the burst base  
    # this will return an array of where those values are more than our min_symbol length
    base_bit_locations = np.where(base_thresh_bit_diffs > min_symbol_length)[0]
    #now get the gaps
    base_gap_locations = np.where(base_thresh_gap_diffs > min_symbol_length)[0]
    # remember that the first location is the beginning of a bit

    ## correcting info lost from the diffs, adding the next bit/gap locations
    #  appending a value that is one more than the previous last value:
    ### worked with GDO -- which hada preamble, not with SolidRemote -- need to analyze its preamble
    #preamb_bit_locations = np.append(preamb_bit_locations,(preamb_bit_locations[0][len(preamb_bit_locations[0])-1])+1)
    #preamb_gap_locations = np.append(preamb_gap_locations,(preamb_gap_locations[0][len(preamb_gap_locations[0])-1])+1)

    ## correcting info lost from the diffs, adding the next bit/gap locations
    #  appending a value that is one more than the previous last value:
    if base_bit_locations.size > 0:
        base_bit_locations = np.append(base_bit_locations,(base_bit_locations[-1])+1)
    if base_gap_locations.size > 0:
        base_gap_locations = np.append(base_gap_locations,(base_gap_locations[-1])+1)


    #and this shows us the "lengths" between bits
    base_bit_lengths = np.diff(base_bit_locations)
    base_gap_lengths = np.diff(base_gap_locations)       

    # all the bitwidths look very normal -- 
    # can now populate threshold values for bits - to determine when assigning 11 vs 1
    ones_threshold = np.average(base_bit_lengths)
    zeros_threshold = np.average(base_gap_lengths)   

    symbol_locs = np.array([])
    bit_list = []
    gap_list = []

    for locs in base_bit_locations:
        val = base_bit_locs[locs]
        symbol_locs = np.append(symbol_locs,val)
        bit_list.append(val)
        
    for locs in base_gap_locations:
        val = base_gap_locs[locs]
        symbol_locs = np.append(symbol_locs,val)
        gap_list.append(val)

    symbol_list = np.sort(symbol_locs)

    ## create containers to capture data on lengths of long and short bits and gaps
    bit_long = []
    bit_short = []
    gap_long = []
    gap_short = []

    prev = 0
    # iterate through the_burst_base, testing the lengths/values of the symbols and assigning symbols
    for pos in symbol_list:
        dat = pos-prev
        test = np.in1d(bit_list,pos)    # test = np.isin(bit_list,pos)
        if len(np.where(test == True)[0]):
            if dat > ones_threshold:
                bit_long.append(dat)
                symbol_vals.append(1)
                symbol_vals.append(1)
            else:
                bit_short.append(dat)
                symbol_vals.append(1)
            prev = pos
            continue
        test = np.in1d(gap_list,pos)    #test = np.isin(gap_list,pos)
        if len(np.where(test == True)[0]):
            if dat > zeros_threshold:
                gap_long.append(dat)
                symbol_vals.append(0)
                symbol_vals.append(0)
            else:
                gap_short.append(dat)
                symbol_vals.append(0)
            prev = pos

    # for good measure....probably lost at least one 1, if not two. need a way to check it
    symbol_vals = symbol_vals + [1]
    # log the measurement data of bit widths
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Long Pulse  avg (s): {0:1.6f}'.format(np.average(bit_long)/ my_burst.samp_rate))
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Short Pulse avg (s): {0:1.6f}'.format(np.average(bit_short)/ my_burst.samp_rate))
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Long Gap    avg (s): {0:1.6f}'.format(np.average(gap_long)/ my_burst.samp_rate))
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Short Gap   avg (s): {0:1.6f}'.format(np.average(gap_short)/ my_burst.samp_rate))
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Symbols Detected   : {0:4d}'.format(len(symbol_vals)))
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' Symbols: '+ str(symbol_vals))
    bits = bitarray(symbol_vals)
    the64 = binascii.b2a_base64(bits.tobytes())
    hexstr = binascii.hexlify(bits.tobytes())
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' base64: '+ str(the64))
    gathered_log = gathered_log + ('\n\t\tBurst:'+str(burst_number)+' hexlify: '+ str(hexstr))
    l.log(gathered_log)
    gathered_log = ''
    ## store this dataset as a .json to be referenced later, perhaps in transmission
    sample_dict['bursts'][str(burst_number)]['bit_long'] = (np.average(bit_long)/ my_burst.samp_rate).tolist()
    sample_dict['bursts'][str(burst_number)]['bit_short'] = (np.average(bit_short)/ my_burst.samp_rate).tolist()
    sample_dict['bursts'][str(burst_number)]['gap_long'] = (np.average(gap_long)/ my_burst.samp_rate).tolist()
    sample_dict['bursts'][str(burst_number)]['gap_short'] = (np.average(gap_short)/ my_burst.samp_rate).tolist()
    sample_dict['bursts'][str(burst_number)]['symbol_nums'] = len(symbol_vals)
    sample_dict['bursts'][str(burst_number)]['symbols'] = symbol_vals
    sample_dict['bursts'][str(burst_number)]['base64'] = str(the64)
    sample_dict['bursts'][str(burst_number)]['hexlify'] = str(hexstr)

end_time = time.time()
sample_dict['run_time'] = end_time-start_time
gathered_log = gathered_log + ('\n\nTotal Run Time (seconds): {0:1.4f}'.format(end_time-start_time))
gathered_log = gathered_log +('\n************************* END **************************\n\n\n')
l.log(gathered_log)
## now store the sample data for future consumption
dumped = json.dumps(sample_dict)
with open(outfile,'w') as my_json_file:
    json.dump(dumped, my_json_file)
exit(0)
#with open(outfile) as infile:
#    data = json.load(infile,parse_float=decimal.Decimal)
#sample_dict = json.loads(data)
#
### From dict to dict example
'''
dumped = json.dumps(sample_dict)
with open('outfile','w') as my_json_file:
    json.dump(dumped, my_json_file)
with open('outfile') as infile:
    data = json.load(infile,parse_float=decimal.Decimal)
loaded = json.loads(data)
'''