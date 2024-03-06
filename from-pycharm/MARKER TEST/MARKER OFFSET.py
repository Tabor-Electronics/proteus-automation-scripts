# Import required Libraries
import os
import sys
import tempfile
import webbrowser

srcpath = os.path.realpath('..//..//SourceFiles')
sys.path.append(srcpath)
from teproteus import TEProteusAdmin as TepAdmin
from teproteus import TEProteusInst as TepInst
from teproteus_functions_v3 import connect
from teproteus_functions_v3 import disconnect
from teproteus_functions_v3 import set_lib_dir_path
from teproteus_functions_v3 import get_cpatured_header
from teproteus_functions_v3 import gauss_env
from teproteus_functions_v3 import iq_kernel
from teproteus_functions_v3 import pack_kernel_data
from teproteus import TEProteusAdmin, TEProteusInst
from tevisainst import TEVisaInst
# matplotlib notebook
import numpy as np
import time
import ipywidgets as widgets
from IPython.core.debugger import set_trace
from scipy.signal import chirp, sweep_poly
import matplotlib.pyplot as plt

plt.style.use('ggplot')
from scipy import signal
import math
import pdb
import pyvisa as visa
from pyvisa.errors import Error
from pyvisa.resources import resource
# from pyvisa.resources import read
###########################################################################################
class marker_offset_test:
    def __init__(self, ip):
        self.proteus_addr = 'TCPIP::' + ip + '::5025::SOCKET'
        print(ip,self.proteus_addr)


    def execute(self):
        inst = TEVisaInst(self.proteus_addr)
        resp = inst.send_scpi_query('*IDN?')
        print('Connected to: ' + resp)
        inst.default_paranoia_level = 2
        # Reset the instrument
        inst.send_scpi_cmd('*CLS; *RST')
        inst.send_scpi_cmd(':TRACe:DELete:ALL')
        print('Restarting '+resp)
        # Get the model name
        model = inst.send_scpi_query(":SYST:iNF:MODel?")
        print("Model: " + model)
        # Get the DAC mode (8 bits or 16 bits)
        resp = inst.send_scpi_query(':SYST:INF:DAC?')
        if resp == 'M0':
            dac_mode = 16
            bpp = 2
            max_dac = 65535
            wpt_type = np.uint16
            offset_factor = 2
            channels_per_dac = 2
        else:
            dac_mode = 8
            bpp = 1
            max_dac = 255
            wpt_type = np.uint8
            offset_factor = 1
            channels_per_dac = 1
        half_dac = max_dac / 2.0
        print('DAC {0} bits'.format(dac_mode))
        resp = inst.send_scpi_query(':SYST:ERR?')
        if not resp.startswith('0'):
            print("ERROR1", resp)
        # Get number of channels
        resp = inst.send_scpi_query(":INST:CHAN? MAX")
        print("Number of channels: " + resp)
        num_channels = int(resp)
        # Get the maximal number of segments
        resp = inst.send_scpi_query(":TRACe:SELect:SEGMent? MAX")
        max_seg_number = int(resp)
        print("Max segment number: {}".format(max_seg_number))

        # Get the available memory in bytes of wavform-data (per DDR):
        resp = inst.send_scpi_query(":TRACe:FREE?")
        arbmem_capacity = int(resp)
        print("Available memory per DDR: {0:,} wave-bytes".format(arbmem_capacity))
########################################################################################################################
        # Build wave-data and markers-data
        print('Build wave-data and markers-data ..')
        test_success = True
        seg_wave_points = 4096
        ncycles = 1
        cyclelen = seg_wave_points / ncycles

        wave = [None]
        mark = [None]
        # Buffer size for the waveform
        if dac_mode == 16:
            seg_wave_bytes = seg_wave_points * 2  # each point is 2 bytes in size
        else:
            seg_wave_bytes = seg_wave_points
        # seg_wave_bytes = seg_wave_points * 2 # each waveform point is represented as uint16 (2 bytes in size)
        sampling_rate = 2500e6
        # Buffer size for Marker
        if dac_mode == 16 and sampling_rate <= 2500000000:
            marker_cycle = 2
            seg_mark_bytes = seg_wave_points // 4
        else:
            marker_cycle = 8
            seg_mark_bytes = seg_wave_points // 8
        x = np.linspace(start=0, stop=seg_wave_points, num=seg_wave_points, endpoint=False)
        yw = np.fmod(x, cyclelen)
        yw = (yw <= cyclelen / 2) * max_dac
        yw = np.round(yw)
        yw = np.clip(yw, 0, max_dac)
        yw = yw.astype(wpt_type)
        yw.reshape(-1)  # = yw.astype(data_type)
        # Build marker
        x = np.linspace(
            start=0, stop=seg_mark_bytes, num=seg_mark_bytes, endpoint=False)
        y = np.fmod(x, seg_mark_bytes)
        y = (y <= seg_mark_bytes / 2) * 0x33
        y = np.round(y)
        y = np.clip(y, 0, 255)
        # y= np.ones(seg_mark_bytes,np.uint8) *0xFF
        mark = y.astype(np.uint8).reshape(-1)
        del x, y
        segnum = 1
        wav = yw
        mrk = mark

########################################################################################################################
        scope_addr = 'USB0::0x2A8D::0x900E::MY55490134::INSTR' # connect to scope via USB
        try:
            resourceManager = visa.ResourceManager()  # Create a connection (session) to the instrument
            #scope = resourceManager.get_instrument(scope_addr2)
            #scope.write('*CLS;:DISPlay:CGRade:LEVels ')
            scope = resourceManager.open_resource(scope_addr)
            print(scope)
            ## scope acquisition
            # Send *IDN? and read the response
            scope.write('*RST')
            scope.write('*IDN?')
            idn = scope.read()
            print('*IDN? returned: %s' % idn.rstrip('\n'))
        except Error as ex2:
            print('Couldn\'t connect to \'%s\', exiting now...' % scope_addr)
            sys.exit()
        scope.write('AUTOscale')
        time.sleep(2)
        scope.write('*OPC')
        scope.write(':MEASure:CLEar')
        scope.write('*CLS;:DISPlay:CGRade:LEVels ')

        ########################################################################################################################

        sampling_freq = 2500e6
        inst.send_scpi_cmd(':FREQ {}'.format(sampling_freq))
        # Build waveforms
        seglen = 1024
        ncycles = 1
        cyclelen = seglen / ncycles
        waves = [None for _ in range(num_channels)]
        marks = [None for _ in range(num_channels)]
        if (dac_mode == 16) and int(sampling_freq) <= 2.5e9:
            seg_mark_bytes = seglen // 4
        else:
            seg_mark_bytes = seglen // 8
            if dac_mode != 16:
                num_channels = int(num_channels // 2)

        for ii in range(num_channels):
            x = np.linspace(start=0, stop=seglen, num=seglen, endpoint=False)
            y = np.fmod(x, cyclelen)
            y = (y <= cyclelen / 2) * max_dac
            y = np.round(y)
            y = np.clip(y, 0, max_dac)
            if dac_mode == 16:
                waves[ii] = y.astype(np.uint16)
            else:
                waves[ii] = y.astype(np.uint8)
            del x, y
        for ii in range(num_channels):
            cycle_bytes = seg_mark_bytes / ncycles
            x = np.linspace(
                start=0, stop=seg_mark_bytes, num=seg_mark_bytes, endpoint=False)
            ym = np.fmod(x, cycle_bytes)
            ym = (ym <= cycle_bytes / 2) * 255
            ym = np.round(ym)
            ym = np.clip(ym, 0, 255)
            marks[ii] = ym.astype(np.uint8)
            del x, ym

        # num_channels = 4
        for ii in range(num_channels):
            ichan = ii

            channb = ichan + 1
            segnum = ichan % 2 + 1
            wav = waves[ichan]
            mrk = marks[ichan]
            # print('Download wave to segment {0} of channel {1}'.format(segnum, channb))

            # Select channel
            inst.send_scpi_cmd(':INST:CHAN {0}'.format(channb))
            # Define segment
            inst.send_scpi_cmd(':TRAC:DEF {0}, {1}'.format(segnum, seglen))
            # Select the segment
            inst.send_scpi_cmd(':TRAC:SEL {0}'.format(segnum))
            # Increase the timeout before writing binary-data:
            inst.timeout = 30000
            # Send the binary-data:
            inst.write_binary_data(':TRAC:DATA', wav)
            resp = inst.send_scpi_query(':SYST:ERR?')
            resp = resp.rstrip()

            if not resp.startswith('0'):
                print('ERROR: "{0}" after writing binary values'.format(resp))
            # print('Download markers to segment {0} of channel {1}'.format(segnum, channb))
            # Increase the timeout before writing binary-data:
            inst.timeout = 10000

            # Send the binary-data with *OPC? added to the beginning of its prefix.
            inst.write_binary_data(':MARK:DATA', mrk)
            resp = inst.send_scpi_query(':SYST:ERR?')
            resp = resp.rstrip()
            if not resp.startswith('0'):
                print('ERROR: "{0}" after writing binary values'.format(resp))

            # Play the specified segment at the selected channel:
            inst.send_scpi_cmd(':SOUR:FUNC:MODE:SEGM {0}'.format(segnum))
            # Turn on the output of the selected channel:
            inst.send_scpi_cmd(':OUTP ON')
            inst.send_scpi_cmd(':MARK:SEL 1;:MARK:STAT ON')
            inst.send_scpi_cmd(':MARK:SEL 2;:MARK:STAT ON')
            # Turn on the markers of the selected channel
            resp = inst.send_scpi_query(':MARK:SEL? MAX')
            #     print(f'Maximun markers = {resp}')
            resp = resp.rstrip()
            markers_per_chan = int(resp)
            #     print(f'Markers_per_chan = {markers_per_chan}')
            for in_marker in range(markers_per_chan):
                test_success = True
                marker_number = in_marker + 1
                inst.send_scpi_cmd(':MARK:SEL {}; :MARK:STAT ON'.format(marker_number))

                print('\nConnect Channel {1} and its Marker {0} to channel 1 and 2 of the oscilloscope'
                      '\nPress enter '.format(marker_number, channb))
                input()
                marker_offset = [0, 0.2, 0.5]
                off_error = [27, 24, 33]
                mrkr_on_scope = 2
                horizontal_scale_per_divison = 0.1e-6
                vertical_scale_per_divison = 200e-3

                for in_marker_volt in range(len(marker_offset)):
                    inst.send_scpi_cmd(':MARK:VOLT:PTOP 0.5;:MARK:VOLT:OFFS {}'.format(marker_offset[in_marker_volt]))
                    scope.write('*RST;:CHAN1:DISP OFF')
                    time.sleep(1)
                    scope.write(':CHAN{0}:DISP ON;:TIMebase:SCALe {1}'.format(mrkr_on_scope, horizontal_scale_per_divison))
                    scope.write(':CHAN{0}:SCAL {1};:CHAN{0}:INP DC50'.format(mrkr_on_scope, vertical_scale_per_divison))
                    scope.write(':CHANnel{0}:OFFSet 0'.format(channb, marker_offset[in_marker_volt]))
                    time.sleep(1)
                    scope.write(':MEASure:VAVerage DISPlay,CHANnel{}'.format(mrkr_on_scope))
                    time.sleep(5)
                    scope.write(':MEASure:RESults?')
                    result = scope.read()
                    result_offset = float(result.split(',')[2])
                    # print(result_offset)
                    diff = abs(result_offset - marker_offset[in_marker_volt])
                    if diff >= off_error[in_marker_volt] * 1e-3:
                        test_success = False
                        print('Test FAIL for offset {}V'.format(marker_offset[in_marker_volt]))
                    else:
                        print('Test pass for offset {}V'.format(marker_offset[in_marker_volt]))
                if (test_success):
                    print('Voltage Offset test pass for channel {0} marker {1} '.format(channb, marker_number))
                else:
                    print('Voltage Offset test Fail for channel {0} marker {1} '.format(channb, marker_number))
        if (test_success):
            print('Test pass for marker voltage offset')
        else:
            print('Test failed for marker voltage offset')
        disconnect()

        print("Test completed")

test = marker_offset_test('192.90.70.22')
test.execute()