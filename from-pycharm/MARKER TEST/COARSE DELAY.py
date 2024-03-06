## Import required Libraries
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
class coarse_delay_test:
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
        scope_addr = 'USB0::0x2A8D::0x900E::MY55490134::INSTR'  # connect to scope via USB
        try:
            resourceManager = visa.ResourceManager()  # Create a connection (session) to the instrument
            # scope = resourceManager.get_instrument(scope_addr2)
            # scope.write('*CLS;:DISPlay:CGRade:LEVels ')
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
        scope.write(':SYSTem:HEADer OFF')
        scope.write('CDIS')


############################################################################################################################

        inst.send_scpi_cmd(':ROSC:SOUR INT')
        sampling_rate = [1.25e9, 2.5e9, 9e9]
        del_range = [28, 32, 128]
        l_range = [16e-9, 9.6e-09, 13.3e-9]
        h_range = [28.8e-9, 16e-9, 15.1e-9]
        # Build waveforms
        seglen = 1024
        ncycles = 1
        cyclelen = seglen / ncycles
        waves = [None for _ in range(num_channels)]
        marks = [None for _ in range(num_channels)]

        resp = inst.send_scpi_query(':MARK:SEL? MAX')
        # print(f'Maximun markers = {resp}')
        resp = resp.rstrip()
        markers_per_chan = int(resp)

        for in_sr in range(len(sampling_rate)):

            if (dac_mode == 16) and int(sampling_rate[in_sr]) <= 2.5e9:
                seg_mark_bytes = seglen // 4
                num_channels = [1, 2, 3, 4]
                num_markers = [1, 2]

            else:
                seg_mark_bytes = seglen // 8
                #         num_channels = int(num_channels // 2)
                num_channels = [1, 3]
                num_markers = [1, 2, 3, 4]

            for ii in num_channels:
                channb = ii
                segnum = 1
                x = np.linspace(start=0, stop=seglen, num=seglen, endpoint=False)
                y = np.fmod(x, cyclelen)
                y = (y <= cyclelen / 2) * max_dac
                y = np.round(y)
                y = np.clip(y, 0, max_dac)
                if dac_mode == 16:
                    waves[0] = y.astype(np.uint16)
                else:
                    waves[0] = y.astype(np.uint8)
                del x, y
                cycle_bytes = seg_mark_bytes / ncycles
                x = np.linspace(start=0, stop=seg_mark_bytes, num=seg_mark_bytes, endpoint=False)
                ym = np.fmod(x, cycle_bytes)
                ym = (ym <= cycle_bytes / 2) * 255
                ym = np.round(ym)
                ym = np.clip(ym, 0, 255)
                marks[0] = ym.astype(np.uint8)

                del x, ym
                wav = waves[0]
                mrk = marks[0]
                inst.send_scpi_cmd(':FREQ {}'.format(sampling_rate[in_sr]))
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
                # Increase the timeout before writing binary-data:
                inst.timeout = 10000
                # Play the specified segment at the selected channel:
                inst.send_scpi_cmd(':SOUR:FUNC:MODE:SEGM {0}'.format(segnum))
                inst.send_scpi_cmd(':VOLT 0.5;:VOLT:OFFS 0')
                # Turn on the output of the selected channel:
                inst.send_scpi_cmd(':OUTP ON')
                # inst.send_scpi_cmd(':INST:CHAN {0}'.format(channb))
                # Send the binary-data with *OPC? added to the beginning of its prefix.
                inst.write_binary_data(':MARK:DATA', mrk)
                resp = inst.send_scpi_query(':SYST:ERR?')
                resp = resp.rstrip()
                if not resp.startswith('0'):
                    print('ERROR: "{0}" after writing binary values'.format(resp))
                print('\t***Test for sampling rate = {} GSas***'.format(int(sampling_rate[in_sr]) / 1e9))

                channb_on_scope = 1
                mrkr_on_scope = 2
                horizontal_scale_per_divison = 0.1e-6
                vertical_scale_per_divison = 200e-3
                for marker_number in num_markers:
                    test_success = True
                    # marker_number = imarker
                    print(
                        'Connect channel {0} and its marker {1} of channel {2} and channel {3} to the scope\nPress enter to continue'.format(
                            channb, marker_number, channb_on_scope, mrkr_on_scope))
                    input()
                    inst.send_scpi_cmd(':INST:CHAN {0}'.format(channb))
                    inst.send_scpi_cmd(':MARK:SEL {};:MARK:STAT ON'.format(marker_number))
                    inst.send_scpi_cmd(':MARK:VOLT:PTOP 1.0;:MARK:VOLT:OFFS 0')
                    inst.send_scpi_cmd(':MARK:DEL:COAR 0')
                    scope.write('*RST;:CHAN1:DISP ON')
                    scope.write(
                        ':CHAN{0}:DISP ON;:TIMebase:SCALe {1}'.format(mrkr_on_scope, horizontal_scale_per_divison))
                    scope.write(':CHAN{0}:SCAL {1};:CHAN{0}:INP DC50'.format(mrkr_on_scope, vertical_scale_per_divison))
                    scope.write(':CHANnel{0}:OFFSet 0')
                    scope.write(':MEASure:DELTAtime:DEF RISing,1,MIDD,RISing,1,MIDDle')
                    time.sleep(1)
                    scope.write(':MEASure:DELTAtime CHANnel{0},CHANnel{1}'.format(channb_on_scope, mrkr_on_scope))
                    time.sleep(5)
                    scope.write(':MEASure:RESults?')
                    result = scope.read()
                    initial_delay = float(result.split(',')[2])
                    inst.send_scpi_cmd(':MARK:SEL {};:MARK:STAT ON'.format(marker_number))
                    inst.send_scpi_cmd(':MARK:VOLT:PTOP 1.0;:MARK:VOLT:OFFS 0')
                    inst.send_scpi_cmd(':MARK:DEL:COAR {}'.format(del_range[in_sr]))
                    time.sleep(1)
                    scope.write(':MEASure:DELTAtime CHANnel{0},CHANnel{1}'.format(channb_on_scope, mrkr_on_scope))
                    time.sleep(5)
                    scope.write(':MEASure:RESults?')
                    result = scope.read()
                    temp_delay = float(result.split(',')[2])
                    diff = temp_delay - initial_delay

                    if diff < l_range[in_sr] or diff > h_range[in_sr]:
                        test_success = False
                        print('Test1 FAIL for coarse delay {} points'.format(del_range[in_sr]))
                    else:
                        print('Test1 pass for coarse delay {} points'.format(del_range[in_sr]))
                    print("Press enter to continue to test 2")
                    input()
                    inst.send_scpi_cmd(':MARK:SEL {};:MARK:STAT ON'.format(marker_number))
                    inst.send_scpi_cmd(':MARK:VOLT:PTOP 1.0;:MARK:VOLT:OFFS 0')
                    inst.send_scpi_cmd(':MARK:DEL:COAR -{}'.format(del_range[in_sr]))
                    scope.write(':MEASure:DELTAtime CHANnel{0},CHANnel{1}'.format(channb_on_scope, mrkr_on_scope))
                    time.sleep(5)
                    scope.write(':MEASure:RESults?')
                    result = scope.read()
                    temp_delay = float(result.split(',')[2])
                    # print(result_offset)
                    diff = temp_delay - initial_delay

                    if diff < -(h_range[in_sr]) or diff > -(l_range[in_sr]):
                        test_success = False
                        print('Test2 FAIL for coarse delay -{} points '.format(del_range[in_sr]))
                    else:
                        print('Test2 pass for coarse delay -{} points'.format(del_range[in_sr]))
                if (test_success):
                    print('Test pass for marker coarse delay at sampling rate {} GSas'.format(int(sampling_rate[in_sr]) / 1e9))
                else:
                    print('Test fail for marker coarse delay at sampling rate {} GSAs'.format(int(sampling_rate[in_sr]) / 1e9))

        disconnect()
        print("Test completed")


test = coarse_delay_test('192.90.70.22')
test.execute()
