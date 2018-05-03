#!/home/dada/frb_detector/mpsr/bin/python

import os
from ctypes import *
#functions = CDLL(os.environ['DADA_ROOT']+'/lib/'+'frb_detector_wrapper.so')
#libc = CDLL("libc.so.6")

import numpy as np
import atexit
from multiprocessing import Queue,Process,Manager,Value,Lock
import threading
import time
import sys
sys.path.append(os.environ['DADA_ROOT']+'/lib/')
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import logging
import socket
import cPickle
import argparse
from helpers import parse_cfg,client_control_monitor,daemonize,sigHandler,saveFBank
from helpers import create_xml_elem,delpid,my_snr
import signal
import datetime

from sklearn.externals import joblib

class CandidateFeatures(Structure):
    _fields_ = [("sn",c_float),("sn_0",c_float),("width",c_int),
            ("nstart",c_int),("nend",c_int),("F1",c_float),
            ("F2",c_float),("F3",c_float),("n_rms_mask",c_int),
            ("sn_rms",c_float),("mod_ind",c_float),
            ("mod_indT",c_float), ("isphonecall",c_int)]


class RFIWriterThread(threading.Thread):
    def __init__(self,bp_numb,rfi_writer_queue,name=None):
        self.rfi_file = None
        self.rfi_writer_queue = rfi_writer_queue
        self.bp_numb = bp_numb
        super(RFIWriterThread,self).__init__(name = name)
    def run(self):
        logging.info("RFI-writer thread initiated")
        while True:
            rfi_str = self.rfi_writer_queue.get()
            if rfi_str is None:
                logging.info("RFI-writer thread exited")
                if self.rfi_file != None:
                    self.rfi_file.close()
                break
            try:
                self.rfi_file.write(rfi_str)
            except ValueError:
                pass
    def change_file_name(self,utc):
        if self.rfi_file == None: #For first obs
            t = time.time()
            while time.time() - t < 20:
                try:
                    f_dir = FIL_FILE_DIR+"/BP"+str(self.bp_numb).zfill(2)+\
                            "/"+utc+"/FB/candidates.list.BP"+str(self.bp_numb).zfill(2)
                    self.rfi_file = open(f_dir,"a+")
                    logging.info("successfully opened '"+f_dir+\
                            "' for rfi logging")
                    return
                except IOError:
                    time.sleep(0.5)
            logging.critical("Couldn't open "+f_dir+" after 20 sec of trying")
        else:
#            self.empty_queue()
#            self.rfi_file.close()
            t = time.time()
            while time.time() - t < 20:
                try:
                    f_dir = FIL_FILE_DIR+"/BP"+str(self.bp_numb).zfill(2)+\
                            "/"+utc+"/FB/candidates.list.BP"+str(self.bp_numb).zfill(2)
                    self.rfi_file = open(f_dir,"a+")
                    logging.info("successfully opened '"+f_dir+\
                            "' for rfi logging")
                    return
                except IOError:
                    time.sleep(0.5)
            logging.critical("Couldn't open "+f_dir+" after 20 sec of trying")
    def terminate_writer(self):
        logging.info("Terminating writer thread")
        time.sleep(0.2) #Give time to flush file
        self.empty_queue() #Empty queue just in case
        self.rfi_writer_queue.put(None) #Send poison pill
    def empty_queue(self):
        while not self.rfi_writer_queue.empty():
            _ = self.rfi_writer_queue.get(timeout=0.1)
    def close_file(self):
        try:
            self.rfi_file.close()
            logging.debug("Candidates log file sucessfully closed")
        except AttributeError:
            if self.rfi_file is None:
                logging.critical("Candidate.list file is a NoneType")
            logging.critical("Candidates.list file couldn't be closed")


def sort_features(ftrs):
    """ Function that takes in CandidateFeatures object, and returns a numpy
    array that serves as an input for the classifier"""
#    sorted_features = np.array([ftrs.width,ftrs.sn/ftrs.sn_0,ftrs.F1,ftrs.F2,
#        ftrs.F3,ftrs.sn_rms,ftrs.n_rms_mask,ftrs.mod_ind,ftrs.mod_indT])
    sorted_features = np.array([ftrs.box,ftrs.F1,ftrs.F2,ftrs.F3,ftrs.event,
        ftrs.left,ftrs.right,ftrs.event_div,ftrs.mean_off,ftrs.std_off,
        ftrs.sig_0,ftrs.sig_1,ftrs.sig_2,ftrs.ks_d,ftrs.ks_p,ftrs.sw_w,
        ftrs.sw_p,ftrs.Mod_ind,ftrs.Mod_indT])
    return sorted_features.reshape(1,-1)


def classify(features,threshold=0.5):
    """ Classifier function

    Args:
        features (np.ndarray): the output of sort_features
        (optional) threshold: the detection probability threshold

    Returns:
        (bool): True if candidate is an FRB
        prob (float): The probability of being an FRB
    """
    y=clf.predict_proba(features)[0]
    ind = np.where(clf.classes_=="PULSES")[0]
    if y[ind]>threshold:
        return True,float(y[ind])
    return False,float(y[ind])

def process_monitor_thread(process_list,refresh_time=10):
    """Function that monitors processes, and reports if any is dead
        
        Args:
            process_list (list): list of process objects to monitor
            refresh_time (float): refresh time in seconds

    """
    logging.debug("Process monitor thread initiated")
    boolian_list = [True for i in range(len(process_list))]
    dead_count = 0
    n_procs = len(process_list)
    while True:
        ii = 0
        for proc in process_list:
            if not proc.is_alive() and boolian_list[ii]:
                logging.critical("%i is dead",proc.pid)
                boolian_list[ii] = False
                ii+=1
                dead_count+=1
        if dead_count == n_procs:
            logging.critical("All processing slaves are dead")
            os.kill(os.getpid(), signal.SIGTERM)
        time.sleep(refresh_time)


def terminate_all(proc_list,in_queue):
    logging.critical("Terminating sub-processes")
    n_proc = len(proc_list)
    for i in range(n_proc):
#        Poison pill approach
        in_queue.put(None)
    time.sleep(0.5)
    counter = 3
    for i in range(3):
        dead = 0
        for proc in proc_list:
            if proc.is_alive():
                dead += 1
        if dead == 0:
            return
        else:
            logging.critical("Some child processes are still alive, sending"\
                    +" sigkill in %s",counter)
            counter -= 1
            time.sleep(1)
    for proc in proc_list:
        if proc.is_alive():
            os.kill(proc.pid,signal.SIGKILL)


def incoherent_beam_proc(inc_beam_port):
    """
    args:
        inc_beam_port (int): Port number to communicate incoherent beam
            statistics
    
    Should receive a pickled dictionary with entries:
        - utc
        - t_sample
        - H_dm
        - H_w
        - source_name
    Add 'EOF' to the end of the serialized dictionary to insure stopping
    broadcast

    Sends back a pickled dictionary of incoherent beam statistics
    """
    logging.debug("Innitialising incoherent beam thread")
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
    host = socket.gethostname().split(".")[0]
    s.bind((host,inc_beam_port))
    s.listen(10)
    logging.debug("Binded to %s %s",host,inc_beam_port)
    while True:
        #Accept connection from any bf node
        conn,addr = s.accept()
        msg = recvall(conn)
        req_info = cPickle.loads(msg)
        logging.debug("Received signal %s from %s",req_info,conn.getpeername())
        search_dir = FIL_FILE_DIR+'/BP00/'+req_info['utc']+\
                '/'+req_info['source_name']+'/BEAM_001/'+\
                req_info['utc']+'.fil'
        t_sample = req_info['t_sample']
        H_dm = req_info['H_dm']
        H_w = req_info['H_w']
        #compute SNR
        logging.debug("Computed SNR, trying to send it back")
        snr = my_snr(t_sample, H_dm, H_w, search_dir)
        response_dict = {'snr':snr}
        #Send the SNR back to respective node
        conn.sendall(cPickle.dumps(response_dict)+"EOF")
        logging.debug("SNR sent")

def get_snr_incoherent(t_sample, H_dm, H_w, utc, search_dir):
    d = {'utc':utc,'t_sample':t_sample,'H_dm':H_dm,'H_w':H_w,
            'source_name':'FB'}
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((INCOHERENT_BEAM_HOST,INCOHERENT_BEAM_PORT))
    s.send(cPickle.dumps(d)+"EOF")
    resp = cPickle.loads(recvall(s))
    return resp['snr']



def pulsar_pulse_classifier(utc,host,port):
    logging.info("pulsar_pulse_classifier => loaded. Listening on (%s, %s)",
            host,port)
    while True:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
        s.bind((host,port))
        s.listen(1)
        conn, addr = s.accept()
        msg = recvall(conn)
        logging.info("pulsar_pulse_classifier => Recieved msg from (%s, %s)",
                host,port)
        pulsar_pulse_candidate,to_save = cPickle.loads(msg)

        sn = float(pulsar_pulse_candidate['SN'])
        beam = int(pulsar_pulse_candidate['beam'])
        H_dm = float(pulsar_pulse_candidate['H_dm'])
        H_w = int(pulsar_pulse_candidate['H_w'])
        sample = int(pulsar_pulse_candidate['sample'])
        base_dir = FIL_FILE_DIR +'/BP'+str(THIS_BPNODE).zfill(2)+'/'+\
                utc.value+'/FB/BEAM_'+str(beam).zfill(3)
        search_dir = base_dir + '/'+utc.value+'.fil'
        #logging.info("pulsar_pulse_classifier => searching directory: %s" %search_dir)
        ftrs = get_features(beam,sample,sn,
                H_dm,H_w,search_dir)
        if not ftrs:
            logging.critical("pulsar_pulse_classifier => Got non features, skipping candidate")
            conn.sendall(str(-2)+"EOF")
            continue

        ftrs.utc = utc.value
        if not ftrs.isphonecall:
            classifier_input = sort_features(ftrs)
            _ , proba = classify(classifier_input,CLASSIFIER_THRESHOLD)
            logging.info("pulsar_pulse_classifier => pulse at sample %i classified, probability: %f",
                    pulsar_pulse_candidate['sample'],proba)
            conn.sendall(str(proba)+"EOF")
        else:
            logging.info("pulsar_pulse_classifier => pulse at sample %i is threshold phonecall",
                    pulsar_pulse_candidate['sample'])
            conn.sendall("-1"+"EOF")
        if to_save and FRB_DETECTOR_CFG["SAVE_FIL"].upper() in ["YES","TRUE"]:
            out_name = base_dir + "/PULSE_" + utc.value + "_" + str(sample) + ".fil"
            saveFBank(search_dir,sample,H_w,H_dm,out_name)

        conn.close()
            



def process_candidate(in_queue,utc,source_name,rfi_writer_queue,
        lock,training_file_dir):
    """ Processing function to be multiprocessed """
    logging.debug("%s Initiated, waiting for candidates" %os.getpid())
    global n_detect
    while True:
        candidate = in_queue.get()
        if candidate is None:
            logging.info("%s recieved a poison pill, terminating...",
                    os.getpid())
            break
        sn = float(candidate['SN'])
        beam = int(candidate['beam'])
        sample = int(candidate['sample'])
        base_dir = FIL_FILE_DIR+'/BP'+str(THIS_BPNODE).zfill(2)+'/'+\
                utc.value+'/'+source_name.value+'/BEAM_'+str(beam).zfill(3)
        search_dir = base_dir + '/'+utc.value+'.fil'
        logging.info('Searching directory: %s, candidate ==> sample: %s, width: %s',
                search_dir,candidate['sample'],candidate['H_w'])
#        file_directory = c_char_p(search_dir)        
        ftrs = get_features(beam,candidate['sample'],sn,
                candidate['H_dm'],candidate['H_w'],search_dir)
        if not ftrs:
            logging.info("Skipping candidate on sample = %s, and beam = %s",
                    candidate['sample'],beam)
            continue
        ftrs.utc = utc.value
#        lock.acquire()
#        logging.info('BP %s trying to write to training file',THIS_BPNODE)
#        c = str(output_l).strip("[]").replace(", "," ")+"\n"
#        logging.info(c)
#        training_file = open(training_file_dir,"a+")
#        training_file.write(c)
#        training_file.close()
#        lock.release()
#        continue
        if not ftrs.isphonecall:
            classifier_input = sort_features(ftrs)
            isFRB, proba = classify(classifier_input,CLASSIFIER_THRESHOLD)
            if isFRB:
                logging.info(str(proba*100)+"%% chance FRB! Beam: %i, "+\
                        "sample: %i",beam,candidate['sample'])
                cand_my_snr = my_snr(candidate['sample'],
                        candidate['H_dm'], candidate['H_w'], search_dir)
                inco_my_snr = get_snr_incoherent(candidate['sample'],
                        candidate['H_dm'],candidate['H_w'], utc.value,
                        search_dir)
                ratio = cand_my_snr / inco_my_snr

                logging.debug("SNR in current beam: %s",cand_my_snr)
                logging.debug("SNR in incoherent beam: %s",inco_my_snr)
                if DUMP_VOLTAGES:
                    if cand_my_snr > inco_my_snr:
                        obs_header = parse_cfg(FIL_FILE_DIR+'/BP'+str(THIS_BPNODE).zfill(2)+'/'+\
                                        utc.value+'/'+source_name.value+'/BEAM_'+str(beam).zfill(3)+\
                                        '/obs.header',['TSAMP'])
                        sampling_time = float(obs_header['TSAMP'])/10**6 # in seconds
                        send_dump_command(utc.value,sampling_time,
                                        candidate,ftrs,proba)
                        out_name = base_dir + "/FRB_" + utc.value + "_" + str(sample) + ".fil"
                        saveFBank(search_dir,sample,candidate['H_w'],
                                candidate['H_dm'],out_name)
                    if cand_my_snr > inco_my_snr:
                            rfi_writer_queue.put(ftrs.str_fmt("PULSES"))
                    else:
                            rfi_writer_queue.put(ftrs.str_fmt("THRSH"))
            else:
                logging.debug("Classified phone call: %i, %i, with probability: %f",
                        beam,candidate['sample'],proba)
                rfi_writer_queue.put(ftrs.str_fmt("RFI"))
                if len(os.listdir(base_dir)) < 6 and\
                        FRB_DETECTOR_CFG['SAVE_FIL'].upper() in ["YES","TRUE"]:
                    out_name = base_dir + "/RFI_" + utc.value + "_" + str(sample) + ".fil"
                    saveFBank(search_dir,sample,candidate['H_w'],
                            candidate['H_dm'],out_name)
        else:
            logging.debug("Phone call: %i, %i",beam,candidate['sample'])
            rfi_writer_queue.put(ftrs.str_fmt("RFI"))


def send_dump_command(utc,sampling_time,candidate,ftrs,proba):
    disp_delay = ((31.25*0.0000083*candidate['H_dm'])/pow(0.840,3))
    padded_delay = max(0.1*((2**ftrs.box/2)*sampling_time + disp_delay),0.2)
    time_sec1 = candidate['sample']*sampling_time -\
            padded_delay
    time_sec2 = candidate['sample']*sampling_time + disp_delay +\
            2*padded_delay
    fmt = "%Y-%m-%d-%H:%M:%S"
    fmtms = "%Y-%m-%d-%H:%M:%S.%f"
    cand_start_utc = datetime.datetime.strptime(utc,fmt) +\
            datetime.timedelta(seconds=time_sec1)
    cand_end_utc = datetime.datetime.strptime(utc,fmt) +\
            datetime.timedelta(seconds=time_sec2)
    cand_utc = datetime.datetime.strptime(utc,fmt) +\
            datetime.timedelta(seconds=candidate['sample']*sampling_time)
    dump_tag = Element('frb_detector_message')
    xml_cmd = SubElement(dump_tag,'cmd')
    xml_cmd.text = 'dump'
    xml_cand_start_utc = SubElement(dump_tag,'cand_start_utc')
    xml_cand_start_utc.text = datetime.datetime.strftime(\
            cand_start_utc,fmtms)[:-5]
    xml_cand_utc = SubElement(dump_tag,'cand_utc')
    xml_cand_utc.text = datetime.datetime.strftime(cand_utc,fmtms)[:-5]
    xml_cand_end_utc = SubElement(dump_tag,'cand_end_utc')
    xml_cand_end_utc.text = datetime.datetime.strftime(cand_end_utc,fmtms)[:-5]
    xml_cand_dm = SubElement(dump_tag,'cand_dm')
    xml_cand_dm.text = str(candidate['H_dm'])
    xml_cand_width = SubElement(dump_tag,'cand_width')
    xml_cand_width.attrib['units'] = 'seconds'
    xml_cand_width.text = str(2**ftrs.box*sampling_time)
    xml_beam_number = SubElement(dump_tag,'beam_number')
    xml_beam_number.text = str(candidate['beam'])
    xml_utc_start = SubElement(dump_tag,'utc_start')
    xml_utc_start.text = utc
    xml_snr = SubElement(dump_tag,'cand_snr')
    xml_snr.text = str(ftrs.sn)
    xml_probability = SubElement(dump_tag,'probability')
    xml_probability.text = str(proba)
    xml_cand_sample = SubElement(dump_tag,'cand_sample')
    xml_cand_sample.text = str(candidate['sample'])
    xml_cand_filter = SubElement(dump_tag,'cand_filter')
    xml_cand_filter.text = str(candidate['H_w'])
    xml_dump_msg = tostring(dump_tag,encoding='ISO-8859-1').replace("\n","")
    logging.info("Trying to send xml dump message to server")
    logging.info(xml_dump_msg)
    n_trials = 0
    for _ in range(3):
        try:
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
            s.connect((SERVER_HOST,DUMPPORT))
            s.send(xml_dump_msg)
            s.close()
            break
        except socket.error:
            logging.critical("Couldn't send message, trying again")
            n_trials += 1
            if n_trials != 3:
                time.sleep(0.2)
    if n_trials == 3:
        logging.critical("Couldn't send the message after 3 trials")
    else:
        logging.info("xml dump message sent")

def recvall(the_conn):
    total_data=[]
    while True:
        data = the_conn.recv(16384)
        if not data: break
        if data[-3:] == "EOF":
            data = data[:-3]
            total_data.append(data)
            break
        total_data.append(data)
    return ''.join(total_data)

def gracefull_file_close(f):
    logging.info("Closing file: %s",f.name)
    f.close()

def update_cfg(): 
    """
    Updates the cfg global variables
    """
    global CLASSIFIER_THRESHOLD, DUMP_VOLTAGES
    updated_cfg = parse_cfg(DADA_ROOT_SHARE+'frb_detector.cfg')
    
    if float(updated_cfg['CLASSIFIER_THRESHOLD']) != CLASSIFIER_THRESHOLD:
        CLASSIFIER_THRESHOLD = float(updated_cfg['CLASSIFIER_THRESHOLD'])
        logging.info("Config file updated, CLASSIFIER_THRESHOLD = %f",
                CLASSIFIER_THRESHOLD)
    if updated_cfg['DUMP_VOLTAGES'] == 'yes' and\
            DUMP_VOLTAGES == False:
        DUMP_VOLTAGES = True
        logging.info("Config file updated, switched dump voltages to true")
    if updated_cfg['DUMP_VOLTAGES'] == 'no' and\
            DUMP_VOLTAGES == True:
        DUMP_VOLTAGES = False
        logging.info("Config file update, switched dump voltages to false")


# GLOBALS
# -------
DADA_ROOT_SHARE = os.environ['DADA_ROOT']+'/share/'


# Loading config files, and initializing loggers
# ----------------------------------------------
FRB_DETECTOR_CFG = parse_cfg(DADA_ROOT_SHARE+'frb_detector.cfg')
#FRB_DETECTOR_CFG = parse_cfg('./frb_detector.cfg')
MOPSR_CFG_DIR = DADA_ROOT_SHARE+'mopsr.cfg'
MOPSR_BP_CFG_DIR = DADA_ROOT_SHARE+'mopsr_bp.cfg'
CORNERTURN_CFG_DIR = DADA_ROOT_SHARE+'mopsr_bp_cornerturn.cfg'
CORNERTURN_CFG = parse_cfg(CORNERTURN_CFG_DIR)
MOPSR_CFG = parse_cfg(MOPSR_CFG_DIR,["CLIENT_CONTROL_DIR","CLIENT_LOG_DIR",
"FRB_DETECTOR_BASEPORT","FRB_DETECTOR_DUMPPORT","SERVER_HOST",
"CLIENT_RECORDING_DIR","FRB_DETECTOR_INCOHERENT"])
SERVER_HOST = MOPSR_CFG["SERVER_HOST"]
BASEPORT = int(MOPSR_CFG["FRB_DETECTOR_BASEPORT"])

DUMPPORT = int(MOPSR_CFG["FRB_DETECTOR_DUMPPORT"])
INCOHERENT_BEAM_PORT = int(MOPSR_CFG[
    "FRB_DETECTOR_INCOHERENT"])
#INCOHERENT_BEAM_HOST = 'mpsr-bf00'
INCOHERENT_BEAM_HOST = CORNERTURN_CFG["RECV_0"]

CLASSIFIER_THRESHOLD = float(FRB_DETECTOR_CFG['CLASSIFIER_THRESHOLD'])
PULSAR_BASEPORT = int(FRB_DETECTOR_CFG['PULSAR_PULSE_CLASS_BASEPORT'])

if FRB_DETECTOR_CFG['DUMP_VOLTAGES'] == 'yes':
    DUMP_VOLTAGES = True
elif FRB_DETECTOR_CFG['DUMP_VOLTAGES'] == 'no':
    DUMP_VOLTAGES = False
FIL_FILE_DIR = MOPSR_CFG["CLIENT_RECORDING_DIR"]

# Wrapper functions initialization
# --------------------------------
#get_features = functions.get_features
#get_features.restype = CandidateFeatures

from helpers import get_features,get_feature_names

def main():
    # Parsing args
    # ------------
    parser = argparse.ArgumentParser(description='BP server that handles signals\
            from main server. Spawns ')
    parser.add_argument('bpnode', type=str, help='Number of the current BP\
            node running this instance')
#    parser.add_argument('--nproc', type=int, help ='Number of processes\
#            to spawn in each BP node for real time searching',
#            required = False, default = 4)
    parser.add_argument('--verbose','-v',action="store_true",help='Verbose\
            mode.')
    parser.add_argument('--test','-t',action='store_true',help='Dry run.\
            Doesn\'t log to system log files.')
    parser.add_argument('--daemonize','-d',action='store_false',
            help='Don\'t daemonize', default = True)
    args = parser.parse_args()

    verbose = args.verbose
    global THIS_BPNODE
    THIS_BPNODE = args.bpnode

    dry_run = args.test
#    n_processes = args.nproc
    daemon = args.daemonize

    n_processes = int(FRB_DETECTOR_CFG['N_PROCS'])

    if dry_run:
        client_ctrl_dir = FRB_DETECTOR_CFG['TEST_DIR']+'/control'
        client_log_dir = FRB_DETECTOR_CFG['TEST_DIR']+'/logs'
    else:
        client_ctrl_dir = MOPSR_CFG['CLIENT_CONTROL_DIR']
        client_log_dir = MOPSR_CFG['CLIENT_LOG_DIR']

    pid = os.getpid()
    script_name = os.path.basename(sys.argv[0]).lstrip("client_").\
            rstrip(".py")
    script_name_suffix = script_name + "_"+THIS_BPNODE
    logfile = client_log_dir+'/'+script_name_suffix+'.log'
    pidfile = client_ctrl_dir+'/'+script_name_suffix+'.pid'
    verbose = True
    if verbose:
        logging.basicConfig(filename=logfile,level=logging.DEBUG,
                format='(%(levelname)s): [%(asctime)s.%(msecs)03d]:'\
                        +'\t%(message)s',
                datefmt='%m-%d-%Y-%H:%M:%S')
    else:
        logging.basicConfig(filename=logfile,level=logging.INFO,
                format='(%(levelname)s) [%(asctime)s.%(msecs)03d]:'\
                        +'\t%(message)s',
                datefmt='%m-%d-%Y-%H:%M:%S')
    
    atexit.register(logging.shutdown)
    logging.info("BP master script initializing")
#    logging.info("Main thread pid: %s",pid)
    logging.info("Classifier threshold: %s",CLASSIFIER_THRESHOLD)
    if daemon:
        logging.info("Daemonizing")
        daemonize(pidfile, logfile)
    else:
        atexit.register(delpid,pidfile)
        pid = str(os.getpid())
        logging.debug("Writing pid file (pid %s)",pid)
        file(pidfile,'w+').write("%s\n" % pid)

    #Special process for the incoherent beam
    #---------------------------------------

    if int(THIS_BPNODE) == 0:
        inc_beam_process = Process(name = 'Incoherent beam process',
                target = incoherent_beam_proc, args = (INCOHERENT_BEAM_PORT,))
        inc_beam_process.start()
        atexit.register(inc_beam_process.terminate)
    
    controlThread = threading.Thread(name = 'controlThread',
            target = client_control_monitor,
            args=(client_ctrl_dir,script_name,THIS_BPNODE))
    controlThread.setDaemon(True)
    controlThread.start()

    lock = Lock()

    # Loading Classifier
    # -----------------
    logging.debug("Loading Classifier")
    global clf
    clf = joblib.load(FRB_DETECTOR_CFG['RANDOM_FOREST_FILE'])
    logging.debug("Classifier Loaded")
    
    # Spawning Processes
    # ------------------
    training_file_dir = "/home/wfarah/highres_test/feature_extractor/"+\
            "online_training_set/BP"+str(THIS_BPNODE)+".txt"
    if not os.path.exists(training_file_dir):
        hdr = get_feature_names()
        training_file = open(training_file_dir,"a+")
        training_file.write(hdr)
        training_file.close()
#    atexit.register(gracefull_file_close,training_file)
    logging.debug("Spawning "+str(n_processes)+" processes")
    in_queue = Queue()
    rfi_writer_queue = Queue()
    manager = Manager()
    utc = manager.Value(c_char_p,"")
    source_name = manager.Value(c_char_p,"")
    process_list = [Process(target = process_candidate, 
        args = (in_queue,utc,source_name,rfi_writer_queue,lock,
            training_file_dir)) for i in range(n_processes)]
    for proc in process_list:
        proc.start()
    time.sleep(0.5)
    
    monitorThread = threading.Thread(name = 'monitorThread',
                target=process_monitor_thread,
                args=(process_list,))
    monitorThread.setDaemon(True)
    monitorThread.start()
    
    writerThread = RFIWriterThread(THIS_BPNODE,rfi_writer_queue,
            name = 'writerThread')
    writerThread.setDaemon(True)
    writerThread.start()

    atexit.register(terminate_all,process_list,in_queue)
    atexit.register(writerThread.terminate_writer)

    # Creating Server Socket
    # ----------------------
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
    
    MOPSR_BP_CFG = parse_cfg(MOPSR_BP_CFG_DIR,['BP_'+str(THIS_BPNODE)])
    if dry_run:
#        host = socket.gethostname()
        host = MOPSR_BP_CFG['BP_'+str(THIS_BPNODE)]
    else:
        host = MOPSR_BP_CFG['BP_'+str(THIS_BPNODE)]
    logging.debug("Host name: %s, This BP node: %s",host,THIS_BPNODE)
    port_no = BASEPORT + 100 + (int(THIS_BPNODE)+1)
    pulsar_port_no = PULSAR_BASEPORT + 1 + (int(THIS_BPNODE)+1) 

    pulsar_classifier_proc = Process(target = pulsar_pulse_classifier,
            args = (utc,host,pulsar_port_no))
    pulsar_classifier_proc.start()
    atexit.register(pulsar_classifier_proc.terminate)
#    assert host == socket.gethostname().split(".")[0]
    s.bind((host,port_no))
    s.listen(10)
    logging.debug("listening for connection from: %s, %i",host,port_no)

    # Entering infinite loop
    # ----------------------
    t_old = time.time()
    while True:
        conn,addr = s.accept()
        from_srv0 = recvall(conn)
        if from_srv0[:3] == 'utc':
            update_cfg()
            if in_queue.qsize() != 0:
                logging.warning("Flushin candidates for new utc")
                while in_queue.qsize () != 0:
                    _ = in_queue.get(timeout=0.1)
            utc_str,source_str = from_srv0.split('/')
#            utc.value = from_srv0[4:]
            utc.value = utc_str[4:]
            source_name.value = source_str.split(':')[1]
            logging.debug("Acquired new utc: %s",utc)
            writerThread.change_file_name(utc_str[4:])
        elif from_srv0 == 'poison_pill':
            logging.debug("Poison_pill received, exiting")
            sys.exit(1)
        elif from_srv0 == 'STOP':
            logging.info('Observation stopped')
            if in_queue.qsize() != 0:
                logging.warning('%i candidates were flushed out after obs ended' %in_queue.qsize())
            while in_queue.qsize() != 0:
                _ = in_queue.get(timeout=0.1)
            writerThread.empty_queue()
            writerThread.close_file()
        else:
            t_new = time.time()
            candidate_list = cPickle.loads(from_srv0)
            if len(candidate_list) != 0:
                logging.debug("Acquired candidates, sending data to processing slaves")
                if t_new - t_old > 4 and in_queue.qsize() != 0:
                    logging.warning("Flushing candidates from last run")
                    while in_queue.qsize() != 0:
                        _ = in_queue.get()
                for candidate in candidate_list:
                    in_queue.put(candidate)
                t_old = t_new
            logging.debug("Candidates parsed, waiting for next round")
        conn.close()
    s.close()

if __name__ == "__main__":
    signal.signal(signal.SIGTERM,sigHandler)
    try:
        main()
    except SystemExit as e:
        if e.code == "Fork #1":
            logging.info("Exited from first parent")
        elif e.code == "Fork #2":
            logging.info("Exited from second parent")
        else:
            logging.exception("SystemExit")
    except:
        logging.exception("Exception")
