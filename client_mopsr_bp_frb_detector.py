#!/home/dada/frb_detector/mpsr/bin/python

import os
from ctypes import *
functions = CDLL(os.environ['DADA_ROOT']+'/lib/'+'frb_detector_wrapper.so')
libc = CDLL("libc.so.6")

import numpy as np
import atexit
from multiprocessing import Queue,Process,Manager,Value
import threading
import time
import sys
sys.path.append(os.environ['DADA_ROOT']+'/lib/')
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import logging
import socket
import cPickle
import argparse
from helpers import parse_cfg,client_control_monitor,daemonize,sigHandler
from helpers import create_xml_elem,delpid
import signal
import datetime

from sklearn.externals import joblib

class CandidateFeatures(Structure):
	_fields_ = [("sn",c_float),("sn_0",c_float),("width",c_int),
			("nstart",c_int),("nend",c_int),("F1",c_float),
			("F2",c_float),("F3",c_float),("n_rms_mask",c_int),
			("sn_rms",c_float),("mod_ind",c_float),
			("mod_indT",c_float), ("isphonecall",c_int)]


def sort_features(ftrs):
	""" Function that takes in CandidateFeatures object, and returns a numpy
	array that serves as an input for the classifier"""
	sorted_features = np.array([ftrs.width,ftrs.sn/ftrs.sn_0,ftrs.F1,ftrs.F2,
		ftrs.F3,ftrs.sn_rms,ftrs.n_rms_mask,ftrs.mod_ind,ftrs.mod_indT])
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
	y=clf.predict_proba(features)
	y=y[0]
	if y[0]<threshold and y[1]>threshold:
		return True,y[1]
	return False,y[1]

def process_monitor_thread(process_list,refresh_time=20):
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
			sys.exit(1)
		time.sleep(refresh_time)


def terminate_all(n_proc,in_queue):
	logging.critical("Terminating sub-processes")
	for i in range(n_proc):
#		Poison pill approach
		in_queue.put(None)


def process_candidate(in_queue,utc,source_name):
	""" Processing function to be multiprocessed """
	logging.debug("%s Initiated, waiting for candidates" %os.getpid())
	global n_detect
	while True:
		candidate = in_queue.get()
		if candidate is None:
			logging.info("%s recieved a poison pill, terminating...",
					os.getpid())
			break
		beam = int(candidate['beam'])
		H_dm = c_float(candidate['H_dm'])
		H_w = c_int(candidate['H_w'])
		time_sample = c_int(candidate['sample'])
		search_dir = FIL_FILE_DIR+'/'+utc.value+'/'+source_name.value+\
				'/BEAM_'+str(beam).zfill(3)+'/'+utc.value+'.fil'
		logging.info('Searching directory: %s',search_dir)
		file_directory = c_char_p(search_dir)
		ftrs = get_features(time_sample,H_dm,H_w,file_directory)
		if not ftrs.isphonecall:
			classifier_input = sort_features(ftrs)
			isFRB, proba = classify(classifier_input,CLASSIFIER_THRESHOLD)
			if isFRB:
				logging.info(str(proba*100)+" chance FRB! Beam: %i, "+\
						"sample: %i",beam,candidate['sample'])
				if DUMP_VOLTAGES:
					obs_header = parse_cfg(FIL_FILE_DIR+'/'+utc.value+'/'+\
							'obs.header',['TSAMP'])
					sampling_time = float(obs_header['TSAMP'])/10**6 # in seconds
					send_dump_command(utc.value,sampling_time,
							candidate,ftrs,proba)
			else:
				logging.debug("Classified phone call: %i, %i",
						beam,candidate[0])
		else:
			logging.debug("Phone call: %i, %i",beam,candidate[0])



def send_dump_command(utc,sampling_time,candidate,ftrs,proba):
	disp_delay = ((31.25*0.0000083*candidate['H_dm'])/pow(0.840,3))
	time_sec1 = candidate['sample']*sampling_time -\
			max(((ftrs.width/2)*sampling_time + disp_delay),0.5)
	time_sec2 = candidate['sample']*sampling_time + disp_delay +\
			max(((ftrs.width/2)*sampling_time + disp_delay),0.5)
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
	xml_cand_width.text = str(ftrs.width*sampling_time)
	xml_beam_number = SubElement(dump_tag,'beam_number')
	xml_beam_number.text = str(candidate['beam'])
	xml_utc_start = SubElement(dump_tag,'utc_start')
	xml_utc_start.text = utc
	xml_snr = SubElement(dump_tag,'cand_snr')
	xml_snr.text = str(ftrs.sn)
	xml_probability = SubElement(dump_tag,'probability')
	xml_probability.text = str(proba)
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
		total_data.append(data)
	return ''.join(total_data)



# GLOBALS
# -------
DADA_ROOT_SHARE = os.environ['DADA_ROOT']+'/share/'


# Loading config files, and initializing loggers
# ----------------------------------------------
FRB_DETECTOR_CFG = parse_cfg(DADA_ROOT_SHARE+'frb_detector.cfg')
#FRB_DETECTOR_CFG = parse_cfg('./frb_detector.cfg')
MOPSR_CFG_DIR = DADA_ROOT_SHARE+'mopsr.cfg'
MOPSR_BP_CFG_DIR = DADA_ROOT_SHARE+'mopsr_bp.cfg'
MOPSR_CFG = parse_cfg(MOPSR_CFG_DIR,["CLIENT_CONTROL_DIR","CLIENT_LOG_DIR",
"FRB_DETECTOR_BASEPORT","FRB_DETECTOR_DUMPPORT","SERVER_HOST",
"CLIENT_RECORDING_DIR"])
SERVER_HOST = MOPSR_CFG["SERVER_HOST"]
BASEPORT = int(MOPSR_CFG["FRB_DETECTOR_BASEPORT"])
DUMPPORT = int(MOPSR_CFG["FRB_DETECTOR_DUMPPORT"])

CLASSIFIER_THRESHOLD = float(FRB_DETECTOR_CFG['CLASSIFIER_THRESHOLD'])
if FRB_DETECTOR_CFG['DUMP_VOLTAGES'] == 'yes':
	DUMP_VOLTAGES = True
elif FRB_DETECTOR_CFG['DUMP_VOLTAGES'] == 'no':
	DUMP_VOLTAGES = False
FIL_FILE_DIR = MOPSR_CFG["CLIENT_RECORDING_DIR"]

# Wrapper functions initialization
# --------------------------------
get_features = functions.get_features
get_features.restype = CandidateFeatures

def main():
	# Parsing args
	# ------------
	parser = argparse.ArgumentParser(description='BF server that handles signals\
			from main server. Spawns ')
	parser.add_argument('bfnode', type=str, help='Number of the current BF\
			node running this instance')
#	parser.add_argument('--nproc', type=int, help ='Number of processes\
#			to spawn in each BF node for real time searching',
#			required = False, default = 4)
	parser.add_argument('--verbose','-v',action="store_true",help='Verbose\
			mode.')
	parser.add_argument('--test','-t',action='store_true',help='Dry run.\
			Doesn\'t log to system log files.')
	parser.add_argument('--daemonize','-d',action='store_false',
			help='Don\'t daemonize', default = True)
	args = parser.parse_args()

	verbose = args.verbose
	bfnode_numb = args.bfnode
	dry_run = args.test
#	n_processes = args.nproc
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
	script_name_suffix = script_name + "_"+bfnode_numb
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
	logging.info("BF master script initializing")
	logging.info("Main thread pid: %s",pid)
	logging.info("Classifier threshold: %s",CLASSIFIER_THRESHOLD)
	if daemon:
		daemonize(pidfile, logfile)
	else:
		atexit.register(delpid,pidfile)
		pid = str(os.getpid())
		logging.debug("Writing pid file (pid %s)",pid)
		file(pidfile,'w+').write("%s\n" % pid)

	controlThread = threading.Thread(name = 'controlThread',
			target = client_control_monitor,
			args=(client_ctrl_dir,script_name,bfnode_numb))
	controlThread.setDaemon(True)
	controlThread.start()
	
	# Loading Classifier
	# -----------------
	logging.debug("Loading Classifier")
	global clf
	clf = joblib.load(FRB_DETECTOR_CFG['RANDOM_FOREST_FILE'])
	logging.debug("Classifier Loaded")
	
	# Spawning Processes
	# ------------------
	logging.debug("Spawning "+str(n_processes)+" processes")
	in_queue = Queue()
	manager = Manager()
	utc = manager.Value(c_char_p,"")
	source_name = manager.Value(c_char_p,"")
	process_list = [Process(target = process_candidate, 
		args = (in_queue,utc,source_name)) for i in range(n_processes)]
	for proc in process_list:
		proc.start()
	time.sleep(0.5)
	monitorThread = threading.Thread(name = 'monitorThread',
                target=process_monitor_thread,
                args=(process_list,))
	monitorThread.setDaemon(True)
	monitorThread.start()
	atexit.register(terminate_all,n_processes,in_queue)

	# Creating Server Socket
	# ----------------------
	logging.debug("Creating Socket")
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
	
	MOPSR_BP_CFG = parse_cfg(MOPSR_BP_CFG_DIR,['BP_'+str(bfnode_numb)])
	if dry_run:
#		host = socket.gethostname()
		host = MOPSR_BP_CFG['BP_'+str(bfnode_numb)]
	else:
		host = MOPSR_BP_CFG['BP_'+str(bfnode_numb)]
	logging.debug("Host name: %s",host)
	port_no = BASEPORT + 100 + (int(bfnode_numb)+1)
	assert host == socket.gethostname().split(".")[0]
	s.bind((host,port_no))
	s.listen(10)
	logging.debug("listening to: %s, %i",host,port_no)

	# Entering infinite loop
	# ----------------------
	t_old = time.time()
	while True:
		conn,addr = s.accept()
		from_srv0 = recvall(conn)
		if from_srv0[:3] == 'utc':
			if in_queue.qsize() != 0:
				logging.warning("Flushin candidates for new utc")
				while in_queue.qsize () != 0:
					_ = in_queue.get()
			utc_str,source_str = from_srv0.split('/')
#			utc.value = from_srv0[4:]
			utc.value = utc_str[4:]
			source_name.value = source_str.split(':')[1]
			logging.debug("Acquired new utc: %s",utc)
		elif from_srv0 == 'poison_pill':
			logging.debug("Poison_pill received, exiting")
			sys.exit(1)
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
	except:
		logging.exception("")
