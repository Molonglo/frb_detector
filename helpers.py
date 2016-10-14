import os
import sys
import atexit
from signal import SIGTERM
import time
import logging
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment, tostring


"""
def parse_cfg(cfg_file,tags):
	config_dict = {}
	for tag in tags:
		config_dict[tag] = ["/home/wfarah/realtime_model"]
	return config_dict
"""

###############################################################################
#
# Turn the calling process into a daemon
#
def daemonize(pidfile, logfile):
# standard input will always be directed to /dev/null
	stdin = "/dev/null"
	stdout = logfile
	stderr = logfile

	try:
		pid = os.fork()
		if pid > 0:
			# exit first parent
			sys.exit(0)
	except OSError, e:
		sys.stderr.write("fork #1 failed: %d (%s)\n" % (e.errno, e.strerror))
		sys.exit(1)

# decouple from parent environment
	os.chdir("/")
	os.setsid()
	os.umask(0)

# do second fork
	try:
		pid = os.fork()
		if pid > 0:
# exit from second parent
			sys.exit(0)
	except OSError, e:
		sys.stderr.write("fork #2 failed: %d (%s)\n" % (e.errno, e.strerror))
		sys.exit(1)

# redirect standard file descriptors
	sys.stdout.flush()
	sys.stderr.flush()
	si = file(stdin, 'r')
	so = file(stdout, 'a+')
	se = file(stderr, 'a+', 0)
	os.dup2(si.fileno(), sys.stdin.fileno())
	os.dup2(so.fileno(), sys.stdout.fileno())
	os.dup2(se.fileno(), sys.stderr.fileno())
	atexit.register(delpid,pidfile)
# write pidfile, enable a function to cleanup pid file upon crash
	pid = str(os.getpid())
	logging.debug("Writing pid file (pid %s)",pid)
	file(pidfile,'w+').write("%s\n" % pid)




def parse_cfg(cfg_file,tags=None):
	"""Function that returns config file with given tags as dictionar

	Args:
		cfg_file (str): full directory to config file
		tags (list): list of tags to search the cgf_file

	Returns:
		config_dict (dict): dictionary with keys given in tags, and values
							extracted from cfg_file. If one tag doesn't exist,
							value corresponded will be None, else value is of
							type str, or list if multiple values exist for 
							same key.
	"""
	if tags == None:
		tags = []
		with open(cfg_file) as o:
			for line in o:
				if line[0] in ["\n","#"]: continue
				tags.append(line.split()[0])
	config_dict = {}
	with open(cfg_file) as o:
		for line in o:
			if line[0] in ["\n","#"]: continue
			for tag in tags:
				if tag in line:
					i = line.split()
					assert tag == i[0]
					config_dict[tag] = []
					for ii in i[1:]:
						if "#" in ii: break
						config_dict[tag].append(ii)
					if len(config_dict[tag]) == 1:
						config_dict[tag] = config_dict[tag][0]
					tags.remove(tag)
	for tag in tags:
		logging.warning("Couldn't parse <"+tag+"> from "+cfg_file)
		config_dict[tag] = None
	return config_dict


def control_monitor(control_dir,script_name):
	""" Function that writes pid into control folder, and constantly lists
		directory for .quit file, and kills the script
		
	Args:
		control_dir (str): control directory specified in mopsr.cfg
		script_name (str): name of script as it should show in ctrl direc
							(without extension)
		pid (int): process id

	"""
	while True:
		time.sleep(0.5)
		lst = os.listdir(control_dir)
		if script_name+".quit" in lst:
			logging.critical("Read .quit file, cleaning and exiting")
			stop_daemon(control_dir+"/"+script_name+".pid")


def client_control_monitor(control_dir,script_name,bf_numb):
	""" Same as control_monitor, but edited to check for the quit file
	    without the bp number as suffix """
	script_name_suffix = script_name+"_"+str(bf_numb)
	while True:
		time.sleep(0.5)
		lst = os.listdir(control_dir)
		if script_name+".quit" in lst or\
				script_name_suffix+".quit" in lst:
			logging.critical("Read .quit file, cleaning and exiting")
			stop_daemon(control_dir+"/"+script_name_suffix+".pid")


def delpid(pidfile):
	logging.critical("Deleting pidfile")
	os.remove(pidfile)

def stop_daemon(pidfile):
	pf = file(pidfile,'r')
	pid = int(pf.read().rstrip())
	pf.close()
	while True:
		logging.critical("Trying to kill %s ",pid)
		os.kill(pid,SIGTERM)
		time.sleep(0.1)

def sigHandler(signo,frame):
	logging.critical("%s Recieved a SIGTERM, cleaning...",os.getpid())
	sys.exit(1)

def create_xml_elem(msg_type,dump_dict=None):
	if msg_type == "ok":
		ok_tag = Element('frb_detector_message')
		child = SubElement(ok_tag, 'reply')
		child.text = "ok"
		return tostring(ok_tag,encoding ='ISO-8859-1').replace("\n","")
	elif msg_type == "fail":
		fail_tag = Element('frb_detector_message')
		child = SubElement(fail_tag, 'reply')
		child.text = "fail"
		return tostring(fail_tag,encoding='utf8').replace("\n","")
	elif msg_type == "dump":
		dump_tag = Element('frb_detector_message')
		child = SubElement(dump_tag,'dump')
		start_utc = SubElement(child,'start_utc')
		start_utc.text = dump_dict['start_utc']
		end_utc = SubElement(child,'end_utc')
		end_utc.text = dump_dict['end_utc']
		dm = SubElement(child,'dm')
		dm.text = dump_dict['dm']
		beam_number = SubElement(child,'beam_number')
		beam_number.text = dump_dict['beam_number']
		utc = SubElement(child,'utc')
		utc.text = dump_dict['utc']
		snr = SubElement(child,'snr')
		snr.text = dump_dict['snr']
		probability = SubElement(child,'probability')
		probability.text = dump_dict['probability']
		return tostring(dump_tag,encoding='ISO-8859-1').replace("\n","")
