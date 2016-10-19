#!/usr/bin/env python

import numpy as np
from numpy.lib.recfunctions import rec_append_fields
import threading,Queue
import logging
import sys
import socket
import os
import signal
import atexit
from time import sleep
import datetime
import argparse
import time
import cPickle
from subprocess import PIPE, Popen
from helpers import parse_cfg,control_monitor,create_xml_elem
from helpers import daemonize,sigHandler,delpid
import ephem
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment, tostring




def year_fraction(date):
	""" Function that returns the year fraction for the input date

	Args:
		date (datetime): Datetime object

	"""
	start = datetime.date(date.year,1,1).toordinal()
	year_length = datetime.date(date.year+1,1,1).toordinal() - start
	return date.year + float(date.toordinal() - start) / year_length


def get_pulsars_db(db_directory):
	""" Function that returns pulsars database, from the existing psrcat.dat
		file.

	Args:
		db_directory (str): Full directory to database

	Returns:
		pulsar_db (np.array): Recarray containing the fields:
							  'NAME','RAJ','DECJ','DM','PREC_RA',
							  'PREC_DEC','HA_COR'
							  Where: prec stands for precessed for current
							  year, and HA_COR is the correction on the HA
	"""
	dtype = {'names':('NAME','RAJ','DECJ','DM'),'formats':('S10','S16','S16','f4')}
	pulsar_db = np.loadtxt(db_directory,dtype=dtype)
	current_year = year_fraction(datetime.datetime.today())
	ha_corrections = []
	prec_ra = []
	prec_dec = []
	for pulsar in pulsar_db:
		# Precessing coordinates
		pulsar_ephem = ephem.readdb(pulsar['NAME']+",f|L,"+pulsar['RAJ']+","+
				pulsar['DECJ']+",2000")
		pulsar_ephem.compute(epoch = str(current_year))
		prec_ra.append(pulsar_ephem.a_ra)
		prec_dec.append(pulsar_ephem.a_dec)
		ha_corrections.append(get_ha_cor(pulsar_ephem.a_dec))
	pulsar_db = rec_append_fields(pulsar_db,'PREC_RA',prec_ra)
	pulsar_db = rec_append_fields(pulsar_db,'PREC_DEC',prec_dec)
	pulsar_db = rec_append_fields(pulsar_db,'HA_COR',ha_corrections)
	return pulsar_db


def get_ha_cor(dec):
	""" Function that returns the correction on the RA,
		due to the skewness and slope of telescope
	Args:
		dec (ephem.degrees): dec of pulsar
	
	Returns:
		added_factor (ephem.degrees): Correction on HA
	"""
	num = np.sin(dec - lmbda)*np.sin(eta) - np.cos(dec-lmbda)*np.cos(eta)*np.sin(zeta)
	denom = -np.sin(dec-lmbda)*np.sin(lmbda)*np.cos(zeta) + np.cos(dec-lmbda)*(np.cos(lmbda)*np.cos(zeta)-np.sin(lmbda)*np.sin(zeta)*np.sin(eta))
	added_factor = np.arctan(num/denom)
	return added_factor


def _cmdline(command):
	"""Function that captures output from screen
	"""
	process = Popen(
			args=command,
			stdout=PIPE,
			shell=True
			)
	return process.communicate()[0]


def _get_nsmd(utc,ra,dec):
	"""Function that returns ns and md coordinates of boresight
	"""
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.connect(("localhost",nsmd_port))
	s.sendall(utc+","+ra+","+dec)
	ns,md = [float(i) for i in (s.recv(32).rstrip("\x00").split(","))]
	s.close()
	return ns,md
	#dec="-- "+dec[:]	   #This is done for the mopsr_getns command. 
	#cmd_getns="mopsr_getns "+utc+" "+ra+" "+dec
	#cmd_getmd="mopsr_getmd "+utc+" "+ra+" "+dec
	#print cmd_getns,cmd_getmd
	#return float(_cmdline(cmd_getns)),float(_cmdline(cmd_getmd)) #captures mopsr_getns output from screen



def refine_pulsar_db(mode,pulsar_db,**kwargs):
	""" Function that refines the bright pulsars according to their
		proximity to boresight. If mode is "Tracking", returns all
		pulsars within 10 deg (in ra and dec). If mode is "Stationary",
		returns all pulsars with dec +- 5 (for any ra).

	Args:
		mode (str): "TRACKING" or "STATIONARY"
		pulsar_db (np.ndarray): Array containing pulsar names, coordinates,
								and DM
		**kwargs: provide boresight_ra, boresight_dec for "Tracking"
				  provide boresight_dec for "Stationary"

	Returns:
		refined_pulsar_db (np.ndarray): A subset of pulsar_db

	"""
	if mode == "TRACKING":
		ra_lim = 10 #deg
		dec_lim = 10 #deg
		e_bs_dec = ephem.degrees(kwargs['boresight_dec'])
		e_bs_ra = ephem.hours(kwargs['boresight_ra'])
		i = 0
		ind_list = []
		for pulsar in pulsar_db:
			e_pulsar_dec = ephem.degrees(pulsar['DECJ'])
			e_pulsar_ra = ephem.hours(pulsar['RAJ'])
			dec_dif = np.abs(e_pulsar_dec - e_bs_dec) * 180/np.pi
			ra_dif = np.abs(e_pulsar_ra - e_bs_ra) * 180/np.pi
			if (ra_dif < ra_lim or ra_dif > (360 - ra_lim)) and dec_dif < dec_lim:
				ind_list.append(i)
			i+=1
		return pulsar_db[ind_list]
	elif mode == "STATIONARY":
		dec_lim = 5
		e_bs_dec = ephem.degrees(kwargs['boresight_dec'])
		i = 0
		ind_list = []
		for pulsar in pulsar_db:
			e_pulsar_dec = ephem.degrees(pulsar['DECJ'])
			dec_dif = np.abs(e_pulsar_dec - e_bs_dec) * 180/np.pi
			if dec_dif < dec_lim:
				ind_list.append(i)
			i+=1
		return pulsar_db[ind_list]



def getPotentialPulsars_tracking(utc,boresight_ra,boresight_dec,n_beams,pulsar_db):
	"""Function that computes a list of potential pulsar in fanbeam space
	
	Args:
		utc (str): Current UTC 
		boresight_ra (str): RA of boresight
		boresight_dec (str): Dec of boresight
		n_beams (int): Number of fanbeams
	
	Returns:
		np.array: Potential Pulsars with added FB column.

	"""
	telescope_ns,telescope_md = _get_nsmd(utc,boresight_ra,boresight_dec)
	index_to_flag = []	#index of pulsars from pulsars_list to flag
	estimated_fb = []		#Fanbeam where a pulsar is estimated to be located
	i=0
	# Pulsar database is: col 0 being name, col 1 being RA, col2
	# being DEC, col 3 being DM
	for pulsar in pulsar_db:
		#print PSRJ,RAJ,DECJ
		pulsar_ns,pulsar_md = _get_nsmd(utc,pulsar['RAJ'],pulsar['DECJ'])
		if np.abs(telescope_ns-pulsar_ns)<ns_threshold and np.abs(telescope_md-pulsar_md)<md_threshold:
			index_to_flag.append(i)
			estimated_fb.append(n_beams*(((pulsar_md-telescope_md)/(2.0/np.cos(np.radians((pulsar_md+telescope_md)/2)))+1)/2.0) + 1)
		i+=1
	l = len(estimated_fb)
	if l == 0:
		return []
	if l == 1:
		estimated_fb = [estimated_fb]
		pulsars_to_flag = [pulsar_db[index_to_flag]]
	else:
		pulsars_to_flag = pulsar_db[index_to_flag]
	pulsar_list = rec_append_fields(pulsars_to_flag,('FB'),estimated_fb,dtypes='f4')
	return pulsar_list


def getPotentialPulsars_stationary(ref_pulsar_db):
	index_to_flag = []
	estimated_fb = []
	i=0
	for pulsar in ref_pulsar_db:
		sec_dif = (molonglo.sidereal_time() - pulsar['PREC_RA'] - pulsar['HA_COR'])*R2S
		FB = sec_dif / (time_cover_beam/np.cos(pulsar['PREC_DEC'])/351) + 177
		if FB < 352+3 and FB > -2:
			index_to_flag.append(i)
			estimated_fb.append(FB)
		i+=1
	l = len(estimated_fb)
	if l == 0:
		return []
	if l == 1:
		estimated_fb = [estimated_fb]
		pulsars_to_flag =[ref_pulsar_db[index_to_flag]]
	else:
		pulsars_to_flag = ref_pulsar_db[index_to_flag]
	pulsar_list = rec_append_fields(pulsars_to_flag,('FB'),estimated_fb,dtypes='f4')
	return pulsar_list

	

def candidateFilter(candidates,PulsarList,threshold_filter=True):
	"""Function that performs stage 1 masking
	
	Returns:
			np.array: Array of the good candidates

	"""
	if threshold_filter:
		good_candidates = thresholdFilter(candidates)
	else:
		good_candidates = candidates
	
	if len(PulsarList) == 0:
		""" No pulsars in FOV """
		return good_candidates
	mask = np.ones(len(good_candidates),np.bool)
	# if pulsar is present, candidate equivalent index in mask will 
	# be replaced by False
	i = 0
	for candidate in good_candidates:
		beam, H_dm = int(candidate['beam']), candidate['H_dm']
		for Pulsar in PulsarList:
			if candidateIsPulsar(beam,H_dm,Pulsar):
				logging.info("%s detected at fanbeam: %i with dm: %f at time: %f and sample: %i",Pulsar['NAME'],
						beam,H_dm,candidate['time'],candidate['sample'])
				mask[i] = False
			else:
				pass
		i+=1
	return good_candidates[mask] #Great Candidates


def thresholdFilter(candidates):
	"""Function that returns candidates HEIMDAL's output above threshold
	
	Args:
		candidates (np.array): output line from HEIMDAL
		(None): If no candidates exist
	
	"""
	logging.info("Received %s candidates",candidates.size)
	a = np.where(candidates['H_w'] < boxcar_threshold)[0]
	b = np.where(candidates['H_dm'] > dm_threshold)[0]
	mask = np.where((candidates['H_w'] < boxcar_threshold) & (candidates['H_dm'] > dm_threshold))[0]
	logging.info("%s events have widths < than %s",len(a),boxcar_threshold)
	logging.info("%s events have dm > %s",len(b),dm_threshold)
	logging.info("%s events that passed both criteria",len(mask))
	return candidates[mask]


def candidateIsPulsar(beam,H_dm,pulsar):
	"""Function that returns bool, whether pulsar in FB or not """
	if pulsar['NAME'] is 'J0835-4510' and H_dm < 100:
		""" Vela alert! Discard all candidates with DM<100"""
		return True
	if (beam >= (pulsar['FB']-2) and beam <= (pulsar['FB']+2)) and (H_dm<1.2*pulsar['DM'] and H_dm>0.8*pulsar['DM']):
		return True
	else:
		return False


def updatePulsarList(pulsar_list_queue,obsInfoQueue,pulsar_refresh_time,
			n_beams):
	""" Function to be executed by a thread, computes whether pulsar might be visible in a fanbeam.
		Goes through an infinite loop, executes every 'pulsar_refresh_time', and puts output in queue,
		(replaces in case queue is not empty)
		Args:
			queue (Queue.Queue instance): Queue that holds pulsar object
			pulsar_refresh_time (int): Wait time for 
			boresight_ra (float): RA of boresight
			boresight_dec (float): Dec of boresight
			n_beams (int): Number of fanbeams
			mode (Str): 'TRACKING' or 'STATIONARY'
		Returns:
			(None)
			
	"""
	logging.debug("Pulsar List Thread spawned")
	obsInfoDict = obsInfoQueue.get() #Blocks for the first obs, then goes in infinite loop
	boresight_ra = obsInfoDict['boresight_ra']
	boresight_dec = obsInfoDict['boresight_dec']
	pulsar_db = obsInfoDict['pulsar_db']
	while True:
		sleep(pulsar_refresh_time)
		logging.debug("Updating pulsar list")
		if not obsInfoQueue.empty():
			obsInfoDict = obsInfoQueue.get()
			boresight_ra = obsInfoDict['boresight_ra']
			boresight_dec = obsInfoDict['boresight_dec']
			pulsar_db = obsInfoDict['pulsar_db']
		utc = datetime.datetime.utcnow()
		utc = datetime.datetimeToStr(utc)
		if queue.empty():   #Ensures only 1 pulsar list in queue
			pulsar_list = getPotentialPulsars_tracking(utc,boresight_ra,boresight_dec,n_beams,pulsar_db)
			queue.put(pulsar_list)
			queue.task_done()
		else:
			pulsar_list = getPotentialPulsars_tracking(utc,boresight_ra,boresight_dec,n_beams,pulsar_db)
			_ = queue.get()
			queue.put(pulsar_list)
			queue.task_done()
		print "task done, data in queue"
		#pulsar_list = getPotentialPulsars(utc,boresight_ra,boresight_dec,n_beams,mode)


def datetimeToStr(datetime_utc):
	"""Function that convert datetime object to str with specific format"""
	fmt = "%Y-%m-%d-%H:%M:%S"
	return datetime.datetime.strftime(datetime_utc,fmt)



def getObsInfo(utc):
	"""Funtion that returns obs.info

	Args:
		utc (str): format yyyy-mm-dd-hh-mm-ss should be the start UTC of observation!

	Returns:
		obs_info (dict): Dictionary with keys being the first column of 
						 obs.info file, and values the second column
	"""
	obs_info={}
	with open(filterbank_directory+utc+"/obs.info") as o:
		for line in o:
			if line[0] not in ["\n","#"]:
				i = line.split(" ")
				if i[-1].rstrip():
					if i[-2] == "FRB" and i[-1]=="Transit\n":
						obs_info[i[0]] = "FRB Transit"
					else:
						obs_info[i[0]] = i[-1].rstrip()
				else:
					obs_info[i[0]] = " "
	return obs_info


def realTimePlot(queue):
	plot_array = []
	eps=10**(-2)
	plb.ion()
	while True:
		sleep(0.05)
		plot_array.append(queue.get()+eps)
		plb.clf()
		plb.bar(range(1,1+len(plot_array)),plot_array)
		plb.title("Candidate Rate")
		#plb.ylim(-1,np.max(PlotArray))
		plb.draw()


def send_cands(addr,candidates):
	"""Sends candidates through socket

	Args:
		addr (tuple): (IP,port)
		candidates (list): list of lists Heimdal candidates

	"""
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	try:
		s.connect(addr)
		p = cPickle.dumps(candidates[['sample','H_dm','H_w','beam']],protocol=2)
		s.sendall(p)
		s.close()
		logging.info("Successfully sent candidates to (%s, %i)",addr[0],addr[1])
	except socket.error:
		logging.critical("Trying to send candidates. Connection to: (%s, %i) refused",addr[0],addr[1])


def send_cands_to_bf(bf_addrs,beam_config,filtered_candidates):
	"""Function that sends candidate list to the corresponding BP nodes

	Args:
		bf_addrs (dict): dictionary with keys as 'RECV_0', 'RECV_1',... and values
							as tuples of (IP,port_number)
		beam_config (dict): dictionary, output of parse_beams_config()
		filtered_candidates (list): List of list of candidates to be processed
									by the BP nodes
	"""
	logging.debug("Sending Candidates")
	candidate_dict = {}
	n_bf = len(bf_addrs)
	for i in range(n_bf):
		bf_node = 'RECV_'+str(i)
		candidate_dict[bf_node] = []
	for index,candidate in zip(xrange(len(filtered_candidates)),filtered_candidates):
		#NOTE: This loop is time inefficient
		beam = candidate['beam']
		for i in range(n_bf):
			if beam >= beam_config["BEAM_FIRST_RECV_"+str(i)] and beam <= beam_config["BEAM_LAST_RECV_"+str(i)]:
				candidate_dict['RECV_'+str(i)].append(index)
				break
	for i in range(n_bf):
		if len(candidate_dict['RECV_'+str(i)]) != 0:
			send_cands(bf_addrs['RECV_'+str(i)],filtered_candidates[candidate_dict['RECV_'+str(i)]])

	
def parse_bf_config(config_dir):
	""" Function that reads configuration file for the beams, in the shared directory

	Args:
		config_dir (str): directory to the config file,
							usually /home/dada/linux_64/share/mopsr_bp_cornerturn.cfg
	
	Returns:
		beam_config (dict): dictionary that defines the start/end beam for each bf node
							keys as 'BEAM_FIRST_RECV_0','BEAM_LAST_RECV_1',...
							values are int 
		bf_ips (dict): dictionary with keys as 'RECV_0', 'RECV_1',...
						and values representing the addresses of bf nodes
						as 'mpsr-bf00', 'mpsr-bf01',...
	NOTE: incremented the beam number by 1. (ex: starting 1 ending 352)
	"""
	beam_config = {}
	bf_ips = {}
	with open(config_dir) as o:
		logging.debug("Opening config file: "+config_dir)
		for line in o:
			if line[:6] == "NBEAM ":
				line = line.split()
				beam_config[line[0]] = int(line[1])
				logging.debug("Parsing in "+line[0]+": "+line[1])
			if line[:15] in ["BEAM_FIRST_RECV","BEAM_LAST_RECV_"]:
				line = line.split()
				beam_config[line[0]] = int(line[1]) + 1 #NOTE: Read docstring for increment
				logging.debug("Parsing in "+line[0]+": "+line[1])
			if line[:5] == "RECV_":
				line = line.split()
				bf_ips[line[0]] = line[1]
				logging.debug("Parsing in "+line[0]+": "+line[1])
	return beam_config,bf_ips


def test_socket_listen(host):
	"""For testing"""
	logging.debug("*Test* Listening to backend (%s, %i)",host[0],host[1])
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
	s.bind(host)
	s.listen(5)
	conn,addr = s.accept()
	t = time.time()
	msg = recvall(conn)
	#msg = conn.recv(16384)
	s.close()
	[flag,data] = cPickle.loads(msg)
	print "time to recv: ",time.time() - t
	return flag,data


def recvall(the_conn):
	total_data=[]
	while True:
		data = the_conn.recv(256)
		logging.debug('recieved %s number of bytes:',len(data))
		logging.debug('%s',data)
		if not data: break
		total_data.append(data)
	return ''.join(total_data)

def socket_listen(host):
	"""Function that listens to socket
		
	Args:
		host (tuple): Tuple containing hostname, and port number
						ex: host = ('localhost',23456)
	
	Return:
		flag (str): flag of what the type of signal is.
		
		data ('any type'): If flag is 'START', data is obsInfo (dict)
							If flag is 'CANDIDATES', data is (np.ndarray)
							If flag is 'STOP', data is (None)
							IF flag is 'UNKOWN_FLAG', data is (None)
	"""
	logging.debug("Listening to backend (%s, %i)",host[0],host[1])
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR,1)
	s.bind(host)
	s.listen(2)
	conn,addr = s.accept()
	tmp = conn.recv(4096)
	msg = tmp
	while not msg.endswith("</frb_detector_message>\r\n"):
		logging.debug(tmp)
		tmp = conn.recv(4096)
		msg += tmp

	logging.debug(msg)
	flag,data = translate_xml(msg)
	if flag == "UNKOWN_FLAG":
		logging.critical("Recieved an UNKOWN FLAG from backend")
		conn.send(fail_tag)
	else:
		conn.send(ok_tag)
	conn.close()
	s.close()
	return flag,data

def create_xml_elem(msg_type,dump_dict=None):
	if msg_type == "ok":
		ok_tag = Element('frb_detector_message')
		child = SubElement(ok_tag, 'reply')
		child.text = "ok"
		msg = "<?xml version='1.0' encoding='ISO-8859-1'?>"+\
				"<frb_detector_message><reply>ok</reply></frb_detector_message>"
		return msg
#		return tostring(ok_tag,encoding ='ISO-8859-1')
	elif msg_type == "fail":
		fail_tag = Element('frb_detector_message')
		child = SubElement(fail_tag, 'reply')
		child.text = "fail"
		return tostring(fail_tag,encoding='ISO-8859-1')


def translate_xml(msg):
	""" Function that reads message received from socket, and translates
		it from xml into python data type
	"""
	root = ET.fromstring(msg)
	assert root.tag == "frb_detector_message"
	flag = root.find('cmd').text.upper()
	if flag == 'START':
		xml_tags = get_xml_tags(flag)
		obsInfo = {}
		for tag,fmt in xml_tags:
			obsInfo[tag.upper()] = fmt(root.find(tag).text)
		return flag,obsInfo
	elif flag == 'CANDIDATES':
		utc = root.find('utc_start').text
		elem_list = root.findall('cand')
		heimdal_candidates = []
		for elem in elem_list:
			line = elem.text.split()
			cand = ()
			for i in [0,1,2,3,5,12]:
				cand += (line[i],)
			heimdal_candidates.append(cand)
		heimdal_candidates = np.array(heimdal_candidates,dtype=cand_dtype)
		return flag,heimdal_candidates
	elif flag == 'STOP':
		return flag,None
	else:
		logging.critical("Unable to translate xml msg recieved")
		flag = 'UNKNOWN_FLAG'
		return flag,None

def get_xml_tags(flag):
	""" Returns a list of lists containing xml tag and data type"""
	if flag == 'START':
		xml_tags = [["utc_start",str],["source",str],["ra",str],["dec",str],
				["md_angle",float],["ns_tilt",float],["pid",str],["mode",str],
				["config",str],["observing_type",str],["nchan",int],
				["nbit",int],["tsamp",float]]
		return xml_tags
	else:
		raise("Unkown Flag")

def send_utc_to_bf(start_utc,source_name,bf_addrs):
	""" Sends utc to all bf nodes.
	
	Args:
		start_utc (str): utc with format 'YYYY-MM-DD-HH:MM:SS'
		bf_addrs (dict): dictionary with keys as 'RECV_0', 'RECV_1',... and values
							as tuples of (IP,port_number)
	
	"""
	n_bf = len(bf_addrs)
	for bf_node,addr in bf_addrs.iteritems():
		s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		try:
			s.connect(addr)
			msg = 'utc:'+start_utc+"/source:"+source_name
			s.sendall(msg)
			s.close()
			logging.debug("Connection established, %s sent to: (%s, %i)",msg,addr[0],addr[1])
		except socket.error:
			logging.critical("Connection to: (%s, %i) refused, couldn't send utc to %s",addr[0],addr[1],bf_node)





# --------
# Globals:
# --------
DADA_ROOT_SHARE = os.environ['DADA_ROOT']+'/share/'
FRB_DETECTOR_CFG = parse_cfg(DADA_ROOT_SHARE+'frb_detector.cfg')
#FRB_DETECTOR_CFG = parse_cfg('./frb_detector.cfg')
MOPSR_CFG_DIR = DADA_ROOT_SHARE+'mopsr.cfg'
MOPSR_BP_CFG_DIR = DADA_ROOT_SHARE+'mopsr_bp.cfg'
CORNERTURN_CFG_DIR = DADA_ROOT_SHARE+'mopsr_bp_cornerturn.cfg'


pulsar_refresh_time = float(FRB_DETECTOR_CFG['PULSAR_REFRESH_TIME'])
filterbank_directory = "/data/mopsr/archives/"
if FRB_DETECTOR_CFG['PLOTTING_THREAD'] == 'yes':
	plotting_thread_on = True
elif FRB_DETECTOR_CFG['PLOTTING_THREAD'] == 'no':
	plotting_thread_on = False

pulsar_db_file = FRB_DETECTOR_CFG['PULSAR_DB_FILE']
ns_threshold = float(FRB_DETECTOR_CFG['NS_THRESHOLD'])
md_threshold = float(FRB_DETECTOR_CFG['MD_THRESHOLD'])
sn_threshold = float(FRB_DETECTOR_CFG['SN_THRESHOLD'])
boxcar_threshold = float(FRB_DETECTOR_CFG['BOXCAR_THRESHOLD'])
dm_threshold = float(FRB_DETECTOR_CFG['DM_THRESHOLD'])
nsmd_port = int(FRB_DETECTOR_CFG['NSMD_PORT'])
nsmd_script = FRB_DETECTOR_CFG['NSMD_SCRIPT']


ok_tag = create_xml_elem('ok')
fail_tag = create_xml_elem('fail')


DD2R = np.pi/180
R2S = 43200/np.pi

molonglo = ephem.Observer()
molonglo.lat = ephem.degrees(-( 35 + 22/60.0 + 14.5452 / 3600.0) * DD2R)
molonglo.long = ephem.degrees((149 + 25/60.0 + 28.7682 / 3600.0) * DD2R)
#molonglo.elev = 400

time_cover_beam = (4*(351./352))/15*3600*(365.242190402/366.242190402)

eta = 2.37558703744e-5
zeta = -1/289.9
lmbda = -0.617335295908

cand_dtype={'names':('SN','sample','time','H_w','H_dm','beam')
		,'formats':('f4','i8','f4','i4','f4','i4')}


def main():
	# Parsing args
	# ------------
	parser = argparse.ArgumentParser(description='Master module that controls\
			real-live detection system. Handles start/stop signals, and\
			parses	HEIMDAL candidates to the corresponding BF nodes for\
			processing.')
	parser.add_argument('server',type=str)
	parser.add_argument('--nbf', type=int, help='Number of bf nodes\
			processing data',required=False)
	parser.add_argument('--verbose','-v',action="store_true",help='Verbose\
			mode.')
	parser.add_argument('--test','-t',action='store_true',help='Dry run.\
			Doesn\'t log to system log files.')
	parser.add_argument('--daemonize','-d',action='store_false',
			help='Don\'t daemonize', default = True)
	args = parser.parse_args()

	dry_run = args.test
	# --------------------------
	#n_bf = args.nbf
	n_bf = 8
	daemon = args.daemonize

	# Specifying level of verbosity
	# -----------------------------
	verbose = args.verbose

	# Spawning control thread and daemonizing
	# ---------------------------------------
	mopsr_cfg = parse_cfg(MOPSR_CFG_DIR,["SERVER_CONTROL_DIR",
		"SERVER_LOG_DIR","SERVER_HOST","FRB_DETECTOR_BASEPORT"])
	srv_host = mopsr_cfg["SERVER_HOST"]
	baseport = int(mopsr_cfg["FRB_DETECTOR_BASEPORT"])

	if dry_run:
		srv_ctrl_dir = FRB_DETECTOR_CFG['TEST_DIR']+'/control'
		srv_log_dir = FRB_DETECTOR_CFG['TEST_DIR']+'/logs'
	else:
		srv_ctrl_dir = mopsr_cfg['SERVER_CONTROL_DIR']
		srv_log_dir = mopsr_cfg['SERVER_LOG_DIR']
	
	pid = os.getpid()
	script_name = os.path.basename(sys.argv[0]).lstrip('server_').\
			rstrip('.py')
	logfile = srv_log_dir+'/'+script_name+'.log'
	pidfile = srv_ctrl_dir+'/'+script_name+'.pid'
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
	logging.info("SRV0 master script initializing")
	logging.info("Width boxcar threshold: %s, DM threshold: %s",
			boxcar_threshold,dm_threshold)
	if daemon:
		logging.info("Daemonizing")
		daemonize(pidfile, logfile)
	else:
		atexit.register(delpid,pidfile)
		pid = str(os.getpid())
		logging.debug("Writing pid file (pid %s)",pid)
		file(pidfile,'w+').write("%s\n" % pid)

	controlThread = threading.Thread(name = 'controlThread',
			target = control_monitor,
			args=(srv_ctrl_dir,script_name))
	controlThread.setDaemon(True)
	controlThread.start()


	# Parsing beam info from mopsr_bp_cornerturn.cfg file
	# ---------------------------------------------------
	logging.debug("Parsing bf configuration file")
	beam_config,bf_ips = parse_bf_config(CORNERTURN_CFG_DIR)
	logging.debug("Configuration files parsed in")

	# Saving IPs and port numbers for each BF node
	# --------------------------------------------
	bf_addrs = {}
	for i in range(n_bf):
		bf_addrs['RECV_'+str(i)] = (bf_ips['RECV_'+str(i)],baseport+100+(i+1))
	
	# Defining listening address
	# --------------------------
	listening_addrs = (srv_host,baseport)

	# Loading pulsar database
	# -----------------------
	if FRB_DETECTOR_CFG['PULSAR_MONITOR'] == 'yes':
		pulsar_monitor_on = True
		logging.info("Pulsar monitor on")
	elif FRB_DETECTOR_CFG['PULSAR_MONITOR'] == 'no':
		pulsar_monitor_on = False
		logging.info("Pulsar monitor off")
	pulsar_db = get_pulsars_db(pulsar_db_file)


	# Initializing pulsar monitor thread
	# ----------------------------------
	proc = Popen(args=[nsmd_script,str(nsmd_port)],shell = False)
	time.sleep(0.05)

	# -----------------------------------------------------------------
	# -----------------------------------------------------------------
	# Entering infinite loop (idle state), and waiting for start signal
	# -----------------------------------------------------------------
	# -----------------------------------------------------------------
	while True:
		logging.debug("Listening to backend for observation start")
		flag,obsInfo = socket_listen(listening_addrs)
		if flag != 'START': #If not a start flag, listen to backend again
			logging.critical('Received a non-starting flag: (%s)',flag)
			logging.critical('Trying again')
			continue

		start_utc = obsInfo['UTC_START']
		logging.info("Received a new start UTC: %s",start_utc)

		# Sending current utc to BF nodes
		# -------------------------------
		logging.debug("Sending UTC to bf_nodes")
		send_utc_to_bf(start_utc,obsInfo['SOURCE'],bf_addrs)
		
		logging.debug("Parsing observation parameters")
		SOURCE = obsInfo['SOURCE']
		boresight_ra , boresight_dec = obsInfo['RA'],obsInfo['DEC']
		mode = obsInfo['MODE']
		utc_start = obsInfo['UTC_START']
		observing_type = obsInfo['OBSERVING_TYPE'].upper()
		pulsar_list = []
		refined_pulsar_db = []
		utc_start_datetime = datetime.datetime.strptime(utc_start,"%Y-%m-%d-%H:%M:%S")
		if observing_type == "TRACKING" and pulsar_monitor_on:
			refined_pulsar_db = refine_pulsar_db(observing_type,pulsar_db,
					boresight_ra = obsInfo['RA'],boresight_dec = obsInfo['DEC'])
			pulsar_list = getPotentialPulsars_tracking(utc_start,boresight_ra,
					boresight_dec,beam_config["NBEAM"],refined_pulsar_db)
			time_tag = time.time()
			logging.debug("Tracking observation, pulsar monitor is on")
		elif observing_type == "STATIONARY" and pulsar_monitor_on:
			refined_pulsar_db = refine_pulsar_db(observing_type,pulsar_db,boresight_dec = obsInfo['DEC'])
			logging.debug("Stationary observation, pulsar monitor is on")

		# ---------------------------------------------------------------------
		# ---------------------------------------------------------------------
		# Observation commenced, waiting for Heimdals candidates/or stop signal
		# ---------------------------------------------------------------------
		# ---------------------------------------------------------------------
		logging.debug("Observation commenced")
		observing = True
		while observing:
			logging.debug("listening for Heimdals' candidates or STOP signal")
			flag, heimdal_candidates = socket_listen(listening_addrs)
			if flag == "CANDIDATES" and len(heimdal_candidates) > 0:
				#Obtained HEIMDAL Candidates
				logging.info("Received candidates from HEIMDAL")
				if observing_type == "TRACKING":
					t = time.time()
					utc_now = utc_start_datetime + datetime.timedelta(
							milliseconds = heimdal_candidates[0]['time']*1000)#NOTE: for testing
					utc_now = datetime.datetime.strftime(
							utc_now,"%Y-%m-%d-%H:%M:%S")	#NOTE: for testing
					if pulsar_monitor_on:
						pulsar_list = getPotentialPulsars_tracking(utc_now,boresight_ra,boresight_dec,
								beam_config["NBEAM"],refined_pulsar_db)
					filtered_candidates = candidateFilter(heimdal_candidates,pulsar_list)
					if filtered_candidates.size != 0:
						send_cands_to_bf(bf_addrs,beam_config,filtered_candidates)
				elif observing_type == "STATIONARY":
					threshold_candidates = thresholdFilter(heimdal_candidates)
					n_cands = threshold_candidates.size
					if n_cands == 0:
						continue
					threshold_candidates.sort(order = 'time')
					low_t = threshold_candidates[0]['time']
					ind = 1
					bulk = [0]
					logging.debug("number of candidates: %i",n_cands)
					while ind < n_cands:
						if threshold_candidates[ind]['time'] > 1+low_t: #NOTE: '1' is for 1 second of candidates per gulp
							if pulsar_monitor_on:
								molonglo.date = utc_start_datetime + datetime.timedelta(milliseconds = low_t*1000)
								pulsar_list = getPotentialPulsars_stationary(refined_pulsar_db)
							send_cands_to_bf(bf_addrs,beam_config,candidateFilter(threshold_candidates[bulk],pulsar_list,False))
							low_t = threshold_candidates[ind]['time']
							bulk=[]
						else:
							bulk.append(ind)
							ind+=1
					# Flushing last loop
					# ------------------
					logging.debug("Flushing last loop")
					molonglo.date = utc_start_datetime + datetime.timedelta(milliseconds = low_t*1000)
					pulsar_list = getPotentialPulsars_stationary(refined_pulsar_db)
					send_cands_to_bf(bf_addrs,beam_config,candidateFilter(threshold_candidates[bulk],pulsar_list,False))
				else:
					logging.critical("Observing type is neither 'STATIONARY' nor 'TRACKING'")
			elif flag == "STOP":
				#Obtained a Stop flag
				observing = False
			else:
				logging.critical("Unkown flag")
				logging.critical(heimdal_candidates)
		logging.info("Stop signal, back to idle mode")

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
			logging.exception("")
	except:
		logging.exception("")
