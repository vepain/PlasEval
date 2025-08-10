import networkx as nx
import itertools
from collections import defaultdict
import copy
from math import factorial

import logging
import psutil
import time
import sys

logger = logging.getLogger(__name__)

def generate_matchings(m, n):
	'''
	Input: Number of copies of a contig in left and right plasmid sets
	Returns: List of matchings where each matching is a pair of lists of indices (int) of the contig copies, one for each side.
	'''
	matchings = []
	if m >= n:
		pmutns = list(itertools.permutations(list(range(m)), n))
		n_list = list(range(n))
		for pmutn in pmutns:
			matchings.append((list(pmutn), n_list))
	else:
		pmutns = list(itertools.permutations(list(range(n)), m))
		m_list = list(range(m))
		for pmutn in pmutns:
			matchings.append((m_list, list(pmutn)))
	return matchings

def get_matching_positions(ctg_copies, matching):
	'''
	Input:
		Dictionary of contig copies: L_copies: list of contig copies in left plasmid
									 R_copies: list of contig copies in right plasmid
									 Each copy is a triple [contig, plasmid index (int), position in plasmid (int)]
		Matching: Pair of lists of indices (int) of the contig copies, one for each side
	Returns: Pair of lists of contig copies, one for each side, according to respective indices in the matching 
	'''	
	L, R = ctg_copies['L_copies'], ctg_copies['R_copies']
	l_posn, r_posn = [], []
	for x in matching[0]:
		l_posn.append(L[x])
	for x in matching[1]:
		r_posn.append(R[x])	
	return l_posn, r_posn

def rename_by_matching(matching_dict):
	'''
	Input: 
		Dictionary of matchings: Key: contig, 
								 Value: matching for copies of contig
								 Each matching is pair of lists of indices (int) of the contig copies, one for each side
	Returns:
	 	Two lists of contig copies, one for each plasmid set, renamed with indices according to the input matching 
		Each contig is a triple [contig, plasmid index (int), position in plasmid (int)]
	'''
	reached_contigs = matching_dict.keys()
	left_copies_renamed, right_copies_renamed = [], []
	for contig in reached_contigs:
		M = matching_dict[contig]
		left_ctgs = M[0]
		right_ctgs = M[1]
		for i in range(len(M[0])):
			lpls, lidx = left_ctgs[i][1], left_ctgs[i][2]
			left_copies_renamed.append([contig+'_'+str(i), lpls, lidx])
			rpls, ridx = right_ctgs[i][1], right_ctgs[i][2]
			right_copies_renamed.append([contig+'_'+str(i), rpls, ridx])
	return left_copies_renamed, right_copies_renamed

def add_nodes(G, left_ctg_list, right_ctg_list, pls_ids_dict):
	'''
	Input:
		Empty graph object G
		List of contig copies in left plasmid set
		List of contig copies in right plasmid set
		Dictionary of plasmids: Keys: L, R, Values: Bidict of plasmid indices <-> names/ids
	Returns:
		Bipartite graph G (no edges added yet)
	'''
	G.add_nodes_from(list(set([pls_ids_dict['L'].inv[x[1]] for x in left_ctg_list])), bipartite=0)
	G.add_nodes_from(list(set([pls_ids_dict['R'].inv[x[1]] for x in right_ctg_list])), bipartite=1)
	return G

def add_edges(G, left_ctg_list, right_ctg_list, pls_ids_dict):
	'''
	Input:
		Graph object G with vertices added
		List of contig copies in left plasmid set
		List of contig copies in right plasmid set
		Dictionary of plasmids: Keys: L, R, Values: Bidict of plasmid indices <-> names/ids
	Returns:
		Bipartite graph G (with edges added)
	'''	
	def list_ctg_ids_by_pls(pls_ids_dict_one_side, ctg_ids_list):
		'''
		Input:
			Bidict of plasmid indices <-> names/ids,
			List of contig copies
		Returns: 
			Dictionary: Key: plasmid id, 
						Value: List of contig ids in the plasmid from given contig list
		'''
		ctg_ids_by_pls = defaultdict(list)
		for ctg in ctg_ids_list:
			ctg_id, pls_id = ctg[0], pls_ids_dict_one_side.inv[ctg[1]]
			ctg_ids_by_pls[pls_id].append(ctg_id)
		return ctg_ids_by_pls
	left_ctg_ids_by_pls = list_ctg_ids_by_pls(pls_ids_dict['L'], left_ctg_list)
	right_ctg_ids_by_pls = list_ctg_ids_by_pls(pls_ids_dict['R'], right_ctg_list)
	
	edges_list = []
	for left_plasmid in left_ctg_ids_by_pls:
		for right_plasmid in right_ctg_ids_by_pls:
			left_keys = set(left_ctg_ids_by_pls[left_plasmid])
			right_keys = set(right_ctg_ids_by_pls[right_plasmid])
			if left_keys.intersection(right_keys) != set():
				edges_list.append((left_plasmid, right_plasmid))
	G.add_edges_from(edges_list)
	return G

def modify_partitions(partitions, common):
	'''
	Input:
		partitions: List of sets of contigs
		common: Set of contigs common to both plasmids 	
	Returns:
		list of partitions, modified by splitting each partition if required, according to the set of common contigs
	'''
	modified_partitions = []
	for S in partitions:
		if len(S.intersection(common)) != 0 or S.intersection(common) != S:
			modified_partitions.append(S.intersection(common))
			modified_partitions.append(S.difference(common))
		elif S.intersection(common) == S:
			modified_partitions.append(S)				
	return modified_partitions		

def get_partition_cost(partitions, contigs_dict, p):
	'''
	Input:
		partitions: List of sets of contigs
		contigs_dict: Key: contig (str), Value: Nested dictionary: 	length (int), 
																	L_copies/R_copies (list of contig copies in plasmid set)
	Returns:
		Total of (length of contig sets)^p and the cost of partitioning
	'''	
	total_len = 0
	largest_part_cost = 0
	for S in partitions:
		S_len = 0
		for contig in S:
			contig_len = contigs_dict[contig.split('_')[0]]['length']
			S_len += contig_len
		S_cost = S_len**p
		total_len += S_cost
		largest_part_cost = max(largest_part_cost, S_cost)
	cost = total_len - largest_part_cost
	return total_len, cost

def compute_splits_cost(pls_ids, side_contig_copies, opp_contig_copies, B, flag, pls_ids_dict, contigs_dict, p):
	'''
	Input:
		pls_ids: List of plasmid ids,
		side_contig_copies: list of contig copies in plasmid set
		opp_contig_copies: list of contig copies in opposite plasmid set
		B: bipartite graph object
		flag (binary): variable to indicate if side is left or right
		Dictionaries of plasmids and contigs
	Returns:
		Total cost of splits (cuts OR joins) for one plasmid set
	'''
	def get_ctg_list_by_pls(contig_copies, pls_ids_dict, side):
		'''
		Returns list of contigs for each plasmid
		'''
		ctgs_by_pls = {}
		for x in contig_copies:
			ctg, pls = x[0], pls_ids_dict[side].inv[x[1]]
			if pls not in ctgs_by_pls:
				ctgs_by_pls[pls] = []
			ctgs_by_pls[pls].append(ctg)
		return ctgs_by_pls
	[s, o] = ['L', 'R'] if flag == 0 else ['R', 'L']
	side_ctgs_by_pls = get_ctg_list_by_pls(side_contig_copies, pls_ids_dict, s)
	opp_ctgs_by_pls = get_ctg_list_by_pls(opp_contig_copies, pls_ids_dict, o)

	side_len, side_cost = 0, 0
	for node in pls_ids:
		partitions = [set(side_ctgs_by_pls[node])]
		for edge in B.edges:
			if edge[flag] == node:
				side_contigs, opp_contigs = set(side_ctgs_by_pls[edge[flag]]), set(opp_ctgs_by_pls[edge[1-flag]])
				common = side_contigs.intersection(opp_contigs)
				partitions = modify_partitions(partitions, common)
		node_len, cost = get_partition_cost(partitions,contigs_dict,p)
		side_len += node_len
		side_cost += cost
	return side_cost

def compute_match_cost(left_contig_copies, right_contig_copies, pls_ids_dict, contigs_dict, p):
	'''
	Input:
		List of contig copies, one for each side
		Dictionaries of plasmids and contigs
	Returns:
		Cost of cuts (left side splits) and joins (right side splits)
	'''
	B = nx.Graph()	#Create graph with vertices named according to matching and obtain connected components
	B = add_nodes(B, left_contig_copies, right_contig_copies, pls_ids_dict)
	B = add_edges(B, left_contig_copies, right_contig_copies, pls_ids_dict)
	A = [B.subgraph(c) for c in nx.connected_components(B)]
	left_pls_set = set([pls_ids_dict['L'].inv[x[1]] for x in left_contig_copies])
	n_conn_comp = len(list(A))
	left_splits_cost,right_splits_cost = 0, 0
	for i in range(n_conn_comp):
		C = list(A)[i]
		left_pls_ids, right_pls_ids = nx.bipartite.sets(C)							#Split the component according to bipartite sets
		if len(list(right_pls_ids)) != 0 and list(right_pls_ids)[0] in left_pls_set: 	#Ensuring proper assignments of bipartite parts
			right_pls_ids,left_pls_ids = left_pls_ids,right_pls_ids
		left_splits_cost += compute_splits_cost(left_pls_ids, left_contig_copies, right_contig_copies, B, 0, pls_ids_dict, contigs_dict, p)
		right_splits_cost += compute_splits_cost(right_pls_ids, right_contig_copies, left_contig_copies, B, 1, pls_ids_dict, contigs_dict, p)
	return left_splits_cost, right_splits_cost	

def compute_current_cost(matching_dict, pls_ids_dict, contigs_dict, p):
	'''
	Input:
		Dictionary of matchings: Key: contig, Value: matching for copies of contig,
		Dictionary of plasmids: Keys: side (L/R), Values: Bidict of plasmid indices <-> names/ids
		Dictionary of contigs: Key: contig (str), Value: Nested dictionary: length (int), 
																			L_copies/R_copies (list of contig copies in plasmid set)
	Returns:
		Cost of current matching
	'''	
	left_contig_copies, right_contig_copies = rename_by_matching(matching_dict)
	left_ctg_ids, right_ctg_ids = set(), set()
	for x in left_contig_copies:
		left_ctg_ids.add(x[0])
	for x in right_contig_copies:
		right_ctg_ids.add(x[0])
	return compute_match_cost(left_contig_copies, right_contig_copies, pls_ids_dict, contigs_dict, p)

def run_compare_plasmids(contigs_dict, pls_ids_dict, p, max_calls, results_file):
	'''
	Input:
		Dictionary of contigs: 
			Key: contig (str), Value: Nested dictionary:length (int), 
														L_copies: list of contig copies in left plasmid set
														R_copies: list of contig copies in right plasmid set
														Each copy is a triple [contig, plasmid index (int), position in plasmid (int)]
		Dictionary of plasmids, 
			Keys: L, R, Values: Bidict of plasmid indices <-> names/ids
	Returns:
		Dissimilarity score and associated costs (cuts, joins, contig copies present on only left or right plasmid sets)
	'''
	#Computing set of common contigs 
	left_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['L_copies']) >= 1])
	right_ctg_ids = set([ctg for ctg in contigs_dict.keys() if len(contigs_dict[ctg]['R_copies']) >= 1])
	common_contigs = left_ctg_ids.intersection(right_ctg_ids)

	#Computing upperbound on number of matchings and final_cost
	max_cost = 0
	n_matchings = {}
	max_n_matchings = 1
	for contig in common_contigs:
		m = len(contigs_dict[contig]['L_copies'])
		n = len(contigs_dict[contig]['R_copies'])
		max_cost += m * contigs_dict[contig]['length']
		max_cost += n * contigs_dict[contig]['length']	
		n_matchings[contig] = int(factorial(n)/factorial(n-m)) if n > m else int(factorial(m)/factorial(m-n))
		max_n_matchings *= n_matchings[contig]
	logger.info(f'Maximum possible matchings: {max_n_matchings}')

	start_time = time.time()
	#if max_n_matchings <= 10000000:
	dummy_var = 1
	if dummy_var == 1:
		### Branch-N-Bound ###
		current_state = {'level': 0, 'total_cost': 0, 'matching': {}, 'cuts_cost': 0, 'joins_cost': 0}
		final_state = {'total_cost': max_cost, 'matching': {}, 'cuts_cost': 0, 'joins_cost': 0}

		contig_list = list(common_contigs)
		sorted_contig_list = sorted(contig_list, key=lambda ctg: n_matchings[ctg])

		count = [0]

		def recursive_compare(current_state, sorted_contig_list, pls_ids_dict, contigs_dict, count):
			'''
			Input:
				Current state dictionary: 
					level: Distance from root of tree (int)
					total_cost: Cost of cuts and joins upto this level (int)
					matching: Nested dictionary withs contig ids (str) as keys and a pair (set) of lists of contigs as values
					cuts_cost, joins_cost: Cost of cuts, joins (respectively) upto this level (int)				
			Updates:
				Current state dictionary
				Final state dictionary (non local variable)
			'''
			nonlocal final_state
			#count[0] += 1
			if current_state['level'] < len(sorted_contig_list):				#Compute cost upto current level
				current_contig = sorted_contig_list[current_state['level']]		#Retrieve contig for current level				
				m = len(contigs_dict[current_contig]['L_copies'])
				n = len(contigs_dict[current_contig]['R_copies'])			
				matchings = generate_matchings(m,n); logger.info(f'Number of matchings: {len(matchings)}')
				for matching in matchings:
					matched_posns = get_matching_positions(contigs_dict[current_contig], matching)
					current_state['matching'][current_contig] = matched_posns
					count[0] += 1
					#if count[0] % 10000 == 0:
					#	print(count[0])
					if count[0] > max_calls:
						logger.info(f'Max number of iterations reached: {max_calls}'); sys.exit(f'Max number of iterations reached: {max_calls}')
					current_state['cuts_cost'], current_state['joins_cost'] \
						= compute_current_cost(current_state['matching'], pls_ids_dict, contigs_dict, p)
					current_state['total_cost'] = current_state['cuts_cost'] + current_state['joins_cost']
					if current_state['total_cost'] < final_state['total_cost']:	
						current_state['level'] += 1 
						recursive_compare(current_state, sorted_contig_list, pls_ids_dict, contigs_dict, count)
						current_state['level'] -= 1
					del current_state['matching'][current_contig]

			else:
				final_state['total_cost'] = current_state['total_cost']
				final_state['cuts_cost'], final_state['joins_cost'] = current_state['cuts_cost'], current_state['joins_cost']
				final_state['matching'] = copy.deepcopy(current_state['matching'])
		recursive_compare(current_state, sorted_contig_list, pls_ids_dict, contigs_dict, count)
		
		end_time = time.time()
		logger.info(f'Time taken: {end_time - start_time}')
		logger.info(f'Number of function calls: {count[0]}')	
		
		total_len, total_denom, unique_left_cost, unique_right_cost = 0, 0, 0, 0
		for c in contigs_dict:
			l_copies, r_copies = len(contigs_dict[c]['L_copies']), len(contigs_dict[c]['R_copies'])
			ctg_len = contigs_dict[c]['length']
			if min(l_copies, r_copies) == 0: #NOTE: This is 0 in order to not penalize extra copies of common contigs
				unique_left_cost += max(l_copies - r_copies, 0) * (ctg_len**p)
				unique_right_cost += max(r_copies - l_copies, 0) * (ctg_len**p)
			total_len += (l_copies + r_copies) * ctg_len
			total_denom += (l_copies + r_copies) * (ctg_len**p)

		dissimilarity_score = (unique_left_cost + unique_right_cost + final_state['total_cost'])
		#print("Total_ctg_length\t", total_len)
		#print("Total_ctg_length_alpha\t", total_denom)
		#print("Cuts_cost\t", final_state['cuts_cost'])
		#print("Joins_cost\t", final_state['joins_cost'])
		#print("Unique_left_ctgs\t", unique_left_cost)
		#print("Unique_right_ctgs\t", unique_right_cost)
		#print("Dissimilarity_score\t", dissimilarity)

		logger.info(f'{final_state["matching"]}')

		results_file.write("Total_ctg_length\t" + str(total_len) + "\n")
		results_file.write("Total_ctg_length_alpha\t" + str(total_denom) + "\n")
		results_file.write("Cuts\t" + str(final_state['cuts_cost']) + "\t" + str(final_state['cuts_cost']/total_denom) + "\n")
		results_file.write("Joins\t" + str(final_state['joins_cost']) + "\t" + str(final_state['joins_cost']/total_denom) + "\n")
		results_file.write("Extra_ctgs\t" + str(unique_left_cost) + "\t" + str(unique_left_cost/total_denom) + "\n")
		results_file.write("Missing_ctgs\t" + str(unique_right_cost) + "\t" + str(unique_right_cost/total_denom) + "\n")
		results_file.write("Dissimilarity\t" + str(dissimilarity_score) + "\t" + str(dissimilarity_score/total_denom) + "\n")
