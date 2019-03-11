#! /usr/bin/python

"""
This is just a python3 interface for the anomaly_finder.py script written
by CW Linkem (https://github.com/cwlinkem/anomaly_zone)

If you use this, please cite the original publication:
Linkem, C. W., Minin, V. N., & Leache, A. D. (2016). Detecting the anomaly zone 
in species trees and evidence for a misleading signal in higher-level skink phylogeny 
(Squamata: Scincidae). Systematic Biology, 65(3), 465-477.

"""

import sys
import os
import math
from math import exp
import dendropy
import getopt
from dendropy.utility import bitprocessing

def main():
	params = parseArgs()

	#read trees
	taxon_set = dendropy.TaxonNamespace()
	if params.btrees:
		trees = dendropy.TreeList.get_from_path(params.btrees, params.ftype, taxon_namespace = taxon_set)
	mrc = dendropy.Tree.get_from_path(src=params.tree, schema=params.ftype, taxon_namespace = taxon_set)

	#initialize some dictionaries
	master_dict=dict()
	mrc_dict=dict()
	split_dict=dict()
	split_occ_dict=dict()

	
	#if bootstrap trees provided, get results for them
	#master_dict:
	#key = pair of edges
	#value = list of anomaly zone results (1=True; 0=False)
	if params.btrees:
		for sp_tree in trees:
			#dendropy.treesplit.encode_splits(sp_tree) #deprecated
			sp_tree.encode_bipartitions()
			#split_bitmasks = sp_tree.split_bitmask_edge_map.keys()
			master_dict = get_nodes(sp_tree, master_dict)
			split_occ_dict=split_freq(sp_tree, taxon_set, split_occ_dict)	

	#sys.exit(0)
	#get results for main tree
	mrc.encode_bipartitions()
	mrc_dict=get_nodes(mrc, mrc_dict) #get AZ results for primary tree
	split_dict=split_mapper(mrc, taxon_set)  #match a list of taxa to the split pattern from the bitmask for the whole tree


	#if bootstrap trees provided:
	#for each internal edge pair:
	if params.btrees:
		for k,v in master_dict.items():
			taxlab1=list()
			taxlab2=list()
			#print(k[0])
			#if edge pair has results for main tree
			if k in mrc_dict.keys():
				print("k:",k, "v:",v)
				#bitstring for parent node 
				parentLabel = bitprocessing.int_as_bitstring(k[0], length=len(taxon_set))
				#print(parentLabel)
				#list of descendent taxa
				parentList = taxon_set.bitmask_taxa_list(k[0])
				#print(parentList)
				#label for descendant edge
				descendantLabel = bitprocessing.int_as_bitstring(k[1], length=len(taxon_set))
				#list of descendant taxa
				descendantList = taxon_set.bitmask_taxa_list(k[1])
				#print(descendantList)
				for i in parentList:
					taxlab1.append(i.label)
				for j in descendantList:
					taxlab2.append(j.label)
				print(parentLabel, descendantLabel, "Prop. anomalous:",sum(v)/len(v))
				print()
	else:
		for k,v in mrc_dict.items():
			taxlab1=list()
			taxlab2=list()
			anom=False
			print("k:",k, "v:",v)
			#bitstring for parent node 
			parentLabel = bitprocessing.int_as_bitstring(k[0], length=len(taxon_set))
			#print(parentLabel)
			#list of descendent taxa
			parentList = taxon_set.bitmask_taxa_list(k[0])
			#print(parentList)
			#label for descendant edge
			descendantLabel = bitprocessing.int_as_bitstring(k[1], length=len(taxon_set))
			#list of descendant taxa
			descendantList = taxon_set.bitmask_taxa_list(k[1])
			#print(descendantList)
			for i in parentList:
				taxlab1.append(i.label)
			for j in descendantList:
				taxlab2.append(j.label)
			print(parentLabel, descendantLabel, "Anomalous:",sum(v)/len(v))
			print()

#For a pair of edges, return 1 if anomalous; return 0 if not anomalous
def anomaly_calc(nodes):
	lambda1=nodes[0]
	lambda2=nodes[1]
	if lambda1 is not None:
		Z1=math.log(2.0/3+((3*exp(2*lambda1)-2)/(18*(exp(3*lambda1)-exp(2*lambda1))))) #function a(x) from Degnan and Rosenberg (2006)
		if lambda2 <= Z1:
			return 1
		else:
			return 0
	else:
		return 0

def split_mapper(tree, taxon_set):
	split_dict={}
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			taxlist=[]
			tax = taxon_set.bitmask_taxa_list(node.edge.split_bitmask)
			for i in tax:
				taxlist.append(i.label)
			split_dict[bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=len(taxon_set))] = taxlist
	return split_dict

#function
def split_freq(tree,taxon_set, split_occ_dict):
	#for each node 
	for node in tree.postorder_node_iter():
		#if not root
		if node.parent_node is None:
			node.value = 1.0
		else:
			#split=dendropy.treesplit.split_as_string(node.edge.split_bitmask, width=len(taxon_set)) #deprecated
			split=bitprocessing.int_as_bitstring(node.edge.split_bitmask, length=len(taxon_set))
			print(split)
			if split not in split_occ_dict.keys():
				split_occ_dict[split] = [1]
			else:
				split_occ_dict[split].append(1)
	return split_occ_dict
	
#Function tests for internal edge pairs under anomaly zone
#modifies a dict where key = pair of edges; calue = list of results (1=true; 0=false)		
def get_nodes(tree, master_dict):
	pair_list=[]
	#for each internal node
	#calculate
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			#root
			node.value = 1.0
		else:
			#if not root, and internal node
			if node.taxon is None:
				#get edge length of parent, and own edge length
				#x = parent; y = current
				node_pair = (node.parent_node.edge_length, node.edge_length)
				#calculate if pair of edges are under anomaly boundary
				anomalous=anomaly_calc(node_pair)
				#print(anomalous)
				edge_pair = (node.parent_node.edge.split_bitmask, node.edge.split_bitmask)
				#for edge pair, keep anomaly zone tests results
				if edge_pair not in master_dict.keys():
					master_dict[edge_pair] = [anomalous]
				else:
					master_dict[edge_pair].append(anomalous)	
	return master_dict		

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 't:f:b:h', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.tree=None
		self.ftype="newick"
		self.btrees=None


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg_raw in options:
			arg = arg_raw.replace(" ","")
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == "t":
				self.tree = arg
			elif opt in ('h', 'help'):
				pass
			elif opt == "f":
				if arg == "newick" or arg == "nexus":
					self.ftype =arg
				else:
					self.display_help("Invalid option for -f",arg)
			elif opt == "b":
				self.btrees=arg
			else:
				assert False, "Unhandled option %r"%opt

		#Check manditory options are set
		if not self.tree:
			self.display_help("Error: Need treefile <-t>")
		if not self.ftype:
			print("No filetype (-f) provided: using default (\"newick\")")



	def display_help(self, message=None):
		if message is not None:
			print ("\n",message)
		print ("\nanomaly_finder.py\n")
		print ("\nUsage: ", sys.argv[0], "-t <treefile> -f <nexus or newick> \n")
		print ("Description: anomaly_finder.py by CW Linkem")

		print("""
	Arguments:
		-t		: Tree file (branches scaled by coalescent units)
		-b		: Tree file containing bootstrap trees, if being used
		-f		: Format of tree file: "nexus" or "newick"
		-h		: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()

		
