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

def main():
	params = parseArgs()

	taxon_set = dendropy.TaxonSet()
	trees = dendropy.TreeList.get_from_path('mpest_boots.tre', "nexus", taxon_set = taxon_set)
	mrc = dendropy.Tree.get_from_path('mpest_eMRC.tre', "nexus", taxon_set = taxon_set)

	master_dict=dict()
	mrc_dict=dict()
	split_occ_dict=dict()
	for sp_tree in trees:
		dendropy.treesplit.encode_splits(sp_tree)
		master_dict = get_nodes(sp_tree, master_dict)
		split_occ_dict=split_freq(sp_tree, taxon_set, split_occ_dict)	

	dendropy.treesplit.encode_splits(mrc)
	split_dict=split_mapper(mrc, taxon_set) #match a list of taxa to the split pattern from the bitmask for the whole tree
	mrc_dict=get_nodes(mrc, mrc_dict)


	for k,v in master_dict.iteritems():
		taxlab1=list()
		taxlab2=list()
		if k in mrc_dict.keys():
			tax1 = dendropy.treesplit.split_as_string(k[0], width = len(taxon_set))
			tax11 = taxon_set.split_taxa_list(k[0])
			tax2 = dendropy.treesplit.split_as_string(k[1], width = len(taxon_set))
			tax22 = taxon_set.split_taxa_list(k[1])
			for i in tax11:
				taxlab1.append(i.label)
			for j in tax22:
				taxlab2.append(j.label)
				
			print(tax1, tax2, len(v), sum(v))


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
			tax = taxon_set.split_taxa_list(node.edge.split_bitmask)
			for i in tax:
				taxlist.append(i.label)
			split_dict[dendropy.treesplit.split_as_string(node.edge.split_bitmask, width=len(taxon_set))] = taxlist
	return split_dict

def split_freq(tree,taxon_set, split_occ_dict):
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			split=dendropy.treesplit.split_as_string(node.edge.split_bitmask, width=len(taxon_set))
			if split not in split_occ_dict.keys():
				split_occ_dict[split] = [1]
			else:
				split_occ_dict[split].append(1)
	return split_occ_dict
			
def get_nodes(tree, master_dict):
	pair_list=[]
	for node in tree.postorder_node_iter():
		if node.parent_node is None:
			node.value = 1.0
		else:
			if node.taxon is None:
				node_pair = (node.parent_node.edge_length, node.edge_length)
				anomalous=anomaly_calc(node_pair)
				edge_pair = (node.parent_node.edge.split_bitmask, node.edge.split_bitmask)
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
			options, remainder = getopt.getopt(sys.argv[1:], 't:f:h', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.tree=None
		self.ftype=None


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
				self.ftype =arg
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
		-f		: Format of tree file: "nexus" or "newick"
		-h		: Displays help menu

""")
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()

		
