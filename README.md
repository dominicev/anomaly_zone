# anomaly_zone

Contact: Tyler K. Chafin (tkchafin@uark.edu)

This is my python3 interface for the anomaly_finder.py script written
by CW Linkem (https://github.com/cwlinkem/anomaly_zone)

If you use this, please cite the original publication:
Linkem, C. W., Minin, V. N., & Leache, A. D. (2016). Detecting the anomaly zone 
in species trees and evidence for a misleading signal in higher-level skink phylogeny 
(Squamata: Scincidae). Systematic Biology, 65(3), 465-477.

## Installation
There is only one prerequisite which won't be available with your standard Python 3 installation:
* [DendroPy](https://dendropy.org/) v. > 4.0.0

The easiest way to install this would be with conda:

```
conda install -c bioconda dendropy
```
or pip:
```
pip install dendropy
```

## Usage
View the help menu by calling anomaly_finder.py with the <-h> flag:

```
$ ./anomaly_finder.py -h

Usage:  ./anomaly_finder.py -t <treefile> -f <nexus or newick> 

Description: anomaly_finder.py by CW Linkem

	Arguments:
		-t	: Tree file 
		-f	: Format of tree file: "nexus" or "newick"
		-b	: Tree file containing bootstrap trees, if using them
		-h	: Displays help menu
```
anomaly_finder.py can be run in one of two ways: 1) With bootstrap trees provided, in which case results for anomalous edge pairs as the proportion of bootstrap trees in which the pair of edges falls under the AZ boundary; or 2) without bootstrap trees, in which edge pairs are simply reported as being under the AZ boundary (1) or not (0). 

## Inputs
The inputs for anomaly_finder.py are very simple. You need to provide a reference tree **with branch lengths scaled in coalescent units**, with an optional additional input of bootstrap trees (also caled in coalescent units)which enables you to account for topological error in branch length estimates.

This can be provided either in [NEWICK](http://evolution.genetics.washington.edu/phylip/newicktree.html), or [NEXUS](https://en.wikipedia.org/wiki/Nexus_file) format. 

## Outputs
anomaly_finder.py will produce output in the following format:
```
tyler:anomaly_zone $ python3 anomaly_finder.py -t test.tre 

####################anomaly_finder.py######################

Reading primary tree... test.tre
No bootstrap trees provided.
Finding anomalous nodes in primary tree...

######################ENCODING#############################

Writing results with splits encoded as:
                                                      /-------------------------- HBC gc          
/-----------------------------------------------------"111111111000"                              
|                                                     |            /------------- RTC U Green     
|                                                     \------------"000000000110"                 
|                                                                  \------------- RTC ECC         
"000000000000"                                                                                    
|-------------------------------------------------------------------------------- BTC             
|                                                                                                 
|            /------------------------------------------------------------------- seminuda virgin 
|            |                                                                                    
\------------"111111110000"/----------------------------------------------------- jordani         
             |             |                                                                      
             |             |                                       /------------- RTC L Verde WCC 
             \-------------"111111100000"--------------------------"000011000000"                 
                           |            |                          \------------- RTC L Gila CCK  
                           |            |                                                         
                           \------------"111111000000"             /------------- GLC uGila EAG   
                                        |             /------------"001100000000"                 
                                        |             |            \------------- GLC AguaFria SYC
                                        \-------------"111100000000"                              
                                                      |            /------------- HWC Salt GOR    
                                                      \------------"110000000000"                 
                                                                   \------------- HWC Verde WEB   
                                                                                                  
                                                                                                  
(If tree is large, import the following into FigTree)

[&U] ((HBC_gc:0.10000000000000009,(RTC_U_Green:0.10000000000000009,RTC_ECC:0.10000000000000009)'"000000000110"':0.3104037138)'"111111111000"':1.201628316,BTC:0.10000000000000009,(seminuda_virgin:0.10000000000000009,(jordani:0.10000000000000009,((RTC_L_Verde_WCC:0.10000000000000009,RTC_L_Gila_CCK:0.10000000000000009)'"000011000000"':0.0192994779500002,((GLC_uGila_EAG:0.10000000000000009,GLC_AguaFria_SYC:0.10000000000000009)'"001100000000"':0.0923875015300002,(HWC_Salt_GOR:0.10000000000000009,HWC_Verde_WEB:0.10000000000000009)'"110000000000"':0.0)'"111100000000"':0.057935799650000064)'"111111000000"':0.2587155266000001)'"111111100000"':0.7787775563000001)'"111111110000"':0.5662261648)'"000000000000"';
```
The above section displays an ascii and newick representation of how splits are labeled in the tree. Splits will take the form "000110000", corresponding to the leaf nodes which are descendents. See [the Dendropy documentation](https://dendropy.org/primer/bipartitions.html?highlight=bitmask) for more information. 

The important part is that the node labels represent the edge subtending the labelled node, and that these labels will be used to communicate the output for edge pairs found to be in the AZ. 

If the tree is too large, you can also paste the newick output into a text file and load it into FigTree.
```
######################RESULTS##############################

Showing anomalous nodes only:

ParentEdge	DescendantEdge
111100000000	001100000000
111100000000	110000000000
```
This table shows the outputs for the AZ tests. In this example, there are two edge pairs found to be in the AZ. If you had provided bootstrap trees, the output would look like this:
```
######################RESULTS##############################

Showing anomalous nodes only:

ParentEdge	DescendantEdge	Prop.Anomalous
111100000000	001100000000	1.0
111100000000	110000000000	1.0
```
The last part of the output shows the results mapped onto the tree as node labels. Here, we see once again the ascii and newick-formatted outputs, with AZ 'pairs' labelled according to a unique integer label, with '0' representing edges not contained within any anomalous divergences. 

Note that edges contributing to multiple AZ events will have identifiers separated by "/"
```
########################OUTPUT############################

Tree with anomalous edge pairs annotated (AZ pairs numbered)
                                                      /-------------------------- HBC gc          
/-----------------------------------------------------0                                           
|                                                     |            /------------- RTC U Green     
|                                                     \------------0                              
|                                                                  \------------- RTC ECC         
0                                                                                                 
|-------------------------------------------------------------------------------- BTC             
|                                                                                                 
|            /------------------------------------------------------------------- seminuda virgin 
|            |                                                                                    
\------------0             /----------------------------------------------------- jordani         
             |             |                                                                      
             |             |                                       /------------- RTC L Verde WCC 
             \-------------0            /--------------------------0                              
                           |            |                          \------------- RTC L Gila CCK  
                           |            |                                                         
                           \------------0                          /------------- GLC uGila EAG   
                                        |             /------------1                              
                                        |             |            \------------- GLC AguaFria SYC
                                        \-------------1/2                                         
                                                      |            /------------- HWC Salt GOR    
                                                      \------------2                              
                                                                   \------------- HWC Verde WEB   
                                                                                                  
                                                                                                  
(If tree is large, import the following into FigTree)

[&U] ((HBC_gc:0.10000000000000009,(RTC_U_Green:0.10000000000000009,RTC_ECC:0.10000000000000009)0:0.3104037138)0:1.201628316,BTC:0.10000000000000009,(seminuda_virgin:0.10000000000000009,(jordani:0.10000000000000009,((RTC_L_Verde_WCC:0.10000000000000009,RTC_L_Gila_CCK:0.10000000000000009)0:0.0192994779500002,((GLC_uGila_EAG:0.10000000000000009,GLC_AguaFria_SYC:0.10000000000000009)1:0.0923875015300002,(HWC_Salt_GOR:0.10000000000000009,HWC_Verde_WEB:0.10000000000000009)2:0.0)'1/2':0.057935799650000064)0:0.2587155266000001)0:0.7787775563000001)0:0.5662261648)0;


########################DONE!############################
```
If you had provided bootstrap trees, you would also get the associated proportions of BS trees found to be in the AZ:
```
########################OUTPUT############################

Tree with anomalous pairs annotated as "AZ Event : Prop. BS"
                                                      /-------------------------- HBC gc          
/-----------------------------------------------------0:0                                         
|                                                     |            /------------- RTC U Green     
|                                                     \------------0:0                            
|                                                                  \------------- RTC ECC         
0:0                                                                                               
|-------------------------------------------------------------------------------- BTC             
|                                                                                                 
|            /------------------------------------------------------------------- seminuda virgin 
|            |                                                                                    
\------------0:0           /----------------------------------------------------- jordani         
             |             |                                                                      
             |             |                                       /------------- RTC L Verde WCC 
             \-------------0:0          /--------------------------0:0                            
                           |            |                          \------------- RTC L Gila CCK  
                           |            |                                                         
                           \------------0:0                        /------------- GLC uGila EAG   
                                        |             /------------1:1.0                          
                                        |             |            \------------- GLC AguaFria SYC
                                        \-------------1/2:1.0/1.0                                 
                                                      |            /------------- HWC Salt GOR    
                                                      \------------2:1.0                          
                                                                   \------------- HWC Verde WEB   
                                                                                                  
                                                                                                  
(If tree is large, import the following into FigTree)

[&U] ((HBC_gc:0.10000000000000009,(RTC_U_Green:0.10000000000000009,RTC_ECC:0.10000000000000009)'0:0':0.3104037138)'0:0':1.201628316,BTC:0.10000000000000009,(seminuda_virgin:0.10000000000000009,(jordani:0.10000000000000009,((RTC_L_Verde_WCC:0.10000000000000009,RTC_L_Gila_CCK:0.10000000000000009)'0:0':0.0192994779500002,((GLC_uGila_EAG:0.10000000000000009,GLC_AguaFria_SYC:0.10000000000000009)'1:1.0':0.0923875015300002,(HWC_Salt_GOR:0.10000000000000009,HWC_Verde_WEB:0.10000000000000009)'2:1.0':0.0)'1/2:1.0/1.0':0.057935799650000064)'0:0':0.2587155266000001)'0:0':0.7787775563000001)'0:0':0.5662261648)'0:0';


########################DONE!############################
```
Once again, you can paste these outputs into FigTree if the tree is too large for viewing via stdout.

## Changelog
11 March 2019
- Ported to Python3 
- Ported DendroPy code to DendroPy 4.0.0 (e.g. converted deprecated function calls)
- Added basic argument parsing
- Made bootstrap tree input optional
- Created some graphical and text outputs (e.g. NEWICK)
