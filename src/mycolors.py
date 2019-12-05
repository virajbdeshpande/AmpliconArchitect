# This software is Copyright 2017 The Regents of the University of California. All Rights Reserved. Permission to copy, modify, and distribute this software and its documentation for educational, research and non-profit purposes, without fee, and without a written agreement is hereby granted, provided that the above copyright notice, this paragraph and the following three paragraphs appear in all copies. Permission to make commercial use of this software may be obtained by contacting:
#
# Office of Innovation and Commercialization
#
# University of California
#
# La Jolla, CA 92093-0910
#
# (858) 534-5815
#
# invent@ucsd.edu
#
# This software program and documentation are copyrighted by The Regents of the University of California. The software program and documentation are supplied "as is", without any accompanying services from The Regents. The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The end-user understands that the program was developed for research purposes and is advised not to rely exclusively on the program for any reason.
#
# IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

#Author: Viraj Deshpande
#Contact: virajbdeshpande@gmail.com



chrcolor = {
'b'	: 'b',
'g'	: 'g',
'r'	: 'r',
'c'	: 'c',
'm'	: 'm',
'y'	: 'y',
'k'	: 'k',
'w'	: 'w',
'chr1' : (153/256.0, 102/256.0, 0/256.0),
'chr2' : (102/256.0, 102/256.0, 0/256.0),
'chr3' : (153/256.0, 153/256.0, 30/256.0),
'chr4' : (204/256.0, 0/256.0, 0/256.0),
'chr5' : (255/256.0, 0/256.0, 0/256.0),
'chr6' : (255/256.0, 0/256.0, 204/256.0),
'chr7' : (255/256.0, 204/256.0, 204/256.0),
'chr8' : (255/256.0, 153/256.0, 0/256.0),
'chr9' : (255/256.0, 204/256.0, 0/256.0),
'chr10': (255/256.0, 255/256.0, 0/256.0),
'chr11': (204/256.0, 255/256.0, 0/256.0),
'chr12': (0/256.0, 255/256.0, 0/256.0),
'chr13': (53/256.0, 128/256.0, 0/256.0),
'chr14': (0/256.0, 0/256.0, 204/256.0),
'chr15': (102/256.0, 153/256.0, 255/256.0),
'chr16': (153/256.0, 204/256.0, 255/256.0),
'chr17': (0/256.0, 255/256.0, 255/256.0),
'chr18': (204/256.0, 255/256.0, 255/256.0),
'chr19': (153/256.0, 0/256.0, 204/256.0),
'chr20': (204/256.0, 51/256.0, 255/256.0),
'chr21': (204/256.0, 153/256.0, 255/256.0),
'chr22': (102/256.0, 102/256.0, 102/256.0),
'chr23': (153/256.0, 153/256.0, 153/256.0),
'chrX' : (153/256.0, 153/256.0, 153/256.0),
'chr24': (204/256.0, 204/256.0, 204/256.0),
'chrY' : (204/256.0, 204/256.0, 204/256.0),
'chrM' : (204/256.0, 204/256.0, 153/256.0),
'chr0' : (204/256.0, 204/256.0, 153/256.0),
'chrUn': (121/256.0, 204/256.0, 61/256.0),
'chrNA': (255/256.0, 255/256.0, 255/256.0),

'lum90chr1'  : (255/256.0,216/256.0,156/256.0),
'lum90chr2'  : (230/256.0,230/256.0,165/256.0),
'lum90chr3'  : (232/256.0,232/256.0,135/256.0),
'lum90chr4'  : (255/256.0,166/256.0,166/256.0),
'lum90chr5'  : (255/256.0,147/256.0,147/256.0),
'lum90chr6'  : (255/256.0,152/256.0,255/256.0),
'lum90chr7'  : (255/256.0,214/256.0,214/256.0),
'lum90chr8'  : (255/256.0,202/256.0,102/256.0),
'lum90chr9'  : (255/256.0,220/256.0,58/256.0),
'lum90chr10' : (234/256.0,234/256.0,0/256.0),
'lum90chr11' : (194/256.0,245/256.0,0/256.0),
'lum90chr12' : (34/256.0,255/256.0,34/256.0),
'lum90chr13' : (174/256.0,244/256.0,155/256.0),
'lum90chr14' : (215/256.0,215/256.0,255/256.0),
'lum90chr15' : (182/256.0,224/256.0,255/256.0),
'lum90chr16' : (182/256.0,231/256.0,255/256.0),
'lum90chr17' : (0/256.0,252/256.0,252/256.0),
'lum90chr18' : (185/256.0,236/256.0,236/256.0),
'lum90chr19' : (255/256.0,191/256.0,255/256.0),
'lum90chr20' : (255/256.0,177/256.0,255/256.0),
'lum90chr21' : (255/256.0,206/256.0,255/256.0),
'lum90chr22' : (198/256.0,198/256.0,198/256.0),
'lum90chr23' : (153/256.0,153/256.0,153/256.0),
'lum90chrX'  : (153/256.0,153/256.0,153/256.0),
'lum90chr24' : (204/256.0,204/256.0,204/256.0),
'lum90chrY'  : (204/256.0,204/256.0,204/256.0),
'lum90chrM'  : (174/256.0,174/256.0,122/256.0),
'lum90chr0'  : (174/256.0,174/256.0,122/256.0),
'lum90chrUn' : (108/256.0,191/256.0,38/256.0),
'lum90chrNA' : (171/256.0,171/256.0,171/256.0),

'lum80chr1'  : (244/256.0,188/256.0,127/256.0),
'lum80chr2'  : (202/256.0,202/256.0,136/256.0),
'lum80chr3'  : (203/256.0,203/256.0,103/256.0),
'lum80chr4'  : (255/256.0,137/256.0,137/256.0),
'lum80chr5'  : (255/256.0,116/256.0,116/256.0),
'lum80chr6'  : (255/256.0,119/256.0,255/256.0),
'lum80chr7'  : (237/256.0,186/256.0,186/256.0),
'lum80chr8'  : (255/256.0,174/256.0,62/256.0),
'lum80chr9'  : (243/256.0,192/256.0,0/256.0),
'lum80chr10' : (206/256.0,206/256.0,0/256.0),
'lum80chr11' : (166/256.0,216/256.0,0/256.0),
'lum80chr12' : (0/256.0,232/256.0,0/256.0),
'lum80chr13' : (146/256.0,216/256.0,126/256.0),
'lum80chr14' : (186/256.0,186/256.0,255/256.0),
'lum80chr15' : (152/256.0,196/256.0,255/256.0),
'lum80chr16' : (152/256.0,203/256.0,254/256.0),
'lum80chr17' : (0/256.0,224/256.0,224/256.0),
'lum80chr18' : (156/256.0,208/256.0,208/256.0),
'lum80chr19' : (250/256.0,161/256.0,255/256.0),
'lum80chr20' : (255/256.0,146/256.0,255/256.0),
'lum80chr21' : (227/256.0,177/256.0,255/256.0),
'lum80chr22' : (198/256.0,198/256.0,198/256.0),
'lum80chr23' : (153/256.0,153/256.0,153/256.0),
'lum80chrX'  : (153/256.0,153/256.0,153/256.0),
'lum80chr24' : (204/256.0,204/256.0,204/256.0),
'lum80chrY'  : (204/256.0,204/256.0,204/256.0),
'lum80chrM'  : (174/256.0,174/256.0,122/256.0),
'lum80chr0'  : (174/256.0,174/256.0,122/256.0),
'lum80chrUn' : (108/256.0,191/256.0,38/256.0),
'lum80chrNA' : (171/256.0,171/256.0,171/256.0),

'vlpurple'   : (218/256.0,218/256.0,235/256.0),
'vlorange'   : (253/256.0,208/256.0,162/256.0),
'vlpgreen'   : (218/256.0,218/256.0,235/256.0)

}





ecolor = {'interchromosomal' : 'blue',
          'concordant' : 'black',
          'everted' : (139/256.0, 69/256.0, 19/256.0), # 'brown', yellow',
          'forward' : 'magenta',
          'reverse' : (0/256.0, 139/256.0, 139/256.0), #'cyan',
          'discordant' : 'red'}

