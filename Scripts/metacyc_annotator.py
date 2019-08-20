#!/usr/bin/python

#Creator: Arnau Vich
#Year: 2016

#Usage: python metacyc_annotator.py ./IBD_WGS/Pathways_summary_2.txt test_annotation.txt
#where Pathways_summary first column contains metacyc pathways ids e.g PWY-3581

#Description: Script to retrieve information from Metacyc webpage. 

from lxml import html
import requests
import sys
import re
#Tipical input-output
input_file= open (sys.argv[1],'r')
output_file= open (sys.argv[2], 'w')
count_line=1
with open(sys.argv[1]) as input_file:
	for i, line in enumerate(input_file):
		#Print header at it is adding extra columns
		if count_line==1:
			count_line=2
			header= line.rstrip('\n')
			header = line.split("\t")
			header_1 = "Taxa with the pathway"
			header_2 = "NCBI IDs"
			header_3 = "Taxa range"
			header_4 = "NCBI IDs"
			for hd in header:
				output_file.write ('%s' '\t'% (hd))
			output_file.write ('%s' '\t' '%s' '\t' '%s' '\t' '%s' % (header_1, header_2, header_3, header_4))
			output_file.write ('\n')
		else: 
			line=line.rstrip('\n') #Here I remove the new line at the end of each line
			row= line.split("\t")
			inter_id= row [0]
			inter_id_2= inter_id.split(":")
			path_id= inter_id_2 [0] #Get the pathway id from the first column of each row
			other_info=inter_id_2[1:]
			url='http://websvc.biocyc.org/getxml?META:%s' %(path_id) #create the url to look at MetaCyc database
			page = requests.get(url) #Obtain the information in html
			tree = html.fromstring(page.content)
			ncbi_id = tree.xpath('//species//dblink-oid/text()') #parse the xml information, I don't know how to do it better
			common_name= tree.xpath('//species//common-name/text()')
			tax_range_id=tree.xpath('//taxonomic-range//dblink-oid/text()')
			tax_range_common_name=tree.xpath('//taxonomic-range//common-name/text()')
			output_file.write ('%s' % (path_id))
			for i in other_info:
				output_file.write ('%s' '\t'% (i))
			for a in common_name:
				output_file.write ('%s' ';'% (a))
			output_file.write ('\t')
			for x in ncbi_id:
				output_file.write ('%s' ';'% (x))
			output_file.write ('\t')
			for y in tax_range_common_name:
				output_file.write ('%s' ';'% (y))
			output_file.write ('\t')
			for z in tax_range_id:
				output_file.write ('%s' ';'% (z))
			output_file.write ('\n')
