from Bio import SeqIO
import xml.etree.ElementTree as ET
import copy
from operator import itemgetter
Filename = "EColi_vs_Salmonella_BLAST_results.xml"
filehandler = open(Filename)
#Query= NC_003197
QueryName = "Salmonella"
#hit = NC_000913
HitName = "EColi"
XML_Tree = ET.parse(filehandler)
the_root =  XML_Tree.getroot()
Array_of_hit_start = []
Array_of_hit_end = []
Array_of_query_start = []
Array_of_query_end = []

print the_root
for layer_one in the_root:
    print (layer_one," 1")
    for layer_two in layer_one:
        print (layer_two, " 2")
        for layer_three in layer_two:
            print (layer_three, " 3")
            for layer_four in layer_three:
                print (layer_four, " 4")
                for layer_five in layer_four:
                    print (layer_five, " 5")
                    for layer_six in layer_five:
                        print (layer_six, " 6")
                        for layer_seven in layer_six:
                            print (layer_seven.tag, " 7") 
                            if (layer_seven.tag == 'Hsp_hit-to'):
                                print ("this is the information I give a care about")
                                print (layer_seven.text)
                                Array_of_hit_end.append(layer_seven.text)
                            if (layer_seven.tag == 'Hsp_hit-from'):
                                print ("this is the information I give a care about")
                                print (layer_seven.text)
                                Array_of_hit_start.append(layer_seven.text)
                                
                            if (layer_seven.tag == 'Hsp_query-to'):
                                print ("this is the information I give a care about")
                                Array_of_query_end.append(layer_seven.text)
                            if (layer_seven.tag == 'Hsp_query-from'):
                                print ("this is the information I give a care about")
                                Array_of_query_start.append(layer_seven.text)

# time to pack elements into an array
StorageArray = [[0 for i in range(4)] for i in range(len(Array_of_hit_start))]
j=0
while j < len(Array_of_hit_start):
    StorageArray[j][0] = int(Array_of_hit_start[j])
    StorageArray[j][1] = int(Array_of_hit_end[j])
    StorageArray[j][2] = int(Array_of_query_start[j])
    StorageArray[j][3] = int(Array_of_query_end[j])
    j = j + 1

print StorageArray[0:20]
QueryFirst = sorted(StorageArray, key = lambda x:x[0], reverse = False)
print("asdfasdfasdf")
print QueryFirst[0:20]

#time to sort
SortedMatches = []
SortedMatches.append(QueryFirst[0])
i=0
j=0
while i <= len(QueryFirst):
    current = SortedMatches[i]
    if i ==1:
        print current
    Distance = 999999
    newDistance = 999998
    j=0
    while j < len(QueryFirst):
        if ((QueryFirst[j][3]>current[3])&(QueryFirst[j][1]>current[1])&(QueryFirst[j][0]>current[0])&(QueryFirst[j][2]>current[2])):
            newDistance = (((QueryFirst[j][1]-current[1])**2)+((QueryFirst[j][3]-current[3])**2))**.5
        if (newDistance < Distance):
            Distance = newDistance
            jDistance = j
        j = j+1
    SortedMatches.append(QueryFirst[jDistance])
    i = i+1
print SortedMatches[0:10][:]

Hit_FASTA = SeqIO.read(open("EColi.fasta"),"fasta")
Query_FASTA = SeqIO.read(open("Salmonella.fasta"),"fasta")
Hit_Strand = str(Hit_FASTA.seq)
Query_Strand = str(Query_FASTA.seq)
file = open("DistanceGene.txt", "w")
k = 0
while k < len(SortedMatches[1]):
    if ((len(Hit_Strand[SortedMatches[k][0]:SortedMatches[k][1]])) & (len(Query_Strand[SortedMatches[k][2]:SortedMatches[k][3]]))&(len(Query_Strand[SortedMatches[k][2]:SortedMatches[k][3]])>10)&(len(Hit_Strand[SortedMatches[k][0]:SortedMatches[k][1]])>10)):
     file.write("\n")
     file.write("\n")
     file.write(">" + HitName + " number_"+str(k+1))
     file.write("_"+HitName" Sequence beginning at "+str(SortedMatches[k][0]) + "_ending at_" + str(SortedMatches[k][1]))
     file.write("\n")
     file.write(Hit_Strand[SortedMatches[k][0]:SortedMatches[k][1]])
     file.write("\n") 
     file.write("\n")
     file.write(">"+QueryName+" number_"+str(k+1))
     file.write("_"+QueryName+" Sequence beginning at "+str(SortedMatches[k][2]) + "_ending at_" + str(SortedMatches[k][3]))
     file.write("\n")
     file.write(Query_Strand[SortedMatches[k][2]:SortedMatches[k][3]])
    k = k +1
file.close()

#add on to write a CSV to R to deal with the data
file= open("genedistance.csv", "w")
k=0
file.write("gapnum,salgapstart,salgapend,salgapseq,ecgapstart,ecgapend,ecgapseq")
file.write("\n")
while k < 2700:
    if ((len(Salmonella_Strand[SortedMatches[k][2]:SortedMatches[k][3]])>10)&(len(EColi_Strand[SortedMatches[k][0]:SortedMatches[k][1]])>10)):
     file.write(str(k+1)+","+str(SortedMatches[k][2])+","+str(SortedMatches[k][3])+","+Salmonella_Strand[SortedMatches[k][2]:SortedMatches[k][3]]+","+str(SortedMatches[k][0])+","+str(SortedMatches[k][1])+","+EColi_Strand[SortedMatches[k][0]:SortedMatches[k][1]])
     file.write("\n")
    k=k+1
file.close()