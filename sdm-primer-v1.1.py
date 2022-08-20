#usage python3.7 sdm-primer-v1.1.py gene.txt mutation.txt
#!/usr/bin/env python
#coding:utf-8


import sys
import os,shutil
import string
from operator import itemgetter
from Bio.Seq import Seq
from Bio.SeqUtils import GC,MeltingTemp


#open mutation amino acid file
currentpath = os.getcwd()
os.chdir(currentpath)
print("Current path is %s." % currentpath)
mutaaFile = open(sys.argv[2])

foldername = "output_primer_of_" + os.path.splitext(sys.argv[1])[0]

mutationlines = mutaaFile.readlines()

for mutaaContent in mutationlines:
    mutaaContent = mutaaContent.strip()
    #print the line length in the mutation amino acid file
    #print("The mutation is %s." % mutaaContent)
    mutaaContentStringLength=len(mutaaContent)
    #print(mutaaContentStringLength)
    #print position of mutation amino acid
    aminoacidPosition = int(mutaaContent[1:mutaaContentStringLength-1])
    print("The amino acid position is %s." % aminoacidPosition)

    #print original amino acid
    originalaminoacid = mutaaContent[0]
    originalaminoacid.lower()
    print("The original amino acid is %s." % originalaminoacid)

    #print new amino acid
    newaminoacid = mutaaContent[mutaaContentStringLength-1]
    print("The new amino acid is %s." % newaminoacid)

    mutationmark = originalaminoacid + str(aminoacidPosition) + newaminoacid

    #read input gene file
    dnaFile = open(sys.argv[1])
    dnaSeq = dnaFile.read()
    dnaFile.close()

    #put into alphabet list
    dnaSeqCompact=[]        

#extract alphabet, convert to uppercase and put into list
    for letter in dnaSeq:
        if letter.isalpha():
            upper_letter=letter.upper()
            dnaSeqCompact.append(upper_letter)
    print("The input gene sequence is")   
    for i in range(0,len(dnaSeqCompact),100):
        DNAstring100 = ''.join(dnaSeqCompact[i:i+100])
        print(DNAstring100)

    print("The lenght of input gene sequence is %s bps." % len(dnaSeqCompact))

#Are gene sequece and length correct?
#for index, letter in enumerate(dnaSeqCompact):
#    print("%i %s" % (index,letter))


#Decide whether the primer design mutation regarding amino acid postion can be done with this program
    if(3*aminoacidPosition < 23 or (len(dnaSeqCompact)-3*aminoacidPosition) < 22):
        print("Position not suitable for computation by this program, please edit the gene sequence by adding the flanking sequence and recalculate amino acid position.")
        sys.exit()
  
#Codon preference table of Escherichia coli
#https://www.genscript.com/tools/codon-frequency-table
#can be edited to other species
    CodonPreferenceTable = {
    'TTT': 'F0.58' , \
    'TTC': 'F0.42' , \
    'TTA': 'L0.14' , \
    'TTG': 'L0.13' , \
    'TAT': 'Y0.59' , \
    'TAC': 'Y0.41' , \
    'TAA': '*0.61' , \
    'TAG': '*0.09' , \
    'CTT': 'L0.12' , \
    'CTC': 'L0.10' , \
    'CTA': 'L0.04' , \
    'CTG': 'L0.47' , \
    'CAT': 'H0.57' , \
    'CAC': 'H0.43' , \
    'CAA': 'Q0.34' , \
    'CAG': 'Q0.66' , \
    'ATT': 'I0.49' , \
    'ATC': 'I0.39' , \
    'ATA': 'I0.11' , \
    'ATG': 'M1.00' , \
    'AAT': 'N0.49' , \
    'AAC': 'N0.51' , \
    'AAA': 'K0.74' , \
    'AAG': 'K0.26' , \
    'GTT': 'V0.28' , \
    'GTC': 'V0.20' , \
    'GTA': 'V0.17' , \
    'GTG': 'V0.35' , \
    'GAT': 'D0.63' , \
    'GAC': 'D0.37' , \
    'GAA': 'E0.68' , \
    'GAG': 'E0.32' , \
    'TCT': 'S0.17' , \
    'TCC': 'S0.15' , \
    'TCA': 'S0.14' , \
    'TCG': 'S0.14' , \
    'TGT': 'C0.46' , \
    'TGC': 'C0.54' , \
    'TGA': '*0.30' , \
    'TGG': 'W1.00' , \
    'CCT': 'P0.18' , \
    'CCC': 'P0.13' , \
    'CCA': 'P0.20' , \
    'CCG': 'P0.49' , \
    'CGT': 'R0.36' , \
    'CGC': 'R0.36' , \
    'CGA': 'R0.07' , \
    'CGG': 'R0.11' , \
    'ACT': 'T0.19' , \
    'ACC': 'T0.40' , \
    'ACA': 'T0.17' , \
    'ACG': 'T0.25' , \
    'AGT': 'S0.16' , \
    'AGC': 'S0.25' , \
    'AGA': 'R0.07' , \
    'AGG': 'R0.04' , \
    'GCT': 'A0.18' , \
    'GCC': 'A0.26' , \
    'GCA': 'A0.23' , \
    'GCG': 'A0.33' , \
    'GGT': 'G0.35' , \
    'GGC': 'G0.37' , \
    'GGA': 'G0.13' , \
    'GGG': 'G0.15'   \
    }

#verify codon preference table
#codon = CodonPreferenceTable['GCT']
#print(codon)

    numberofCodons = 0
    for v in CodonPreferenceTable.items():
        print(v)
        numberofCodons = numberofCodons +1

    print(numberofCodons)
    print('Total codon number is %d.' % numberofCodons)


#Search the dictionary with the first alphabet, put into list and sort by preference
    potentialCodonsList = []               
    numberofPotentialCodons = 0
    potentialCodonString = ""

    for item in CodonPreferenceTable.items():
        tupleValuetoList = list(item)
        aiminoacidPercent = tupleValuetoList[1]
    
        if aiminoacidPercent[0] == newaminoacid:
            potentalCodonString = tupleValuetoList[1][1:] + tupleValuetoList[1][0]+tupleValuetoList[0]
            potentialCodonsList.append(potentalCodonString)
            potentalCodonString = ""

    potentialCodonsList.sort(reverse=True)
    print("The number of potential codons is %i." % len(potentialCodonsList))
    print(potentialCodonsList)
    print("Select the most frequent codon %s for primer design." % potentialCodonsList[0])

        
#Primer selection
    middle = potentialCodonsList[0][5:]

#Get left side and right side
    primers5_GC_within = []
    primers5_GC_out = []
    primers5_Tm_out = []
    steps = [11,12,13,14,15,16,17,18,19,20,21]

    for stepleft in steps:
        for stepright in steps:
            leftstart = 3*aminoacidPosition-4-stepleft+1           #decide interface
            leftend = 3*aminoacidPosition-3
            rightstart = 3*aminoacidPosition
            rightend = 3*aminoacidPosition+stepright
        
            left = dnaSeqCompact[leftstart:leftend]                #left and right slices
            right = dnaSeqCompact[rightstart:rightend]
            leftstring =''.join(left)
            rightstring =''.join(right)
            leftseq = Seq(leftstring)
            rightseq = Seq(rightstring)

#Decide if the left and right slices are okay with GC/AT thumb of rule
            leftTm = 4*(leftseq.count("G")+leftseq.count("C"))+2*(leftseq.count("A")+leftseq.count("T"))
            rightTm = 4*(rightseq.count("G")+rightseq.count("C"))+2*(rightseq.count("A")+rightseq.count("T"))

            if((leftTm >= 45) and (rightTm >= 45)):
                addPrimer5 = leftstring+middle+rightstring
                addPrimerSeq = Seq(addPrimer5)
                addPrimerSeq_GCcontent = round(GC(addPrimerSeq),2)
                addPrimerSeq_GCcontent_string = str(addPrimerSeq_GCcontent)
                if(addPrimerSeq_GCcontent >= 40 and addPrimerSeq_GCcontent <= 60):
                    primers5_GC_within.append([addPrimer5,stepleft,stepright,round(MeltingTemp.Tm_NN(addPrimerSeq), 2)])
                else:
                    primers5_GC_out.append([addPrimer5,stepleft,stepright,round(MeltingTemp.Tm_NN(addPrimerSeq), 2)])

            else:
                addPrimer5 = leftstring+middle+rightstring
                addPrimerSeq = Seq(addPrimer5)
                primers5_Tm_out.append([addPrimer5,stepleft,stepright,round(MeltingTemp.Tm_NN(addPrimerSeq), 2)])
            
          
    primers5_GC_within_sorted = sorted(primers5_GC_within, key = itemgetter(3,0))
    primers5_GC_out_sorted = sorted(primers5_GC_out, key = itemgetter(3,0))
    primers5_Tm_out_sorted = sorted(primers5_Tm_out, key = itemgetter(3,0))


    
    for i in range(len(primers5_GC_within_sorted)):
        print(primers5_GC_within_sorted[i])
    print("The number of optional forward primers with suitable 4*GC+2*AT and GC content is %i." % len(primers5_GC_within_sorted))
    for i in range(len(primers5_GC_out_sorted)):
        print(primers5_GC_out_sorted[i])
    print("The number of sub-optional forward primers with suitable 4*GC+2*AT but not GC content is %i." % len(primers5_GC_out_sorted))
    for i in range(len(primers5_Tm_out_sorted)):
        print(primers5_Tm_out_sorted[i])
    print("The number of forward primers with 4*GC+2*AT out of range is %i." % len(primers5_Tm_out_sorted))


#Collect, sort primers in accordance with 2 conditions, ATGC rule of thumb and GC content

    row = 2*len(primers5_GC_within_sorted)
    column = 10

    primers53 =[[[] for j in range(column)] for i in range(row)]

    for i in range(row):

        primer_forward_String = primers5_GC_within_sorted[int(i/2)][0]
        primer_forward_Seq = Seq(primer_forward_String)

        primer_reversecomplement = primer_forward_Seq.reverse_complement()
        primer_reversecomplement_String = ''.join(primer_reversecomplement)
        primer_reversecomplement_Seq = Seq(primer_reversecomplement_String)

        primers53[i][2] = "Length:"
        primers53[i][4] = "GC:"
        primers53[i][6] = "Tm:"
        primers53[i][8] = primers5_GC_within_sorted[int(i/2)][1]
        primers53[i][9] = primers5_GC_within_sorted[int(i/2)][2]
        if(i%2 == 0):
            primers53[i][0] = mutationmark+"_forward_"+str(i)+"_5: "
            primers53[i][1] = primer_forward_String
            primers53[i][3] = len(primer_forward_String)                                                       
            primers53[i][5] = round(GC(primer_forward_Seq),2)                                                                       
            primers53[i][7] = round(MeltingTemp.Tm_NN(primer_forward_Seq), 2)                                                   
       
        else:
            primers53[i][0] = mutationmark+"_reverse_"+str(i)+"_3: "
            primers53[i][1] = primer_reversecomplement_String
            primers53[i][3] = len(primer_reversecomplement_String)
            primers53[i][5] = round(GC(primer_reversecomplement_Seq),2)
            primers53[i][7] = round(MeltingTemp.Tm_NN(primer_reversecomplement_Seq), 2)
             

    for i in range(len(primers53)):
        print(primers53[i])
    print("The number of optional primer pairs is %i." % int(len(primers53)/2))

#write file
    file=open(mutationmark+'_output_primer_optional.txt','w')

    for i in range(len(primers53)):
        file.write(str(primers53[i])+'\r\n')
    file.close()


##Collect, sort primers in accordance with ATGC rule of thumb, but not GC content

    rowsub = 2*len(primers5_GC_out_sorted)
    columnsub = 10

    primers53suboptional = [[[] for j in range(columnsub)] for i in range(rowsub)]

    for i in range(rowsub):

        primer_forward_String = primers5_GC_out_sorted[int(i/2)][0]
        primer_forward_Seq = Seq(primer_forward_String)

        primer_reversecomplement = primer_forward_Seq.reverse_complement()
        primer_reversecomplement_String = ''.join(primer_reversecomplement)
        primer_reversecomplement_Seq = Seq(primer_reversecomplement_String)

        primers53suboptional[i][2] = "Length:"
        primers53suboptional[i][4] = "GC:"
        primers53suboptional[i][6] = "Tm:"
        primers53suboptional[i][8] = primers5_GC_out_sorted[int(i/2)][1]
        primers53suboptional[i][9] = primers5_GC_out_sorted[int(i/2)][2]
    
        if(i%2 == 0):
            primers53suboptional[i][0] = mutationmark+"_forward_"+str(i)+"_5: "
            primers53suboptional[i][1] = primer_forward_String
            primers53suboptional[i][3] = len(primer_forward_String)                                                       
            primers53suboptional[i][5] = round(GC(primer_forward_Seq),2)                                                                       
            primers53suboptional[i][7] = round(MeltingTemp.Tm_NN(primer_forward_Seq), 2)                                                  
        
        else:
            primers53suboptional[i][0] = mutationmark+"_reverse_"+str(i)+"_3: "
            primers53suboptional[i][1] = primer_reversecomplement_String
            primers53suboptional[i][3] = len(primer_reversecomplement_String)
            primers53suboptional[i][5] = round(GC(primer_reversecomplement_Seq),2)
            primers53suboptional[i][7] = round(MeltingTemp.Tm_NN(primer_reversecomplement_Seq), 2)
     
    for i in range(len(primers53suboptional)):
        print(primers53suboptional[i])
    print("The number of sub-optional primer pairs is %i." % int(len(primers53suboptional)/2))

#write file
    file=open(mutationmark+'_output_primer_suboptional.txt','w')

    for i in range(len(primers53suboptional)):
        file.write(str(primers53suboptional[i])+'\r\n')
    file.close()


dirs = os.listdir(currentpath)

if(foldername not in dirs):
    os.system('mkdir temp_foldername')
    os.system('mv *optional* temp_foldername')
    os.rename("temp_foldername", foldername)

else:
    shutil.rmtree(foldername)
    os.system('mkdir temp_foldername')
    os.system('mv *optional* temp_foldername')
    os.rename("temp_foldername", foldername)


    

