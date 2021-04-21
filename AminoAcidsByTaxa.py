import numpy as np 
import Bio
import pandas as pd
from Bio import SeqIO
import requests
from Bio import Seq
from Bio.Align.Applications import MuscleCommandline
from io import StringIO
from Bio import AlignIO
import os.path
from os import path
mtSSU_proteins = ["S2","S24","S5","S6","S7","S9","S10","S11","S12","S14","S15","S16","S17","S18C","S21","S22","S23","S25","S26","S27","S28","S29","S31","S33","S34","S35","S37","S38","S39","S18B"]
mtLSU_proteins = ["L2","L3","L4","L9","L10","L11","L13","L14","L15","L16","L17","L18","L19","L20","L21","L22","L23","L24","L30","S30","L37","L38","L39","L40","L41","L42","L43","L44","L45","L46","L48","L49","L50","L51","L52","L53","L54","L57","L58","L59","L18A"]
mt_name_map = {"S2":"US2M","S24":"US3M","S5":"US5M","S6":"BS6M","S7":"US7M","S9":"US9M","S10":"US10M","S11":"US11M","S12":"US12M","S14":"US14M","S15":"US15M","S16":"BS16M","S17":"US17M","S18C":"BS18M","S21":"BS21M","S22":"MS22","S23":"MS23","S25":"MS25","S26":"MS26","S27":"MS27","S28":"MS28","S29":"MS29","S31":"MS31","S33":"MS33","S34":"MS34","S35":"MS35","S37":"MS37","S38":"MS38","S39":"MS39","S18B":"MS40","L2":"UL2M","L3":"UL3M","L4":"UL4M","L9":"BL9M","L10":"UL10M","L11":"UL11M","L13":"UL13M","L14":"UL14M","L15":"UL15M","L16":"UL16M","L17":"BL17M","L18":"UL18M","L19":"BL19M","L20":"BL20M","L21":"BL21M","L22":"UL22M","L23":"UL23M","L24":"UL24M","L30":"UL30M","S30":"ML65","L37":"ML37","L38":"ML38","L39":"ML39","L40":"ML40","L41":"ML41","L42":"ML42","L43":"ML43","L44":"ML44","L45":"ML45","L46":"ML46","L48":"ML48","L49":"ML49","L50":"ML50","L51":"ML51","L52":"ML52","L53":"ML53","L54":"ML54","L57":"ML63","L58":"ML62","L59":"ML64","L18A":"ML66"}
# map in the format of {human name:new name}, not applicable to mitochondrial
cytosolic_name_map = {"S3A":"eS1","SA":"uS2","S3":"uS2","S9":"uS4","S4":"eS4","S2":"uS5","eS6":"S6","S5":"uS7","S7":"eS7","S15A":"uS8","S8":"eS8","S16":"uS9","S20":"uS10","S10":"eS10","S14":"uS11","S23":"uS12","S12":"eS12","S18":"uS13","S29":"uS14","S13":"uS15","S11":"uS17","S17":"eS17","S15":"uS19","S19":"eS19","S21":"eS21","S24":"eS24","S25":"eS25","S26":"eS26","S27":"eS27","S28":"eS28","S30":"eS30","S27A":"eS31","RACK1":"RACK1","L10A":"uL1","L8":"uL2","L3":"uL3","L4":"uL4","L11":"uL5","L9":"uL6","L6":"eL6","L7A":"eL8","P0":"uL10","L12":"uL11","L13A":"uL13","L13":"eL13","L23":"uL14","L14":"eL14","L27A":"uL15","L15":"el15","L10":"uL16","L5":"uL18","L18":"eL18","L19":"eL19","L18A":"eL20","L21":"eL21","L17":"uL22","L22":"eL22","L23A":"uL23","L26":"uL24","L24":"eL24","L27":"eL27","L28":"eL28","L35":"uL29","L29":"eL29","L7":"uL30","L30":"eL30","L31":"eL31","L32":"eL32","L35A":"eL33","L34":"eL34","L36":"eL36","L37":"eL37","L38":"eL38","L39":"eL39","L40":"eL40","L41":"eL41","L36A":"eL42","L37A":"eL43"}

all_organisms = [('mammalian','40674'),('birds','8782'),('marsupials','(38605+OR+38606+OR+38607+OR+38608+OR+38610+OR+38611+OR+38609)'),('plants','33090'),('reptiles','(8459+OR+8504+OR+8505+OR+8509+OR+1294634)')]

#get sequences from Uniprot by their taxonomy ID and protein (S2, S3, etc) return a dataframe 
def get_seqs(taxid, entry_name1, entry_name2=None,bacterial=False):
   
    if entry_name2 is None:
        additional_str = ""
    else:
        additional_str = "+OR+" + entry_name2 
    
    if not bacterial:
        uniprot_url = f"https://www.uniprot.org/uniprot/?query=name:{entry_name1}{additional_str}+AND+taxonomy:{taxid}&sort=score&columns=id,protein names,organism-id,length,reviewed,organism&format=tab"
    else:
        uniprot_url = f"https://www.uniprot.org/uniprot/?query=name:{entry_name1}{additional_str}+AND+taxonomy:{taxid}+AND+reviewed:yes&sort=score&columns=id,protein names,organism-id,length,reviewed,organism&format=tab"
    s = requests.get(uniprot_url).content
    s = s.decode()
    lst = [sl.split('\t') for sl in s.split('\n')]
    header = lst[0]
    values = lst[1:]
    df = pd.DataFrame(values, columns=header).dropna()
    return df


def seqs_to_csv(seqs_df,taxonomy, protein_name, subunit, filtered=False, is_cytosolic=False):
    if is_cytosolic:
        seqs_filename = f"{cytosolic_name_map[protein_name]}_{str(taxonomy)}_cytosolic_{subunit}.csv"
    else:
        seqs_filename = f"{mt_name_map[protein_name]}_{str(taxonomy)}_mitochondrial_{subunit}.csv"
    if filtered:
        seqs_filename = f"filtered{seqs_filename}"
    if not path.exists(seqs_filename):
        seqs_df.to_csv(seqs_filename)
    

# filter data in a dataframe of those sequences, returning only the trustworthy sequences
def filter_seqs(df, protein_num, is_cytosolic_ribosome=False):
    # filter by seq length and the entry name
    df['Protein names'] = df['Protein names'].str.lower()
    protein_num = protein_num.lower()
    if (is_cytosolic_ribosome):
        condition = (df['Length'].astype(int) > 10) & (df['Length'].astype(int) < 1000) &\
        (df['Protein names'].str.contains("fragment") == False) &\
        (df['Protein names'].str.contains("like") == False) &\
        ((df['Protein names'].str.contains("40s") == True) | (df['Protein names'].str.contains("60s") == True) |\
        (df['Protein names'].str.contains("30s") == True) | (df['Protein names'].str.contains("50s") == True)) &\
        (df['Protein names'].str.contains(protein_num) == True)
    else:
        condition = (df['Length'].astype(int) > 10) & (df['Length'].astype(int) < 1000) &\
        (df['Protein names'].str.contains("fragment") == False) &\
        (df['Protein names'].str.contains("isoform") == False) &\
        (df['Protein names'].str.contains("like") == False) &\
        (df['Protein names'].str.contains("chloroplastic") == False) &\
        ((df['Protein names'].str.contains("mitochondrial r") == True) | (df['Protein names'].str.contains("28s") == True) | (df['Protein names'].str.contains("39s") == True))  &\
        (df['Protein names'].str.contains(protein_num) == True)
    df2 = df[condition].reset_index(drop=True)

    # For the entries with the same organism                                                                                  
    df2['Length'] = df2['Length'].astype(int)
    df2 = df2.sort_values(['Status','Length'],ascending= [False,True])
    df3 = df2.drop_duplicates(subset=['Organism'],keep='last').reset_index(drop=True) # drop entries with same organism, keephighest length                                                                                                               
    df3 = df3.sort_values('Status').reset_index(drop=True)
    return df3;



# create a fasta file of the retrieved Uniprot sequences
def create_fasta(taxonomy, protein_name, cellular_location,  seqs_df, cytosolic=False, bacterial=False):
    nrows, ncols = seqs_df.shape 
    #create a fasta file 
    if not cytosolic and not bacterial:
        ofile = open(mt_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".fasta", "w")
    elif cytosolic and not bacterial:
        ofile = open(cytosolic_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".fasta", "w")
    elif bacterial and not cytosolic:
         ofile = open(protein_name + "_" + str(taxonomy) + "_" + cellular_location + ".fasta", "w")
    for i in range(nrows):
        new_seq_url = "https://uniprot.org/uniprot/" + seqs_df.iat[i,0] + ".fasta"
        seq_url_get = requests.get(new_seq_url)
        protein_seq = seq_url_get.content
        ofile.write(protein_seq.decode())
    ofile.close()

def create_alignment(taxonomy, protein_name, cellular_location,  seqs_df,cytosolic=False):
    #create alignment file with muscle
    if not cytosolic:
        muscle_cline = MuscleCommandline("./muscle",input=mt_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".fasta")
        stdout, stderr = muscle_cline()
        alignment = AlignIO.read(StringIO(stdout),"fasta")
        AlignIO.write(alignment, mt_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".sto", "stockholm")
    else:
        muscle_cline = MuscleCommandline("./muscle",input=cytosolic_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".fasta")
        stdout, stderr = muscle_cline()
        alignment = AlignIO.read(StringIO(stdout),"fasta")
        AlignIO.write(alignment, cytosolic_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".sto", "stockholm")


def parse_to_seqrecords(taxonomy, protein_name, cellular_location,  cytosolic=False,bacterial=False):
    #parse the file to seq records
    if not cytosolic and not bacterial:
        protein_seq_records = list(SeqIO.parse(mt_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".fasta","fasta"))
    elif cytosolic and not bacterial:
        protein_seq_records = list(SeqIO.parse(cytosolic_name_map[protein_name] + "_" + str(taxonomy) + "_" + cellular_location + ".fasta","fasta"))
    elif bacterial and not cytosolic:
        protein_seq_records = list(SeqIO.parse(protein_name + "_" + str(taxonomy) + "_" + cellular_location + ".fasta","fasta"))
    return protein_seq_records

def parse_seqs(taxonomy, protein_name, cellular_location,  seqs_df, cytosolic=False,bacterial=False):
    if not cytosolic and not bacterial:
        create_fasta(taxonomy, protein_name, cellular_location,  seqs_df)
        if seqs_df.shape[0] != 1:
            create_alignment(taxonomy, protein_name, cellular_location,  seqs_df)
        seq_records = parse_to_seqrecords(taxonomy, protein_name, cellular_location)
    
    elif cytosolic and not bacterial:
        create_fasta(taxonomy, protein_name, cellular_location,  seqs_df,cytosolic=True)
        if seqs_df.shape[0] != 1:
            create_alignment(taxonomy, protein_name, cellular_location,  seqs_df,cytosolic=False,bacterial=True)
        seq_records = parse_to_seqrecords(taxonomy, protein_name, cellular_location, cytosolic=True)
    
    elif bacterial and not cytosolic:
        create_fasta(taxonomy, protein_name, cellular_location,  seqs_df,cytosolic=False,bacterial=True)
        #if seqs_df.shape[0] != 1:
        #    create_alignment(taxonomy, protein_name, cellular_location,  seqs_df,cytosolic=False,bacterial=True)
        seq_records = parse_to_seqrecords(taxonomy, protein_name, cellular_location, cytosolic=False,bacterial=True)
    return seq_records

def stdev(lst):
    if len(lst) == 1:
        return 0
    mean = sum(lst) / len(lst)
    sd = (sum((x-mean)**2.0 for x in lst) / float(len(lst)-1))**0.5
    return sd

# the total sum of electron-transferring amino acids (cysteine, tyrptophan, methionine, and tyrosine) in the list of seq records
def sum_CYWM(seqs_records):
    num_aminoacids = 0
    cys_lst = []
    trp_lst = []
    met_lst = []
    tyr_lst = []
    met = 0
    cys = 0
    trp = 0
    tyr = 0
    sumofCYWM = 0
    met
    num_seqs = len(seqs_records)
    sumofCYWM_lst = []
    for i in range(num_seqs):
        cys_lst.append(seqs_records[i].seq.count("C") / len(seqs_records[i])*100)
        trp_lst.append(seqs_records[i].seq.count("W") / len(seqs_records[i])*100)
        met_lst.append(seqs_records[i].seq.count("M") / len(seqs_records[i])*100)
        tyr_lst.append(seqs_records[i].seq.count("Y") / len(seqs_records[i])*100)
        sumofCYWM_lst.append((seqs_records[i].seq.count("C") + seqs_records[i].seq.count("Y")  + seqs_records[i].seq.count("W") +seqs_records[i].seq.count("M")) / len(seqs_records[i])*100)
        num_aminoacids += len(seqs_records[i])
      
    sumofCYWM = round(sum(sumofCYWM_lst)/num_seqs,1)
    cys = round(sum(cys_lst)/num_seqs,1)
    tyr = round(sum(tyr_lst)/num_seqs,1)
    trp = round(sum(trp_lst)/num_seqs,1)
    met = round(sum(met_lst)/num_seqs,1)

    return cys,round(stdev(cys_lst),1),met,round(stdev(met_lst),1),tyr,round(stdev(tyr_lst),1),trp,round(stdev(trp_lst),1),sumofCYWM,round(stdev(sumofCYWM_lst),1),num_seqs,num_aminoacids

# the total sum of electron-transferring amino acids (cysteine, tyrptophan, methionine, and tyrosine) in the list of seq records
def sum_CYWM2(seqs_records):
    num_aminoacids = 0
    cys_lst = []
    trp_lst = []
    met_lst = []
    tyr_lst = []
    met = 0
    cys = 0
    trp = 0
    tyr = 0
    sumofCYWM = 0
    met
    num_seqs = len(seqs_records)
    sumofCYWM_lst = []
    for i in range(num_seqs):
        cys_lst.append(seqs_records[i].seq.count("C") / len(seqs_records[i])*100)
        trp_lst.append(seqs_records[i].seq.count("W") / len(seqs_records[i])*100)
        met_lst.append(seqs_records[i].seq.count("M") / len(seqs_records[i])*100)
        tyr_lst.append(seqs_records[i].seq.count("Y") / len(seqs_records[i])*100)
        num_aminoacids += len(seqs_records[i])
      
    sumofCYWM_lst = cys_lst + trp_lst + met_lst + tyr_lst    
    sumofCYWM = round(sum(sumofCYWM_lst)/num_seqs,1)
    cys = round(sum(cys_lst)/num_seqs,1)
    tyr = round(sum(tyr_lst)/num_seqs,1)
    trp = round(sum(trp_lst)/num_seqs,1)
    met = round(sum(met_lst)/num_seqs,1)

    return cys,round(stdev(cys_lst),1),met,round(stdev(met_lst),1),tyr,round(stdev(tyr_lst),1),trp,round(stdev(trp_lst),1),sumofCYWM,round(stdev(sumofCYWM_lst),1),num_seqs,num_aminoacids

# the total sum of electron-transferring amino acids (cysteine, tyrptophan, methionine, and tyrosine) in the list of seq records
def sum_CYWM3(seqs_records):
    num_aminoacids = 0
    cys_lst = []
    trp_lst = []
    met_lst = []
    tyr_lst = []
    met = 0
    cys = 0
    trp = 0
    tyr = 0
    sumofCYWM = 0
    met
    num_seqs = len(seqs_records)
    sumofCYWM_lst = []
    for i in range(num_seqs):
        sumofCYWM_lst.append((seqs_records[i].seq.count("C") + seqs_records[i].seq.count("Y")  + seqs_records[i].seq.count("W") +seqs_records[i].seq.count("M")) / len(seqs_records[i])*100)
        num_aminoacids += len(seqs_records[i])

    return sumofCYWM_lst,round(stdev(sumofCYWM_lst),1),num_seqs


def append_final_df(final_df, taxonomy, cellular_location, cys,cys_sd,met,met_sd,tyr,tyr_sd,trp,trp_sd,sumCYWM,sumCYWMsd,num_seqs,num_aminoacids):
    new_data = {"Organisms":taxonomy,"Cellular location":cellular_location,"Cysteine":str(cys) +"%", "CYS SD":u"\u00B1" + str(cys_sd), "Tyrosine":str(tyr) + "%","TYR SD":u"\u00B1" + str(tyr_sd), "Tryptophan":str(trp) + "%","TRP SD":u"\u00B1" + str(trp_sd),"Methionine":str(met) + "%","MET SD":u"\u00B1" + str(met_sd),"Sum of C,Y,W,M":str(sumCYWM) + "%","Sum SD":u"\u00B1" + str(sumCYWMsd),"Total Number of Sequences":num_seqs,"Total Number of Amino Acids":num_aminoacids}
    final_df = final_df.append(new_data, ignore_index=True)
    return final_df
def aminoacidcomp(seqs_records):
    num_aminoacids = 0
    cys_lst = []
    trp_lst = []
    met_lst = []
    tyr_lst = []
    ala_lst = []
    arg_lst = []
    asn_lst = []
    asx_lst = []
    glu_lst = []
    glx_lst = []
    gly_lst = []
    his_lst = []
    lle_lst = []
    leu_lst = []
    lys_lst = []
    thr_lst = []
    phe_lst = []
    pro_lst = []
    ser_lst = []
    val_lst = []
    sumofCYWM = 0
    met
    num_seqs = len(seqs_records)
    sumofCYWM_lst = []
    for i in range(num_seqs):
        cys_lst.append(seqs_records[i].seq.count("C") / len(seqs_records[i])*100)
        trp_lst.append(seqs_records[i].seq.count("W") / len(seqs_records[i])*100)
        met_lst.append(seqs_records[i].seq.count("M") / len(seqs_records[i])*100)
        tyr_lst.append(seqs_records[i].seq.count("Y") / len(seqs_records[i])*100)
        cys_lst.append(seqs_records[i].seq.count("C") / len(seqs_records[i])*100)
        trp_lst.append(seqs_records[i].seq.count("W") / len(seqs_records[i])*100)
        met_lst.append(seqs_records[i].seq.count("M") / len(seqs_records[i])*100)
        tyr_lst.append(seqs_records[i].seq.count("Y") / len(seqs_records[i])*100)
        cys_lst.append(seqs_records[i].seq.count("C") / len(seqs_records[i])*100)
        trp_lst.append(seqs_records[i].seq.count("W") / len(seqs_records[i])*100)
        met_lst.append(seqs_records[i].seq.count("M") / len(seqs_records[i])*100)
        tyr_lst.append(seqs_records[i].seq.count("Y") / len(seqs_records[i])*100)
        cys_lst.append(seqs_records[i].seq.count("C") / len(seqs_records[i])*100)
        trp_lst.append(seqs_records[i].seq.count("W") / len(seqs_records[i])*100)
        met_lst.append(seqs_records[i].seq.count("M") / len(seqs_records[i])*100)
        tyr_lst.append(seqs_records[i].seq.count("Y") / len(seqs_records[i])*100)
        cys_lst.append(seqs_records[i].seq.count("C") / len(seqs_records[i])*100)
        trp_lst.append(seqs_records[i].seq.count("W") / len(seqs_records[i])*100)
        met_lst.append(seqs_records[i].seq.count("M") / len(seqs_records[i])*100)
        tyr_lst.append(seqs_records[i].seq.count("Y") / len(seqs_records[i])*100)
        
        num_aminoacids += len(seqs_records[i])
      
    sumofCYWM_lst = cys_lst + trp_lst + met_lst + tyr_lst    
    sumofCYWM = round(sum(sumofCYWM_lst)/num_seqs,1)
    cys = round(sum(cys_lst)/num_seqs,1)
    tyr = round(sum(tyr_lst)/num_seqs,1)
    trp = round(sum(trp_lst)/num_seqs,1)
    met = round(sum(met_lst)/num_seqs,1)

    return cys,round(stdev(cys_lst),1),met,round(stdev(met_lst),1),tyr,round(stdev(tyr_lst),1),trp,round(stdev(trp_lst),1),sumofCYWM,round(stdev(sumofCYWM_lst),1),num_seqs,num_aminoacids
def append_final_df(final_df, taxonomy, cellular_location, cys,cys_sd,met,met_sd,tyr,tyr_sd,trp,trp_sd,phe,phe_sd,sumCYWM,sumCYWMsd,num_seqs,num_aminoacids):
    new_data = {"Organisms":taxonomy,"Cellular location":cellular_location,"Cysteine":str(cys) +"%", "CYS SD":u"\u00B1" + str(cys_sd), "Tyrosine":str(tyr) + "%","TYR SD":u"\u00B1" + str(tyr_sd), "Tryptophan":str(trp) + "%","TRP SD":u"\u00B1" + str(trp_sd),"Methionine":str(met) + "%","MET SD":u"\u00B1" + str(met_sd),"Sum of C,Y,W,M":str(sumCYWM) + "%","Sum SD":u"\u00B1" + str(sumCYWMsd),"Total Number of Sequences":num_seqs,"Total Number of Amino Acids":num_aminoacids}
    final_df = final_df.append(new_data, ignore_index=True)
    return final_df
def append_final_df2(final_df, protein, cys,cys_sd,met,met_sd,tyr,tyr_sd,trp,trp_sd,sumCYWM,sumCYWMsd,num_seqs,num_aa):
    new_data = {"Protein":protein,"Cysteine":str(cys) +"% " + u"\u00B1" + str(cys_sd), "Tyrosine":str(tyr) + "% "+ u"\u00B1" + str(tyr_sd), "Tryptophan":str(trp) + "% " + u"\u00B1" + str(trp_sd),"Methionine":str(met) + "% "+ u"\u00B1" + str(met_sd),"Sum of C,Y,W,M":str(sumCYWM) + "% "+ u"\u00B1" + str(sumCYWMsd),"Total Number of Sequences":num_seqs,"Total Number of Amino Acids":num_aa}
    final_df = final_df.append(new_data, ignore_index=True)
    return final_df
