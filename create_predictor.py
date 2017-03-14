#Description of program: Train a classifier for exposed/buried prediction
#Written by: Tanja Normark
#Last Modified:	14.03.2017

#Note: Works for filenames with a name like '>d1a26_1.a.41.1.1.fasta.pssm'
from sklearn import svm


protein_data_file = "./datasets/buried_exposed_alpha+beta.3line.txt"
pssm_folder = './input/pssm'


f=open(protein_data_file)
list_cont = f.readlines()
window_size=11
half_window = int((window_size-1)/2)

#Put the protein data into three lists
list_name = list()
list_seq = list()
list_feat = list()

for index in range(0,len(list_cont),3):
	list_name.append(list_cont[index].strip())
for index in range(1,len(list_cont)+1,3):
	list_seq.append(list_cont[index].strip())
for index in range(2,len(list_cont)+1,3):
	list_feat.append(list_cont[index].strip())

#create a dictionary
d_data = dict()
for index in range(len(list_name)):
	d_data[list_name[index]]=(list_seq[index], list_feat[index])

#Get psiblast file names
import os
files = os.listdir(pssm_folder)
	
#create PSSM list (all_prot): containing a list for each protein. 
all_prot = list()
all_labels = list()
for file_ in files:
	prot_pssm = list()

	#get labels for protein
	prot_name = file_[:-11]	
	labels = d_data[prot_name][1]
	all_labels.extend(labels)

	#open and read pssm file
	f_path = './input/pssm/' + file_
	pssm_f = open(f_path)
	pssm_lines = pssm_f.read().split('\n')
		
	#making a vector for each line	
	for index in range(3,len(pssm_lines)-7):
		all_items_vector = list()			
		line = pssm_lines[index]
		items_line = line.split()
		for item in items_line:
			all_items_vector.append(item)		
		vector = all_items_vector[22:42]
			
		#normalization by dividing each by 100		
		for i in range(len(vector)):
			vector[i] = int(vector[i])/100

		#add vector to protein vector
		prot_pssm.append(vector)
	
	#add protein vector to all vector
	all_prot.append(prot_pssm)

#convert all_prot to vectors -> all_words:
all_words = list()
for prot in all_prot:
	for index0 in range(half_window):			#Beginning of sequence
		numb=half_window-index0
		word_vector= [0.0]*20*numb	
		for pos in range(0,index0+half_window+1):
			word_vector.extend(prot[pos])		
		all_words.append(word_vector)
	for index1 in range(half_window,len(prot)-half_window):	#Middle of sequence
		word_vector=list()
		for pos in range(index1-half_window,index1+half_window+1):
			word_vector.extend(prot[pos])	
		all_words.append(word_vector)
	for index2 in range(len(prot)-half_window, len(prot)):	#End of sequence	
		word_vector=list()
		for pos in range(index2-half_window, len(prot)):
			word_vector.extend(prot[pos])
		numb=20*window_size-len(word_vector)
		word_vector.extend([0.0]*numb)
		all_words.append(word_vector)


#iterate through each feature-sequence and express as numbers ->all_feat
all_feat = list()
map_f = {'E':0,'B':1}
all_feat = [map_f[char] for char in all_labels]

#input into SVM
clf = svm.SVC(kernel='rbf',C=2,gamma=0.125).fit(all_words,all_feat)

#Save classifier
from sklearn.externals import joblib
joblib.dump(clf,'./classifier')
