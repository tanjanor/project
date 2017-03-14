#Description of program: Load classifier and predict label of input data
#Written by: Tanja Normark
#Last Modified:	14.03.2017


location_pssm = input('Enter path for pssm folder: ')

window_size=11
half_window = int((window_size-1)/2)

############################Process input###############################################################

#Get psiblast file names
import os
files = os.listdir(location_pssm)

#create PSSM list (all_prot): containing a list for each protein. 
all_prot = list()
for file_ in files:
	prot_pssm = list()

	#open and read pssm file(s)
	f_path = location_pssm + file_
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


####################################Make windows##########################################################

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

			
####################################Load classifier#######################################################

from sklearn.externals import joblib
clf = joblib.load('./classifier')
	
	
###################################Predict####################################################

labels = clf.predict(all_words)

#Output
map_f = {0:'E',1:'B'}
all_feat = [map_f[char] for char in labels]
print('Predicted labels: ')
print(all_feat)

	
