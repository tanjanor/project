from sklearn import svm
from sklearn.model_selection import train_test_split

f=open("./datasets/buried_exposed_alpha+beta.3line.txt")
list_cont = f.readlines()
window_size=3 ###shorter than each seq

#Put the data into three lists
list_name = list()
list_seq = list()
list_feat = list()

for index in range(0,len(list_cont),3):
	list_name.append(list_cont[index].strip())
for index in range(1,len(list_cont)+1,3):
	list_seq.append(list_cont[index].strip())
for index in range(2,len(list_cont)+1,3):
	list_feat.append(list_cont[index].strip())


#iterate through each aa and express as vector(with neighbours) ->all_words
all_words = list()
aa_dict = {'A':0,'R':1,'N':2,'D':3,'C':4,'E':5,'Q':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
half_window = int((window_size-1)/2)

for seq in list_seq:
	for index0 in range(half_window):			#Beginning of sequence
		numb=half_window-index0
		word_vector= [0]*20*numb
		for pos in range(0,index0+half_window+1):
			aa_num =aa_dict[seq[pos]]
			vector = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
			vector[aa_num]=1
			word_vector.extend(vector)	
		all_words.append(word_vector)
	for index1 in range(half_window,len(seq)-half_window):	#Middle of sequence
		word_vector=list()
		for pos in range(index1-half_window,index1+half_window+1):
			aa_num =aa_dict[seq[pos]]
			vector = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
			vector[aa_num]=1
			word_vector.extend(vector)
		all_words.append(word_vector)
	for index2 in range(len(seq)-half_window, len(seq)):	#End of sequence
		word_vector=list()
		for pos in range(index2-half_window, len(seq)):
			aa_num =aa_dict[seq[pos]]
			vector = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
			vector[aa_num]=1
			word_vector.extend(vector)
		numb=20*window_size-len(word_vector)
		word_vector.extend([0]*numb)
		all_words.append(word_vector)


#iterate through each feature-sequence and express as numbers ->all_feat
all_feat = list()
map_f = {'E':0,'B':1}
for features in list_feat:
	features = [map_f[char] for char in features]
	for feat in features:			
		all_feat.append(feat)


#input into SVM
clf= svm.SVC().fit(all_words,all_feat)
