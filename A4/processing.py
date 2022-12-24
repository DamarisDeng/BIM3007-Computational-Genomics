import pandas as pd

def extract(seq, pos):
    '''Given a sequence and a position, return a 13-mer centered on the position.'''
    if pos < 7:
        seq = 'X'*(7-pos) + seq
        pos = 7
    elif len(seq) - pos < 6:
        seq = seq + 'X'*(6 - (len(seq) - pos))  
    return seq[pos-7:pos+6]

data = pd.read_csv('Ubiquitination_sites.txt', sep='\t')
data1 = data.drop(['PLMD ID', 'Species', 'PMIDs', 'Type', 'Uniprot Accession'], axis=1)
position = data1['Position']

positive = list()
for i in range(8000):
    sequence = data1['Sequence'][i]
    pos = position[i] 
    positive.append(extract(sequence, pos))

data1['positive'] = positive

positive = open('positive.txt', 'w')
for i in range(8000):
    positive.write('> '+str(i) + '\n')
    positive.write(data1['positive'][i] + '\n')
positive.close()

# merge the identical sequences
data2 = pd.DataFrame(data1.groupby('Sequence').apply(lambda x: x['Position'].tolist()))
data2.reset_index(inplace=True)
data2.columns = ['sequence', 'positive'] # change names
sequences = data2['sequence']
all_K_position = []
for sequence in sequences:
    # identify all the positions of K
    K_position = [i for i in range(len(sequence)) if sequence[i] == 'K']
    all_K_position.append(K_position)
data2['all'] = all_K_position

negatives = []
for i in range(len(data2)):
    positive = data2['positive'][i]
    all_position = data2['all'][i]
    negative = [a+1 for a in all_position if all(abs(a - b) > 6 for b in positive)] # remove negatives that are too close to the positive
    negatives.append(negative)
data2['negative'] = negatives
negative = open('negative.txt', 'w')
# output positive and negative sequences
i = 0
for row in range(len(sequences)):
    # write negatives
    sequence = sequences[row]
    for pos in data2['negative'][row]: 
        negative.write('> '+ str(i) + '\n')
        negative.write(extract(sequence, pos) + '\n')
        i += 1
negative.close()