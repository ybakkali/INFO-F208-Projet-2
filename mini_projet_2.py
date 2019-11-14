def matrixParser(filename):
    f = open(filename, encoding='utf-8')
    matrix = SubstitutionMatrix()
    for line in f:
        if line[0] not in ["#", " ", "\n"]:
            matrix.Append(line[0])
            matrix.addRow(line[1:].split())
    f.close()
    return matrix

def sequencesParser(filename):
    f = open(filename, encoding='utf-8')
    mylist = []
    sequence = ''
    for line in f:
        if line[0] not in [">", " ", "\n"]:
            sequence += line.strip("\n")
        elif len(sequence) :
            mylist.append(ADTsequence(sequence))
            sequence = ''
    if len(sequence) :
        mylist.append(ADTsequence(sequence))
    f.close()
    return mylist

class ADTsequence(object):

    def __init__(self, sequence):
        self.sequence = sequence

    def __getitem__(self, key):
        return self.sequence[key]

    def __len__(self):
        return len(self.sequence)

    def __str__(self):
        output = ""
        for i in range(len(self) // 60 + 1):
            output += self.sequence[i * 60:(i + 1) * 60]
            output += "\n"
        return output[:-1]


class ADTmatrix:

    def __init__(self, n=None, m=None):

        if n and m:
            self.matrix = [[0 for i in range(m)] for j in range(n)]
        else:
            self.matrix = []

    def getItem(self, i, j):
        return self.matrix[i][j]

    def setItem(self, i, j, value):
        self.matrix[i][j] = value

    def addRow(self, row):
        self.matrix.append(row)

    def __str__(self):
        """
        The method for printing the matrix in a pretty form
        """
        lst = []
        for line in self.matrix :#for each line
        		lst.append([str(elem) for elem in line]) #we take the elements and change them to str format
        dist = []
        for i in zip(*lst) : #adjusting space
        	dist.append(max(map(len, i))) #adjusting space
        s = '   '.join('{{:{}}}'.format(e) for e in dist) #adjusting space
        out = []
        for elem in lst: #for each line that we modified in first step
        	out.append(s.format(*elem)) #preparing the out put
        print('\n'.join(out))
        return ""
    """
    def __str__(self):
        output = ""
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                output += str(self.matrix[i][j])
                output += " "
            output += "\n"
        return output[:-1]
    """


class SubstitutionMatrix(ADTmatrix):

    def __init__(self, n=None, m=None):
        super().__init__(n, m)
        self.myList = []

    def __getitem__(self, key):
        return self.myList.index(key)

    def Append(self, value):
        self.myList.append(value)

class NeedlemanWunsch(object):
    def __init__(self, S, V, W, MATSUB, k, I, E):

        self.S = S
        self.V = V
        self.W = W
        self.MATSUB = MATSUB
        self.k = k
        self.I = I
        self.E = E
        self.solutionCounter = 0
        self.score = 0
        self.similarityRate = 0

        for i in range(1, self.S.n):
            for j in range(1, self.S.m):
                self.fill(i, j)
        self.score = self.S.getItem(self.S.n - 1, self.S.m - 1)
        self.findAlignment()

    def findAlignment(self, i=None, j=None,seq = "",seq2 = "",seq3 = ""):
        if i is None and j is None:
            i = self.S.n - 1
            j = self.S.m - 1
        if i >= 0 and j >= 0 and self.solutionCounter < self.k :
            if i == 0 and j == 0:
                self.solutionCounter += 1
                self.similarityRate = round((seq3.count("|") + seq3.count(":")) / len(seq3) * 100,2)
            else:

                if self.S.getItem(i, j) == self.S.getItem(i - 1, j - 1) + self.t(i, j):
                    seq += self.S.sequence1[i - 1]
                    seq2 += self.S.sequence2[j - 1]
                    if self.S.sequence1[i - 1] == self.S.sequence2[j - 1]:
                        seq3 += "|"
                    else:
                        seq3 += ":" if self.t(i, j) >= 0 else "."
                    self.findAlignment(i - 1, j - 1,seq,seq2,seq3)
                    seq = seq[:-1]
                    seq2 = seq2[:-1]
                    seq3 = seq3[:-1]

                if self.S.getItem(i, j) == self.V.getItem(i, j) or j == 0:
                    seq += self.S.sequence1[i - 1]
                    seq2 += "-"
                    seq3 += " "
                    self.findAlignment(i - 1, j,seq,seq2,seq3)
                    seq = seq[:-1]
                    seq2 = seq2[:-1]
                    seq3 = seq3[:-1]

                if self.S.getItem(i, j) == self.W.getItem(i, j) or i == 0:
                    seq += "-"
                    seq2 += self.S.sequence2[j - 1]
                    seq3 += " "
                    self.findAlignment(i, j - 1,seq,seq2,seq3)
                    seq = seq[:-1]
                    seq2 = seq2[:-1]
                    seq3 = seq3[:-1]

    def t(self, i, j):
        x = self.MATSUB[self.S.sequence2[j - 1]]
        y = self.MATSUB[self.S.sequence1[i - 1]]
        return int(self.MATSUB.getItem(x, y))

    def fill(self, i, j):

        self.V.setItem(i, j, max(self.S.getItem(i - 1, j) - self.I, self.V.getItem(i - 1, j) - self.E))
        self.W.setItem(i, j, max(self.S.getItem(i, j - 1) - self.I, self.W.getItem(i, j - 1) - self.E))
        self.S.setItem(i, j, max(self.S.getItem(i - 1, j - 1) + self.t(i, j), self.V.getItem(i,j), self.W.getItem(i,j)))

    def getSimilarityRate(self) :
        return self.similarityRate

def findSequeneces():
    MATSUB = matrixParser("blosum80.txt")
    I = 6
    E = 1
    sequencesName, sequences = sequencesParser2("to-be-aligned.fasta")
    n = len(sequences)
    sequencesList = [i for i in range(n)]
    sequencesSetList = []
    for i in range(n) :
        print(i)
        sequenceNumber = choice(sequencesList)
        j = 0
        inserted = False
        while not inserted and j < len(sequencesSetList):
            inserted = canBeInsert(sequenceNumber,sequencesSetList[j],sequences,I,E,MATSUB)
            if inserted :
                sequencesSetList[j].append(sequenceNumber)
            j += 1
        if not inserted:
            sequencesSetList.append([sequenceNumber])
        sequencesList.remove(sequenceNumber)

    maxLen = [(len(elem),elem) for elem in sequencesSetList]
    maxLen.sort(reverse=True)
    solution = maxLen[0][1]
    generateReducedFile("to-be-aligned-reduced.fasta",solution,sequencesName,sequences)

def canBeInsert(sequenceNumber,set,sequences,I,E,MATSUB):
    seq1 = ADTsequence(sequences[sequenceNumber])
    for j in set:
        seq2 = ADTsequence(sequences[j])
        S = ScoringMatrix("S","GLOBAL",seq1,seq2,I,E)
        V = ScoringMatrix("V","GLOBAL",seq1,seq2)
        W = ScoringMatrix("W","GLOBAL",seq1,seq2)
        if NeedlemanWunsch(S,V,W,MATSUB,1,I,E).getSimilarityRate() > 60:
            return False
    return True

def sequencesParser2(filename):
	f = open(filename, encoding='utf-8')
	sequencesName = []
	sequencesList = []

	for line in f:
		if line[0] == ">":
			sequencesName.append(line.strip("\n"))
		elif line[0].isalnum() :
			sequencesList.append(line.strip("\n"))
	f.close()
	return (sequencesName,sequencesList)

def generateReducedFile(filename,reducedListIndex,sequencesName,sequencesList):
    f = open(filename, "w", encoding='utf-8')

    for index in reducedListIndex:
        sequence = sequencesName[index] + "\n" + sequencesList[index] + "\n"
        print(sequencesList[index])
        f.write(sequence)
    f.close

class PositionSpecificScoringMatrix(ADTmatrix):
    def __init__(self,sequences,I,n=20, m=None):
        self.l = n
        if m == None:
            self.k = len(sequences[0])
        else:
            self.k = m

        super().__init__(self.l, self.k)
        self.Nseq = len(sequences)
        #self.aminoAcidList = ["A","Q","L","S","R","E","K","T","N","G","M","W","D","H","F","Y","C","I","P","V","-"]
        self.aminoAcidList = "ACDEFGHIKLMNPQRSTVWY-"
        self.dico = {}
        self.consensus = [(0,0) for i in range(self.k)]
        self.createDico(sequences)
        self.p =  { "A" : 8.28, "Q" : 3.94, "L" : 9.67, "S" : 6.50,
                    "R" : 5.53, "E" : 6.76, "K" : 5.85, "T" : 5.32,
                    "N" : 4.05, "G" : 7.09, "M" : 2.43, "W" : 1.07,
                    "D" : 5.45, "H" : 2.27, "F" : 3.86, "Y" : 2.91,
                    "C" : 1.36, "I" : 5.99, "P" : 4.68, "V" : 6.87}

        self.initPSSM()
        self.matrix.append([-I for n in range(self.k)])

    def initPSSM(self):
        for i in range(self.l) :
            for j in range(self.k):
                self.setItem(i,j,self.m(j,self.aminoAcidList[i]))

        self.findConsensus()

    def m(self,u,a):
        return round(log(self.q(u,a)/(self.p[a]/100),10),3)

    def q(self,u,a,beta=1):
        return (self.n(u,a) + beta * (self.p[a]/100)) / (self.Nseq + beta)

    def n(self,u,a):
        return self.dico[(u,a)]

    def getItemValue(self, char, j):
        return self.matrix[self.aminoAcidList.index(char)][j]

    def createDico(self,sequences):
        for aminoAcid in self.aminoAcidList:
            for pos in range(self.k):
                self.dico[(pos,aminoAcid)] = 0
                for sequence in sequences:
                    if sequence[pos] == aminoAcid:
                        self.dico[(pos,aminoAcid)] += 1

    def findConsensus(self) :
        for i in range(self.k):
            for j in range(self.l):
                value = self.getItem(j,i)
                if value > self.consensus[i][0]:
                    self.consensus[i] = (value,self.aminoAcidList[j])

class ScoringMatrix(ADTmatrix):

    def __init__(self, sequence, profil, I=None, E=None):
        super().__init__(len(sequence) + 1, len(profil[0]) + 1)
        self.profil = profil
        self.sequence = sequence
        self.m = len(profil[0]) + 1
        self.n = len(sequence) + 1
        self.initS(I, E)

    def initS(self, I, E):

        for i in range(1, self.m):
            gap = - I - (i - 1) * E
            self.setItem(0, i, max(0, gap))


        for j in range(1, self.n):
            gap = - I - (j - 1) * E
            self.setItem(j, 0, max(0, gap))

class SmithWaterman(object):
    def __init__(self, S, PSSM, l, I, E):

        self.S = S
        self.PSSM = PSSM
        self.I = I
        self.E = E

        for i in range(1, self.S.n):
            for j in range(1, self.S.m):
                self.fill(i, j)
        i = 0
        while i < l and sum([sum(row) for row in self.S.matrix]):
            indexList = []
            self.topDown(i+1,indexList)
            self.fillZero(indexList)
            x, y = indexList[-1]
            self.computeAgain(x, y)
            i += 1

    def fill(self, i, j):
        case1 = self.S.getItem(i - 1, j - 1) + self.PSSM.getItemValue(self.S.sequence[i-1], j-1)
        case2 = self.S.getItem(i - 1, j) + self.PSSM.getItemValue("-", j-1)
        case3 = self.S.getItem(i, j - 1) + self.PSSM.getItemValue("-", j-2)
        self.S.setItem(i, j, round(max(0, case1, case2, case3),3))

    def maxValue(self):

        x, y = 1, 1
        maximum = self.S.getItem(x, y)
        for i in range(1, self.S.n):
            for j in range(1, self.S.m):
                if maximum < self.S.getItem(i, j):
                    maximum = self.S.getItem(i, j)
                    x, y = i, j
        return x, y

    def topDown(self, number, indexList):

        i, j = self.maxValue()
        indexList.append((i, j))
        score = self.S.getItem(i, j)
        seq = ""
        while self.S.getItem(i, j):

            value = self.S.getItem(i, j)
            case1 = self.S.getItem(i - 1, j - 1) + self.PSSM.getItemValue(self.S.sequence[i-1], j-1)
            case2 = self.S.getItem(i - 1, j) + self.PSSM.getItemValue("-", j-1)
            case3 = self.S.getItem(i, j - 1) + self.PSSM.getItemValue("-", j-2)
            startPos = i

            if value == round(case1,3):
                seq += self.S.sequence[i - 1]
                i, j = i - 1, j - 1

            elif value == round(case2,3):
                seq += "-"
                i = i - 1

            elif value == round(case3,3):
                seq += "-"
                j = j - 1

            indexList.append((i, j))

        show(number, score, seq[::-1], startPos)
        indexList.pop()

    def fillZero(self, indexList):

        for elem in indexList:
            self.S.setItem(elem[0], elem[1], 0)

    def computeAgain(self, x, y):
        i, j = x, y
        while i < self.S.n:
            while j < self.S.m:
                if self.S.getItem(i,j) > 0:
                    self.fill(i, j)
                j += 1
            i += 1
            j = y

def show(number, score, sequence, start):
    print("Solution n°: {0} ".format(number))
    w = (len(sequence) // 60 + 1) if len(sequence) % 60 != 0 else len(sequence) // 60
    from_ = start
    for i in range(w):
        taille = len(sequence[i * 60:(i + 1) * 60]) - sequence[i * 60:(i + 1) * 60].count("-") - 1
        to_  = from_ + taille
        print(sequence[i * 60:(i + 1) * 60] + " - from {0} to {1}".format(from_,to_) + "\n")
        from_ = to_

    print("Score :", score)
    print("-" * 60)

if __name__ == "__main__":
    from math import log
    from random import choice

    I = E = 2
    l = 1
    profil = sequencesParser("msaresults-ClustalOmega.fasta")
    sequences = sequencesParser("protein-sequences.fasta")
    PSSM = PositionSpecificScoringMatrix(profil,I)
    for i in range(len(sequences)):
        print("Alignment between sequence {0} and profil \n".format(i+1))
        sequence = ADTsequence(sequences[i])
        S = ScoringMatrix(sequence,profil,I,E)
        SmithWaterman(S, PSSM, l, I, E)
    """
    findSequeneces()
    """
