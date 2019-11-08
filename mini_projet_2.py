def matrixParser(filename):
    f = open(filename, encoding='utf-8')
    matrix = SubstitutionMatrix()
    for line in f:
        if line[0] not in ["#", " ", "\n"]:
            matrix.Append(line[0])
            matrix.addRow(line[1:].split())
    f.close()
    return matrix

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
        output = ""
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix[0])):
                output += str(self.matrix[i][j])
                output += " "
            output += "\n"
        return output[:-1]


class SubstitutionMatrix(ADTmatrix):

    def __init__(self, n=None, m=None):
        super().__init__(n, m)
        self.myList = []

    def __getitem__(self, key):
        return self.myList.index(key)

    def Append(self, value):
        self.myList.append(value)


class ScoringMatrix(ADTmatrix):

    def __init__(self, character, mode, sequence1, sequence2, I=None, E=None):
        super().__init__(len(sequence1) + 1, len(sequence2) + 1)
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.n = len(sequence1) + 1
        self.m = len(sequence2) + 1

        if character == "V":
            self.initV()
        elif character == "W":
            self.initW()
        else:
            self.initS(mode, I, E)

    def initS(self, mode, I, E):

        for i in range(1, self.m):
            gap = - I - (i - 1) * E
            if mode == "LOCAL":
                self.setItem(0, i, max(0, gap))
            else:
                self.setItem(0, i, gap)

        for j in range(1, self.n):
            gap = - I - (j - 1) * E
            if mode == "LOCAL":
                self.setItem(j, 0, max(0, gap))
            else:
                self.setItem(j, 0, gap)

    def initV(self):
        for i in range(self.m):
            self.setItem(0, i, - float("inf"))

        for j in range(1, self.n):
            self.setItem(j, 0, 0)

    def initW(self):
        for i in range(1, self.m):
            self.setItem(0, i, 0)

        for j in range(self.n):
            self.setItem(j, 0, - float("inf"))


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


class SmithWaterman(object):
    def __init__(self, S, V, W, MATSUB, l, I, E):

        self.S = S
        self.V = V
        self.W = W
        self.MATSUB = MATSUB
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

    def t(self, i, j):
        x = self.MATSUB[self.S.sequence2[j - 1]]
        y = self.MATSUB[self.S.sequence1[i - 1]]
        return int(self.MATSUB.getItem(x, y))

    def fill(self, i, j):

        self.V.setItem(i, j, max(0, self.S.getItem(i - 1, j) - self.I, self.V.getItem(i - 1, j) - self.E))
        self.W.setItem(i, j, max(0, self.S.getItem(i, j - 1) - self.I, self.W.getItem(i, j - 1) - self.E))
        self.S.setItem(i, j, max(0, self.S.getItem(i - 1, j - 1) + self.t(i, j), self.V.getItem(i,j), self.W.getItem(i,j)))

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
        seq2 = ""
        seq3 = ""
        pos = (i,j)
        while self.S.getItem(i, j):

            value = self.S.getItem(i, j)
            pos2 = (i,j)

            if value == self.S.getItem(i - 1, j - 1) + self.t(i, j):
                seq += self.S.sequence1[i - 1]
                seq2 += self.S.sequence2[j - 1]
                if self.S.sequence1[i - 1] == self.S.sequence2[j - 1]:
                    seq3 += "|"
                else:
                    seq3 += ":" if self.t(i, j) >= 0 else "."
                i, j = i - 1, j - 1

            elif value == self.V.getItem(i, j):
                seq += self.S.sequence1[i - 1]
                seq2 += "-"
                seq3 += " "
                i = i - 1

            elif value == self.W.getItem(i, j):
                seq += "-"
                seq2 += self.S.sequence2[j - 1]
                seq3 += " "
                j = j - 1

            indexList.append((i, j))

        show(number, score, seq[::-1], seq2[::-1], seq3[::-1],pos,pos2)
        indexList.pop()

    def fillZero(self, indexList):

        for elem in indexList:
            self.S.setItem(elem[0], elem[1], 0)
            self.V.setItem(elem[0], elem[1], 0)
            self.W.setItem(elem[0], elem[1], 0)

    def computeAgain(self, x, y):
        i, j = x, y
        while i < self.S.n:
            while j < self.S.m:
                if self.S.getItem(i,j) > 0 or self.V.getItem(i,j) > 0 or self.W.getItem(i,j) > 0:
                    self.fill(i, j)
                j += 1
            i += 1
            j = y

def show(number, score, sequence, sequence2, sequence3, pos=None, pos2=None):
    print("Solution n°: {0} ".format(number))
    w = (len(sequence) // 60 + 1) if len(sequence) % 60 != 0 else len(sequence) // 60
    if pos is not None and pos2 is not None :
        from_ = pos2[0] - 1
        from__ = pos2[1] - 1
    for i in range(w):
        if pos is None and pos2 is None :
            print(sequence[i * 60:(i + 1) * 60] \
                + "\n" + sequence3[i * 60:(i + 1) * 60] \
                + "\n" + sequence2[i * 60:(i + 1) * 60] + "\n")
        else :
            taille = len(sequence[i * 60:(i + 1) * 60]) - sequence[i * 60:(i + 1) * 60].count("-")
            taille_ = len(sequence2[i * 60:(i + 1) * 60]) - sequence2[i * 60:(i + 1) * 60].count("-")
            to_  = from_ + taille
            to__  = from__ + taille_
            print(sequence[i * 60:(i + 1) * 60] \
                + " - from {0} to {1}".format(from_,to_) \
                + "\n" + sequence3[i * 60:(i + 1) * 60] \
                + "\n" + sequence2[i * 60:(i + 1) * 60] \
                + " - from {0} to {1}".format(from__,to__) + "\n")

            from_ = to_
            from__ = to__

    #print("Score :", score)
    #print("Identity percentage :", round(sequence3.count("|") / len(sequence3) * 100,2), "%")
    #print("Similarity rate:", round((sequence3.count("|") + sequence3.count(":")) / len(sequence3) * 100,2), "%")
    idt = round(sequence3.count("|") / len(sequence3) * 100,2)
    sim = round((sequence3.count("|") + sequence3.count(":")) / len(sequence3) * 100,2)
    print("Scr =  {0} | Idt =  {1}% | Sim =  {2}%".format(score,idt,sim))
    print("-" * 60)

def subset(similarityRateMatrix):
	lst = [set() for i in range(len(similarityRateMatrix))]
	for i in range(len(similarityRateMatrix)):
		for j in range(len(similarityRateMatrix)):
			if i !=j and similarityRateMatrix[j][i] <= 60:
				temp = lst[i] | {i,j}
				if generateSubSetOf2(list(temp), similarityRateMatrix):
					lst[i] |= {i,j}

	sortingSet = [(len(lst[i]),i) for i in range(len(lst))]
	sortingSet.sort(reverse=True)
	found = False
	i = 0
	while not found and  i < len(lst):
		set_ = lst[sortingSet[i][1]]
		found = generateSubSetOf2(list(set_), similarityRateMatrix)
		i +=1
	return set_

def generateSubSetOf2(lst, similarityRateMatrix, k=2, setOfK=[], i=0,condition = True):
	if condition :
		if k < 0:
			return condition
		if k == 0:
			return True if similarityRateMatrix[setOfK[0]][setOfK[1]] <= 60 else False
		elif i + k <= len(lst):

			setOfK.append(lst[i])
			condition = generateSubSetOf2(lst,similarityRateMatrix,k-1,setOfK,i+1,condition)

			setOfK.pop()
			condition = generateSubSetOf2(lst,similarityRateMatrix,k,setOfK,i+1,condition)
	return condition

def sequencesParser(filename):
	f = open(filename, encoding='utf-8')
	sequencesName = []
	sequencesList = []

	for line in f:
		if line[0] == ">":
			sequencesName.append(line.strip("\n"))
		else :
			sequencesList.append(line.strip("\n"))
	f.close()
	return (sequencesName,sequencesList)

if __name__ == "__main__":

	MATSUB = matrixParser("pam120.txt")
	sequencesName, sequences = sequencesParser("to-be-aligned.fasta")
	n = len(sequences)
	similarityRateMatrix = [[0 for i in range(n)] for j in range(n)]

	for i in range(len(sequences)):
		for j in range(i+1,len(sequences)) :
			seq1 = ADTsequence(sequences[i])
			seq2 = ADTsequence(sequences[j])
			S = ScoringMatrix("S","GLOBAL",seq1,seq2,12,2)
			V = ScoringMatrix("V","GLOBAL",seq1,seq2)
			W = ScoringMatrix("W","GLOBAL",seq1,seq2)
			similarityRateMatrix[i][j] = similarityRateMatrix[j][i] = NeedlemanWunsch(S,V,W,MATSUB,1,12,2).getSimilarityRate()
	solution = subset(similarityRateMatrix)
	print(solution)
