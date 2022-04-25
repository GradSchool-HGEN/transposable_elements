#!/bin/python
# This script will parse the text file output of nhmmer. 
#
#
#
# Abin Abraham
# created on: 2017-12-24 15:16:04



# store data as a nested dictionary 
    # outer key = one sequecne 
        # model name: <model>
        # model start: ,start.
        # model end : <end> 
        # model seq : <seq> 
        # match bases: <seq> 
        # elem seq Name: 
        # elem start
        # elem end 
        # elem seq
        # elem base score  


# -----------
# functions
# ----------- 

def checkCharLength(inputStr, reqCharLength, lineNumber):
    if len(inputStr) != reqCharLength:
        print("Number of character in line " + str(lineNumber) + "in file does not match what was expected")

def parseTermOut(fileName):
    storeAlign = dict()  
    recordCounter = 0
    recordFlag = False
    duplicateFlag = False
    with open(fileName, 'r') as infile:
        for num,line in enumerate(infile):
        
            if line.find('>>') != -1: 
                line = line.split()
                thisSeq = line[1]
                recordFlag = True
                recordCounter = 1
                continue

            if recordFlag: 
                recordCounter += 1

                if recordCounter ==4: 
                    if line.split()[0] == "!":
                        statSignif_FLAG = True
                    else: 
                        statSignif_FLAG = False
                    
                if recordCounter ==8: 
                    numBufferChar = line.find('x')
                    numChar = len(line.split()[0])

                if recordCounter == 9: 
                    line = line.split()
                    model_name = line[0]
                    model_start = line[1]
                    model_end = line[3]
                    model_seq = line[2]
                    checkCharLength(model_seq, numChar, num)

                if recordCounter == 10: 
                    templine = line.split()[0]
                    matchedBases = line[numBufferChar:numBufferChar+numChar]

                if recordCounter == 11: 
                    line = line.split()
                    elemSeqName = line[0]
                    elemSeqStart = line[1]
                    elemSeqEnd = line[3]
                    elemSeq = line[2]
                    checkCharLength(elemSeq, numChar, num)

                if recordCounter == 12: 
                    line = line.split()
                    elemBaseScore = line[0]
                    checkCharLength(elemBaseScore, numChar, num)
                    recordFlag = False
                    recordCounter = 0
    
                    if statSignif_FLAG: 
                        if thisSeq in storeAlign: 
                            duplicateFlag = True
                            # print("duplicate key found for: " + thisSeq)
                        else: 
                            storeAlign[thisSeq] = {'model_name':model_name, 'model_start':model_start, 'model_end':model_end, 'model_seq':model_seq,
                                'matchedBases':matchedBases, 'elemSeqName':elemSeqName, 'elemSeqStart':elemSeqStart,
                                'elemSeqEnd':elemSeqEnd, 'elemSeq':elemSeq, 'elemBaseScore':elemBaseScore}

    if duplicateFlag: 
        print("\n")
        print("While parsing " +fileName+ " there were duplicates that were not included in analysis")
        print("\n")
        
    return storeAlign

def checkParseFile(parsedOutput):
    for seq in parsedOutput: 

        tempstore = parsedOutput[seq]['model_seq'].replace('-','')
        tempstore = tempstore.replace('.','')
        if len(tempstore) != ((int(parsedOutput[seq]['model_end']) - int(parsedOutput[seq]['model_start']))+1):
            print("Model length does not match for: "+ seq)
            
        tempstore = parsedOutput[seq]['elemSeq'].replace('-','')
        tempstore = tempstore.replace('.','')
        if len(tempstore) != ((int(parsedOutput[seq]['elemSeqEnd']) - int(parsedOutput[seq]['elemSeqStart']))+1):
            print("Element length does not match given start and end coordinatesfor: "+ seq)
            print(tempstore)
            print(len(tempstore))
            print(int(parsedOutput[seq]['elemSeqEnd']))
            print(int(parsedOutput[seq]['elemSeqStart']))
            print(((int(parsedOutput[seq]['elemSeqEnd']) - int(parsedOutput[seq]['elemSeqStart']))+1))
        
        if len(parsedOutput[seq]['model_seq']) == len(parsedOutput[seq]['elemSeq']) == len(parsedOutput[seq]['elemBaseScore']) == len(parsedOutput[seq]['matchedBases']):
            print("Checking length of alignments all match!")
        else: 
            print("For " +seq+ " length of alignments read in do not equal each other:")
            print("Number of char in model sequence: " + str(len(parsedOutput[seq]['model_seq'])))
            print("Number of char in element sequence: " + str(len(parsedOutput[seq]['elemSeq'])))
            print("Number of char in Element Base Score: " + str(len(parsedOutput[seq]['elemBaseScore'])))
            print("Number of char in Matched Bases sequence: " + str(len(parsedOutput[seq]['matchedBases'])))
            print(parsedOutput[seq]['matchedBases'])

# -----------
# main
# ----------- 
if __name__ == '__main__': 
    fileName = "nhmmer-output-terminal_allMER20"
    print("Runnning default settings....")
    print(fileName)
    # fileName = "nhmmer-output-terminal"
    output = parseTermOut(fileName)
    # print(output)
    checkParseFile(output)

                