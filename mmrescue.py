##4.26.2021
##Lexi Morrissey
##Re-implementation of MuMRescueLite

#Used to get counts in regions around multimapped reads
def genSearch(read,winSize,genLand):
    count = 0
    for k in range (int(read[3])-winSize,int(read[3])+winSize):
        key = read[2] + ";" + str(k)
        #Seeing if current pos in genetic landscape
        if key in genLand:
            count = count + genLand[key]
    return count

                    
def readAssign(rBlock, genLand, samOut):
    import random
    ##Uniquely mapped reads##
    if len(rBlock) == 1:
        #Adding to file
        samOut.write('\t'.join(rBlock[0]))
        #Adding to genetic landscape
        key = rBlock[0][2] + ";" + str(int(rBlock[0][3]))
        if key in genLand:
            genLand[key] = genLand[key] + 1
        else:
            genLand[key] = 1
        return genLand
    
    ##Multi-mapped reads##
    #Getting read counts in regions around multimapped reads
    counts = []
    for i in rBlock:
        counts.append(genSearch(i,winSize,genLand))

    #Not rescuing reads that do not have any reads in the window around it
    totalCount = sum(counts)
    if totalCount == 0:
        return genLand

    for i in range(0,len(counts)):
      counts[i] = round(counts[i]/totalCount*100)
    #Choosing which read to keep based on probabilities from nearby read counts
    found = 0
    x = 1       #Used for starting point to make sure percentages add to 100
    randN = random.randint(1,sum(counts))       #Used to randomly allocate reads
    for i in range(0, len(rBlock)):
        if randN in range(x,x+counts[i]):
            #Adding to file
            found = 1
            samOut.write('*'+'\t'.join(rBlock[i]))
            return genLand
        x = x + counts[i]
    return genLand


def parseUniq(tempFile):
    import os
    
    genLand = {}
    rBlock = []
    #File with uniquely mapped reads
    UM = open(tempFile + "UM","w+")
    #File with unallocated multimapped reads
    MM = open(tempFile + "MM","w+")
    
    with open(tempFile) as f:
        for line in f:

            #Exception
            if not line.strip():
                break
            
            #First file will have sam header in it
            if line.startswith("@"):
                continue

            #Splitting columns
            r = line.split('\t')

            #Some seq exceptions
            if r[2] == "*":
                continue
            if "N" in r[9]:
                continue

            #Appending if the block is empty and going to next line
            if len(rBlock) == 0:
                rBlock.append(r)
                continue
            else:
                #Adding read to block if alignment score is equal to last
                leadr = rBlock[0]
                rScore = int(r[11].split(':')[2])
                leadScore = int(leadr[11].split(':')[2])
                
                if r[0] == leadr[0] and rScore == leadScore:
                    rBlock.append(r)
                    continue
                
                #Deleting old block if read found with better score
                if r[0] == leadr[0] and rScore > leadScore:
                    rBlock = []
                    rBlock.append(r)
                    continue
                
                #Allocating reads to genomic landscape for previous block before
                #starting new one
                if not r[0] == leadr[0]:
                    #Put uniquely mapped reads into a file
                    if len(rBlock) == 1:
                        genLand = readAssign(rBlock, genLand, UM)
                    #Put multimapped reads into a file
                    if len(rBlock) > 1:
                        for i in range(0,len(rBlock)):
                            MM.write('\t'.join(rBlock[i]))
                    #Creating a new read block for the next read
                    rBlock = []
                    rBlock.append(r)

    #Running function one more time for last read
    if len(rBlock) == 1:
        genLand = readAssign(rBlock, genLand, UM)
    #Put multimapped reads into a file
    if len(rBlock) > 1:
        for i in range(0,len(rBlock)):
            MM.write('\t'.join(rBlock[i]))

    UM.close()
    MM.close()
    os.remove(tempFile) #Removing old temp
    return genLand


def parseMulti(tempFile, genLand):
    import os
    if os.stat(tempFile+"MM").st_size == 0:
        return
    #File with MM allocations
    AL = open(tempFile + "AL","w+")
    tempFile = tempFile + "MM"
    rBlock = []
    with open(tempFile) as f:
        for line in f:

            #Exception
            if not line.strip():
                break
            
            #Splitting columns
            r = line.split('\t')
            
            #Some seq exceptions
            if r[2] == "*":
                continue
            if "N" in r[9]:
                continue
            #Appending if the block is empty and going to next line
            if len(rBlock) == 0:
                rBlock.append(r)
                continue
            else:
                #Adding line to block if same read
                leadr = rBlock[0]
                rScore = int(r[11].split(':')[2])
                leadScore = int(leadr[11].split(':')[2])
                
                if r[0] == leadr[0]:
                    rBlock.append(r)
                    continue
                
                #Allocating reads to genomic landscape for previous block before
                #starting new one
                if not r[0] == leadr[0]:
                    genLand = readAssign(rBlock, genLand, AL)
                    #Creating a new read block for the next read
                    rBlock = []
                    rBlock.append(r)

    #For last read
    genLand = readAssign(rBlock, genLand, AL)

    os.remove(tempFile) #Removing old temp
    AL.close()



##Main Method##          
if __name__ == '__main__':
    
    import sys
    from joblib import Parallel, delayed
    import multiprocessing
    import subprocess
    import os
    import math

    
    #Taking in arguments
    samfile = sys.argv[1]
    outfile = sys.argv[2]
    if "-w" in sys.argv:
        winSize = int(sys.argv[sys.argv.index("-w")+1])
    else:
        winSize = 200
    if "-t" in sys.argv:
        thr = int(sys.argv[sys.argv.index("-t")+1])
    else:
        thr = multiprocessing.cpu_count()
    

    #Making temporary smaller files based on the number of threads
    file = str(sys.argv[1])
    
    s = subprocess.Popen(['wc','-l', file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    size = str(s.communicate()[0])
    size = str(size)
    size = size[2:size.find(" ")]
    size = str(math.ceil(int(size)/int(thr)))
    
    
    p = subprocess.Popen(['split', '-d', '-l', size, file, 'temp'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()

   
    tempList = []
    #Need list of file names
    for i in range (0,thr):
        tempList.append("temp" + str(i).zfill(2))


    output = Parallel(n_jobs=thr, prefer="threads")(delayed(parseUniq)(i) for i in tempList)
    genLand = {}
    for d in output:
        genLand.update(d)

    output = []
    output = Parallel(n_jobs=thr, prefer="threads")(delayed(parseMulti)(i, genLand) for i in tempList)

    cmd = []
    cmd.append('cat')
    cmd.append('tempHeader')
    for i in tempList:
        if os.path.exists(str(i)+"AL"):
            cmd.append(i + "AL")
        cmd.append(i + "UM")
    cmd.append('>')
    cmd.append(outfile)
    cmd = ' '.join(cmd)
    os.system(cmd)

    cmd = []
    cmd.append('rm')
    for i in tempList:
        if os.path.exists(str(i)+"AL"):
            cmd.append(i + "AL")
        cmd.append(i + "UM")
    cmd = ' '.join(cmd)
    
    os.system(cmd)
    os.system("grep \"@\" " + samfile+ " > tempHead")
    os.system("cat tempHead " + outfile + " >" + outfile + "_temp")
    os.system("mv " + outfile + "_temp " + outfile)
    os.system("rm temp*")
    