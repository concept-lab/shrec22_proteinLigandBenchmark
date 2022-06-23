#       Evaluation script of Shrec2022 protein-ligand detection contest
#       Metrics and ranking are described in the associated publication

#         **    MIT LICENSE    **
#       Copyright (C) 2022  Luca Gagliardi
#       Affiliation: Istituto Italiano di Tecnologia
    
    
#           +DESCRIPTION+
#
#   0. Read from teminal the folder name
#   1. Automatically detects format (OFF or PQR). For OFF the folder testData with all structure triangulations must be present in the working directory.
#   2. Map each participant pocket file to the correspondent ligand (for putative site evaluation) using the provided file "testMap.txt"
#   3. Load ligand atoms from the folder containing ligand files in xyz format
#   4. Load putative pockets of the participant
#   5. Establish success according to LC and PC metrics



import numpy as np
import re
from  C_functs import Pdist_C,getIndex
from time import time, strftime, localtime

thresholdForTRIANG = 4.
thresholdForPQR = 5.

pCoverageTH = 0.2
lCoverageTH = 0.5

print("LC threshold = ",lCoverageTH)
print("PC threshold = ",pCoverageTH)

def secondsToStr(elapsed=0):
    if elapsed ==0 :
        return strftime("LOCAL TIME = %Y-%m-%d %H:%M:%S", localtime())
    else:
        return "ELAPSED TIME: "+ str(elapsed)#str(timedelta(seconds=elapsed))

class Crono(object):
    def init(self):
        self._elapsed=0
        self._start=time()
        return secondsToStr(self._elapsed)
    def get(self):
        end=time()
        self._elapsed=end-self._start
        return secondsToStr(self._elapsed)


def get_protein(structure):
    '''
    Input: PQR file
    Output: List of dictionary. Each entry is a line in the PQR file.
    '''
    try:
        # print(structure+'.pdb')
        inFile = open(structure+'.pqr','r')
    except Exception:
        raise NameError("Cannot load PQR file")
    # try:
    #     # print(structure+'.pdb')
    #     _check = open(structure+'.pdb','r')
    # except Exception:
    #     raise NameError("Cannot load PDB file")
    comment =['#', 'CRYST[0-9]?']
    remark = ['REMARK']
    termination = ['TER', 'END', '\n']
    skip = comment+remark+termination
    skip = '(?:% s)' % '|'.join(skip)
    for line in inFile: 
        if(re.match(skip,line)): 
            pass 
        else:
            linegNOChain=re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            linegChain = re.match("(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)",line)
            break

    if(linegChain):
        isChainID=1 
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+([\w0-9]+)\s*(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
    elif(linegNOChain):
        print('no chain')
        matchPattern = "(ATOM)\s*(\d+)\s*(\S+)\s+([A-Z]+)\s+(\-?\d+[A-Z]?)\s+(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s*(\-?\d*\.?\d+)\s+(\-?\d*\.?\d+)\s+(\d*\.?\d+)"
        isChainID =0 
    else:
        raise NameError("Incorrect pqr file formatting")
    if(isChainID):
        resInd = 5
        chargeInd = 9
        rInd = 10
    else:
        resInd = 4
        chargeInd = 8
        rInd = 9
    nameInd = 3
    atomInd = 2
    coordInd = resInd +1
    chainInd=4
    
    inFile.seek(0)

    resMap=[]
    for line in inFile:
        if(re.match(skip,line)): 
            pass 
        else: 
            # mline=line.split() Not robust when the field number and chain touch..
            mline=re.match(matchPattern,line).groups()
            if(isChainID):
                # try:
                content = {'resName':mline[nameInd],'resNum':mline[resInd],'atomNumber':mline[1],'resAtom': mline[atomInd],'resChain':mline[chainInd],
                'charge':float(mline[chargeInd]),'coord':list(map(float, mline[coordInd:coordInd+3])),'radius':float(mline[rInd])}
                # except:
                #     print(line)
                #     print(mline)
                #     exit()
            else:
                content = {'resName':mline[nameInd],'resNum':mline[resInd],'atomNumber':mline[1],'resAtom': mline[atomInd],
                'charge':float(mline[chargeInd]),'coord':list(map(float, mline[coordInd:coordInd+3])),'radius':float(mline[rInd])}
            
            resMap.append(content)
    
    return resMap

def getStructureIterator(mapFile):
        
    comment = ['^#','\n']
    comment = '(?:% s)' % '|'.join(comment)
    # mapFile = 'ligandMap.txt'
    structures=[]
    try:
        inFile = open(mapFile,'r')
    except Exception:
        raise NameError("Cannot load mapFile")

    content = inFile.readlines()

    structures = {}
    s=0
    while s <(len(content)):
        current_line=content[s]
        if(re.match(comment,current_line)):
            s+=1
            continue
        current_line = current_line.split()
        n_ligands = int(current_line[0])
        name = current_line[1]
        folderNumber = [int(s) for s in current_line[2:] if s.isdigit()][0]
        ligand_names=[]
        for i in range(1,n_ligands+1):
            following_line=content[s+i].split()[0]
            if(re.match(comment,following_line)):
                #skipping ligand
                continue
            ligand_names.append(following_line)
        structures[folderNumber]={'pqr':name,'ligands': ligand_names}
        s+=n_ligands +1
            
    return structures



def readLigand(lname,structureName,filter=True):
    '''
    Return ligand atoms read from lname within 5 angstrom from protein atoms if filter on.
    For that function one needs to provide the whole protein pqr..
    In alternative, use the filtered ligand folder
    '''
    if filter:
        resMap = get_protein('allStructures/'+structureName)
        proteinCoord = np.empty((0,4))
        for i in resMap:
            c = np.append(np.asarray(i['coord']),i['radius'])
            proteinCoord = np.vstack([proteinCoord,c])
        inFile = open('ligands/'+lname+".xyz",'r')
        ligand_coord = np.loadtxt(inFile)

        d,flag = Pdist_C(ligand_coord[:,0:3],proteinCoord[:,0:3])
        index = np.where(d<=5)[0]
        lindex,_pindex=getIndex(flag,index,proteinCoord.shape[0])
        # print(lindex)
        # print(np.unique(lindex))
        # print(ligand_coord[np.unique(lindex)])
        ligand_coord = ligand_coord[np.unique(lindex)]
    else:
        inFile = open('filtered_ligands/'+lname+".xyz",'r')
        ligand_coord = np.loadtxt(inFile)

    return ligand_coord


def closeVert(ligCoord,vert,thresholdDistance):
    """
    Input ligand coordinates and single vertex.
    Returns true if any of of vertices closer or equal than threshold to ligand atom
    Returns indexex of the ligand close to that vertex
    Ligand heavy atoms pre-obtained by MOAD_ligand finder script..
    """
    # not sure passing through C is interesting here
    d,_flag = Pdist_C(ligCoord,vert)
    index = np.where(d<=thresholdDistance)[0]
    lindex,_pindex=getIndex(_flag,index,vert.shape[0])
    if (index.size > 0):
        close = True
    else:
        close = False
    
    # print(close)
    # print(index)
    return close,lindex

def loadTriang(triangName):   
    try:
        triangulationFile = open(triangName,'r')
    except FileNotFoundError:
        print('Cannot find triangulation file (.off).\nBREAKING')
        exit()
    triangulationFile.readline()
    triangulationFile.readline()
    triangulationFile.readline()
    infoLine=triangulationFile.readline()
    nverts=int(infoLine.split()[0])
    print('number of vertices in the triangulation =',nverts)
    triangLines = triangulationFile.readlines()
    
    return triangLines[0:nverts]


#LEGGI POCKET BOOLEAN MAP CHE SIA VERTICE O PQR E RITORNA LINEE CORRISPONDENTI DA FILES ORIGINALI (o pqr è gia riempito?..)


def getClose(ligandCoord,pocket,isPQR):
    '''
    Input a given ligand coordinate and a given pocket lines. Different pockets and ligands must be read externally.
    Output Ligand Coverage and Volume Coverage Score
    LC = how many ligand atoms close divided by total ligan atoms
    PC = homw many pocket vertex close divided by total number of vertices
    '''
    nPocket = len(pocket) # Container which is a PQR line list or a vertices 
    nLig = len(ligandCoord) #previously read and accordinlgy filtered
    if isPQR:
        #Assuming PQR DOES NOT contains the CHAIN field (if missing one can add a dummy A chain in pocket reading phase..)
        # xInd = 5
        # yInd = 6
        # zInd = 7
        thresholdD = thresholdForPQR

    else: 
        # xInd=0
        # yInd=1
        # zInd=2
        thresholdD = thresholdForTRIANG
    pocketCoords = np.array(pocket)
    d,_flag = Pdist_C(ligandCoord,pocketCoords)
    rowmap = np.ones(nPocket,bool)*False
    n_inLig=0
    for i in range(nLig):
        rw = d[i*nPocket:(i+1)*nPocket]<=thresholdD #= row: distance relations between atom ligand and all pocket atoms"
        rowmap = np.logical_or(rw,rowmap)
        n_inLig+=np.any(rw) #at least one per row = one hit
    n_inP= np.sum(rowmap)#can be seen as a mask where a true entry is a hit, a pocket atom within the distance threshold

    LC = n_inLig/nLig
    PC = n_inP/nPocket
    
    return LC,PC



def readP(filename,isPQR,testFolder):
    '''
    I'm assuming in the following a numbering scheme from smaller to larger for the different pocket flags
    '''
    triangName = testFolder+'/triangulatedSurf.off'
    rawList=[]
    seenFlag = set()
    inFile=open(filename,'r')
    if isPQR:
        xInd = 5
        yInd = 6
        zInd = 7
        flagIndex = 8
        for line in inFile:
            entries=line.split()
            # print(entries)
            try:
                flag = int(float(entries[flagIndex]))
            except:
                # print('skipping line:', line)
                continue
            if(flag<=0):
                #skip
                pass
            else:
                seenFlag.add(flag) 
                x = float(line.split()[xInd])
                y = float(line.split()[yInd])
                z = float(line.split()[zInd])

                rawList.append((flag,[x,y,z]))
    else:
        flagIndex = 0 
        xInd=0
        yInd=1
        zInd=2
        triangLines=loadTriang(triangName)
        for index,line in enumerate(inFile):
            entries=line.split()
            if(len(entries)>1):
                print("unexpected. CHECK")
                exit()
            line = triangLines[index]
            try:
                flag = int(float(entries[flagIndex]))
            except:
                # print('skipping line:', line)
                continue
            if(flag<=0):
                #skip
                pass
            else:
                seenFlag.add(flag) 
                x = float(line.split()[xInd])
                y = float(line.split()[yInd])
                z = float(line.split()[zInd])
                rawList.append((flag,[x,y,z]))


    print('tags:',seenFlag)

    #Not guarantees lines per pocket are ordered, that's why I proceed like this
    pocketList =[]
    for rank in sorted(seenFlag):
        cleanList = list(filter(lambda x: x[0]==rank,rawList))
        cleanList = [x[1] for x in cleanList]
        pocketList.append(cleanList)
    
    print('number of putative pockets: ', len(pocketList))

    return pocketList

import sys
def main():
    
    import glob
    args= sys.argv
    
    if(len(args)>1):
        candidateFolder = args[1] #str(sys.argv[1])
    else:
        print("Provide a candidate folder inline!")
        sys.exit()
    mapInfo = getStructureIterator('testMap.txt')
    print('Number of expected test structures:',len(mapInfo))
    # print(mapInfo)

    infileList=[n for n in glob.glob(candidateFolder+'/*')]
    # Establish type (PQR or TXT=triangulation) and assign number for map..
    enquire = infileList[0]
    if(re.match('^.*[.](pqr)$',enquire)):
        isPQR=True
        format = '.*[.](pqr)$'
        print('**pqr format**')
    elif(re.match('^.*[.](txt)$',enquire)):
        isPQR=False
        format = '.*[.](txt)$'
        print('** Triangulation vertices format**')
    else:
        print('Input format not recognized')
        sys.exit()
    input()
    folderNumber=[]
    for line in infileList:
        folderNumber.append([re.match('[^\d]*([\d]+)'+format,line).groups()[0],line])

    #MAIN LOOP
    #COUNTERS
    norm = 0 
    hitTop1 = 0
    hitTop3 = 0
    hitTop10 = 0
    counterLC = 0
    counterPC = 0
    nPockets = 0

    noHitMap=[]
    structureCounter = 0
    for fn in folderNumber:
        
        # #DEBUG
        # if(structureCounter==3):
        #     break
        
        structureCounter +=1

        number = int(fn[0])
        inFile =fn[1]
        #ASSES ORIGINAL STRUCTURE NAME AND CORRESPONDING LIGANDS
        pqrName = mapInfo[number]['pqr']
        print('\n STRUCTURE '+pqrName)
        ligands = mapInfo[number]['ligands']
        print("LIGANDS:"+str(ligands))
        print('File:',inFile)
        ligandCoords = []
        # Filter ligand (heavy atoms) close 5 AA from any protein atom
        for ln in ligands:
            ligandCoords.append(readLigand(ln,pqrName))
        nLigands = len(ligandCoords)
        pList = readP(inFile,isPQR,'testData/'+str(number)) #last argument useful only if triangulation has to be read
        nPockets += len(pList)
        r = 0 #RANKING POSITION
        norm += nLigands #UPDATE NORMALIZATION

        #RANKING LOOP WITHIN A STRUCTURE
        nHits = 0
        hitLig = np.zeros(nLigands) #keeps trace of which ligand has been hit
        for pn,pcoord in enumerate(pList):
            gotHIT=False
            for ind,ln in enumerate(ligands):
                # print("LIGAND:  ",ln)
                lcoord = ligandCoords[ind]
                LC,PC = getClose(lcoord,pcoord,isPQR) 
                if((np.round(LC,2)>=lCoverageTH)and(np.round(PC,2)>=pCoverageTH)):
                    print("pocket %d hit!"%(pn+1))
                    print("LC=%.3f\tPC=%.3f"%(LC,PC))
                    gotHIT = True
                    nHits+=1
                    if(r<10):
                        hitLig[ind] = 1
                        hitTop10+=1
                        counterLC += LC
                        counterPC += PC
                        if(r<3):
                            hitTop3+=1
                            if(r==0):
                                hitTop1 += 1
            if (r>9):
                #Above rank 10 limit
                break
            if (not gotHIT):
                r+=1 #COUNTS MISSING 
            if (nHits==nLigands):
                print("All ligands found")
                break
        if (nHits<nLigands):
            print("NOT all ligands found")
            for ln,hitMap in enumerate(hitLig):
                if(hitMap==0):
                    noHitMap.append([pqrName+': '+str(ligands[ln])])
        
        #COMPUTE (current) AVERAGES
        avHitTop1 = hitTop1/norm
        avHitTop3 = hitTop3/norm
        avHitTop10 = hitTop10/norm
        avPC = counterPC/hitTop10
        avLC = counterLC/hitTop10
        avNpockets = nPockets/structureCounter
        # print("Current averages: top1 = %.2f\ttop3 = %.2f\ttop10 = %.2f\taverageLC = %.2f\taveragePC = %.2f     "
        # %(np.round(avHitTop1,4)*100,np.round(avHitTop3,4)*100,np.round(avHitTop10,4)*100,np.round(avLC,4)*100,np.round(avPC,4)*100), end='\r')
    
    #OUT OF MAIN LOOP
    print("Analysis over")
    print("Averages: top1 = %.2f\ttop3 = %.2f\ttop10 = %.2f\taverageLC = %.2f\taveragePC = %.2f"
        %(np.round(avHitTop1,4)*100,np.round(avHitTop3,4)*100,np.round(avHitTop10,4)*100,np.round(avLC,4)*100,np.round(avPC,4)*100))
    print("NORM = ", norm)
    print("Strucures considered: ",structureCounter)
    #PRINT ON FILES
    outStat = open("rankStats.txt","w")

    outStat.write("#date: "+Crono().init())
    outStat.write("\n#top1\ttop3\ttop10\tLC\tPC\tNP\n")
    outStat.write("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f"
    %(np.round(avHitTop1,4)*100,np.round(avHitTop3,4)*100,np.round(avHitTop10,4)*100,np.round(avLC,4)*100,np.round(avPC,4)*100,np.round(avNpockets,1)))
    outStat.write('\n#NORM = '+str(norm))
    outStat.write('\n#Structures analyzed = '+str(structureCounter))
    outStat.write('\n#Thresholds: LC_Th=%.1f\t PC_Th=%.1f'%(np.round(lCoverageTH,3)*100,np.round(pCoverageTH,3)*100))

    outStat.close()

    failFile = open("failureList.txt","w")
    failFile.write("#Total failures: %d/%d\n"%(len(noHitMap),structureCounter))
    failFile.write('#Thresholds: LC_Th=%.1f\t PC_Th=%.1f\n'%(np.round(lCoverageTH,3)*100,np.round(pCoverageTH,3)*100))
    failFile.write("#List of missed ligands:\n")
    for string in noHitMap:
        failFile.write(repr(string).replace("[","\n").replace("]",""))
    
    failFile.close()

    print("**DONE**")

    return


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("\nUser exit")
        sys.exit()
