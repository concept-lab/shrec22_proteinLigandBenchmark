import numpy as np
import re
import os 
import subprocess
from  C_functs import Pdist_C,getIndex

THR = 5
OUTFOLDER = 'filtered_ligands/'

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

    structures = []
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
        structures.append({'pqr':name,'ligands': ligand_names})
        s+=n_ligands +1
            
    return structures

def filterLigand(lname,structureName):
    
    resMap = get_protein('allStructures/'+structureName)
    proteinCoord = np.empty((0,4))
    for i in resMap:
        c = np.append(np.asarray(i['coord']),i['radius'])
        proteinCoord = np.vstack([proteinCoord,c])
    inFile = open('ligands/'+lname+".xyz",'r')
    ligand_coord = np.loadtxt(inFile)

    d,flag = Pdist_C(ligand_coord[:,0:3],proteinCoord[:,0:3])
    index = np.where(d<=THR)[0]
    lindex,_pindex=getIndex(flag,index,proteinCoord.shape[0])
    # print(lindex)
    # print(np.unique(lindex))
    # print(ligand_coord[np.unique(lindex)])
    ligand_coord = ligand_coord[np.unique(lindex)]

    # Create filtered ligand

    np.savetxt(OUTFOLDER+lname+".xyz",ligand_coord, fmt ='%.5f', delimiter='\t' )
    
    return 


isFolder = os.path.isdir(OUTFOLDER)
if not isFolder:
    print("creating output folder.\nThis contains ligands filtered keeping atoms within "+str(THR)+"A from any protein atom.")
    subprocess.run(['mkdir',OUTFOLDER])
mapInfo = getStructureIterator('testMap.txt')
print(len(mapInfo))
counter = 0
ln=0
for m in mapInfo:
    print(m)
    structName = m['pqr']
    ligands = m['ligands']
    for lname in ligands:
        filterLigand(lname,structName)
        ln+=1
    counter+=1

print("DONE")
print("Number of structures = ",counter)
print("Number of ligands =",ln )