import sys
sys.path.append("../")
import utils_reconcILS.readWrite as readWrite
from ete3 import PhyloTree, Tree

file = open("/N/u/lkonig/Quartz/Desktop/dupcoal-main/example_4_18_2025/trees.tre")
red = readWrite.readWrite()

def copy_function(refTree, nTree): 
    refTreeObj = red.parse(refTree)
    refTreeObj.label_internal()
    nTree.label_internal()
    stack = [(refTreeObj, nTree)]
    while stack: 
        current_ref, current_nTree = stack.pop()
        current_ref.interval = current_nTree.interval
        if not current_ref.isLeaf and not current_nTree.isLeaf:
            if current_ref.leftChild.numbered_taxa == current_nTree.leftChild.numbered_taxa: 
                stack.append((current_ref.leftChild, current_nTree.leftChild)) 
                stack.append((current_ref.rightChild,current_nTree.rightChild))
            else: 
                stack.append((current_ref.leftChild, current_nTree.rightChild))
                stack.append((current_ref.rightChild, current_nTree.leftChild))
    return refTreeObj

def Traverse_key(myTree):
    if myTree == None: 
        return 
    else:
        if myTree.isLeaf: 
            myTree.key = str(1)
            return str(1) 
        else: 
            if myTree.leftChild != None and myTree.rightChild != None: 
                myTree.key = Traverse_key(myTree.leftChild) + Traverse_key(myTree.rightChild) + str(myTree.interval + 1)
                return myTree.key 

def Traverse(root, intervalSize): 
    stack = [root]
    while stack: 
        current_node = stack.pop()
        if current_node: 
            current_node.interval = int(current_node.branchlength // intervalSize)
            stack.append(current_node.leftChild)
            stack.append(current_node.rightChild)

def isequal(tree1, tree2):
        if tree1[-1] != ";":
            tree1 = tree1 + ";"
        if tree2[-1] != ";": 
            tree2 = tree2 + ";"
        firstTree = Tree(tree1)
        secondTree = Tree(tree2)
        node_set_firstTree = {node.name for node in firstTree.traverse() if node.name}
        node_set_secondTree = {node.name for node in secondTree.traverse() if node.name}
        if len(tree1) <= 4 or len(tree2) <= 4: 
            if tree1 == tree2: 
                return True
            else:
                return False 
        else: 
            rf, _, _, _, _, _, _ = firstTree.robinson_foulds(secondTree)
            
            if rf == 0 and node_set_firstTree==node_set_secondTree: 
                return True 
            else: 
                return False 

def isPresent(topology, topologyLst): 
    for value in topologyLst:
        if isequal(value, topology): 
            return value
    return None 

result = {}
visited = []

for line in file: 
    print(line)
    if line == "None\n": 
        continue
    line = line.strip().replace('e-', '0')
    varname = red.parse_bio_with_branch_length(line)
    Traverse(varname, 0.2)
    currentTree = red.to_newick2(varname) 
    if visited == []:
        visited.append(currentTree) 
    else: 
        flag=False
        for tr in visited: 
            if isequal(tr, currentTree): 
                if tr != currentTree: 
                    varname = copy_function(tr, varname)
                    flag=True
                    
        if not flag: 
            visited.append(currentTree)
    Traverse_key(varname)
    thisTree = PhyloTree(red.to_newick1(varname))
    # print(red.to_newick1(varname))
    nodeLst = [node.write(format = 9) for node in thisTree.traverse("postorder") if not node.is_leaf()]
    nodeLst1 = [node.dist for node in thisTree.traverse("postorder") if not node.is_leaf()]
    for i in range(len(nodeLst)): 
        if nodeLst1[i] not in result.keys(): 
            result[nodeLst1[i]] = {nodeLst[i]: 1}
        else: 
            temp = isPresent(nodeLst[i], result[nodeLst1[i]].keys()) # This is the topology
            if temp: 
                result[nodeLst1[i]][temp] = result[nodeLst1[i]][temp]+1
            else:
                result[nodeLst1[i]][nodeLst[i]] = 1

with open("newDictionary_4_26.txt", "w") as file1: 
    file1.write(str(result))

# print(result)

file.close()

