from collections import Counter
from itertools import permutations
from ete3 import Tree
import math 
import ast
import pandas as pd
import re
import matplotlib.pyplot as plt
import statistics


def find_duplication(i):
   
    with open("rep_"+str(i)+".log", "r") as file:
        content = file.read()

    
    match = re.search(r"Total duplications:\s*(\d+)", content)

    if match:
        total_duplications = int(match.group(1))
        return total_duplications
    else:
        return 0

# Brings trees to simplest form 
def hide_integers(inlst):
    newlst = inlst[0].strip()
    output = ""
    i = 0 
    while i < len(newlst): 
        if (newlst[i] != '_') and not newlst[i].isdigit() and (newlst[i] != '.') and (newlst[i] != ':'):
            output = output + newlst[i]
        i = i + 1
    return output 

def assignValues(inputString):
    temp = {}
    output = ""
    for char in inputString: 
        if char.isalpha() and char not in temp.keys(): 
            temp[char] = 1
    for char in inputString: 
        if not char.isalpha(): # if it's not a letter
            output = output + char 
        else: # if it is a letter
            output = output + char # add letter 
            output = output + "_" + str(temp[char])
            temp[char] = temp[char] + 1
    return output

def replacenodes(inTree, newordering): 
    currindex = 0 
    i = 0
    letters = [item.split('_')[0] for item in newordering]
    output = ""
    while i < len(inTree): 
        if inTree[i].isalpha(): # if current index is a letter
            if inTree[i] in letters: 
                output = output + newordering[currindex]
                currindex = currindex + 1 # go to next item in newordering 
                i = i + 3
            else: 
                output = output + inTree[i]
                i = i + 1
        else: 
            output = output + inTree[i]
            i = i + 1
    return output 

# Removes duplicate trees from the generated alternatives
def removeDuplicateTrees(treelst): 
    if len(treelst) == 1: 
        return treelst
    output = []
    visited = []
    for tree in treelst: 
        if tree not in visited: 
            sublst = [tree]
            tree1 = Tree(tree)
            visited.append(tree)  # Don't want to visit this tree anymore
            for subtree in treelst: 
                tree2 = Tree(subtree)
                rf, _, _, _, _, _, _ = tree1.robinson_foulds(tree2)
                if rf == 0 and subtree not in visited:  # Ensure we're not adding already visited trees
                    sublst.append(subtree)
                    visited.append(subtree)
            output.append(sublst)
    return [elem[0] for elem in output]

def generate_alt(sampleTree): #output = list containing alternative trees
    letter_list = [char for char in sampleTree if char.isalpha()]
    letter_counts = Counter(letter_list)
    duplicatedletters = []
    result = []
    if len(letter_list) <= 2 or max(list(letter_counts.values())) < 2: 
        result.append(sampleTree)
        return result 
    for index, char in enumerate(sampleTree): 
        if char.isalpha(): 
            if letter_counts[char] > 1: #Letter appears more than once
                duplicatedletters.append(char + sampleTree[index+1] + sampleTree[index+2])
    possible_ordering = permutations(duplicatedletters)
    for ordering in possible_ordering: 
        result.append(replacenodes(sampleTree, ordering))
    return result

# Removes trees which are completely new, not similar to previous
def find_valid(alt_tree_lst, fileTree): 
    outlst = []
    for tree in alt_tree_lst: 
        new_tree = hide_integers([tree])
        original_tree = hide_integers(fileTree)
        if original_tree == new_tree: 
            outlst.append(tree)
    return outlst 

def isequal(tree1, tree2):

    '''
    Changes made: added nodeset matching.

    This is because rf of  (B1,C1) and (B2,C1)==0 (because there arent any internal branches)
    Hence we take the set of all the node names and also match them.

    I have also changed the size for single taxa trees because len(A_1;)= 4 not 3. 
    '''
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
        
def isequal_value(tree1, tree2):
    '''
    This is similar to the function you wrote(isequal) but instead of True and False this outputs the rf distance.

    --< for single taxa tree for eg.(A_1;), if they are different then the distance is 1.

    
    '''
    if tree2[-1] != ";": 
        tree2 = tree2 + ";"

    firstTree = Tree(tree1)
    secondTree = Tree(tree2)
    
    if len(tree1) <= 4 or len(tree2) <= 4: 
        if tree1 == tree2: 
            return 0
        else:
            return 1 
    else: 
        #print('-------------------------->',firstTree.robinson_foulds(secondTree))
        rf, _, _, _, _, _, _ = firstTree.robinson_foulds(secondTree)
        #print('rf--------------->',firstTree,secondTree,rf)
        return rf

savedValues = {'True':[],'Predicted':[],'RF':[],'Total Duplication':[],'Mapped Interval':[],'Predictive Distribution':[]}
file2 = ast.literal_eval(open('myDictionary.txt', 'r').read()) # Converts string dictionary to dictionary
for i in range(1000):
    print(i)
    try:
        file = open("rep_" + str(i) + ".tre" , 'r').readlines()

        '''
        Read total duplication from the logs and only run the experiments
        if total duplication==1
        '''
        total_duplication=find_duplication(i)

        if total_duplication==0 or total_duplication>1:
            continue


        '''
        If the size of the gene tree is greater than 300(arbitary) then dont run it. It will be slow
        '''

        if len(file[0])>300:
            continue


        currentTree = Tree(file[0]) #file[0] is the true tree
        treeHeight = currentTree.get_farthest_node()[1]
        interval_quantity = math.ceil(treeHeight/0.2) #Finds top interval 
        
        #Generates alternative trees 
        cleaned_tree = removeDuplicateTrees(find_valid(generate_alt(assignValues(hide_integers(file))), file))
        
        resultDictionary = {}
        visited=[]

        '''
        Check if we find a match (xtree==ytree)

        If we dont find a match within an interval than we make the values zero 

        ------------------------------------------------------------

        It seems like we have duplicates in our training dictionary. To mitigate it:

        We check if xtree is in the resultDictionary. If we already have it we add the values
        if we dont have it we equate the values.

        Example where it works:

        For instance xtree= (A1,(B1,C1))
        file2[interval_quantity] = {(A1,(B1,C1)):2,(A1,(C1,B1)):3}

        Both (A1,(B1,C1)) and (A1,(C1,B1)) are equal so we add the values and make it 5.


        '''
        for xtree in cleaned_tree:
            found=False
            for ytree in file2[interval_quantity]: #find if alternative tree is in result dictionary
                if isequal(xtree, ytree): 
                    if xtree not in resultDictionary:
                        resultDictionary[xtree] = file2[interval_quantity][ytree]
                        
                    else:
                        resultDictionary[xtree]+= file2[interval_quantity][ytree]

                    found=True
                    visited.append(ytree)
                        

            if not found:
                resultDictionary[xtree] = 0


        '''
        This is the part where we are trying to get rid of the process where we propose a new alternative tree.
        Rather we are trying to match the tree tree(after removing intergers) with the keys in the interval dictionary(file2[interval_quantity]) 
        And give all the keys that matches.

        Problem is rf distance does not work with duplicates. We can use TreeKO.

        I found this workaround (may not work for all) where if you give ete3 two trees and print them in a format it print is same way. So we can simply equate it.
        
        '''

        for ytree in file2[interval_quantity]:
            tr=Tree(hide_integers(file))
            tr1=Tree(hide_integers([ytree])+';')
            if tr.write(format=1) == tr1.write(format=1) and ytree not in visited:
                resultDictionary[ytree] = file2[interval_quantity][ytree]
                found=True



            
        max_key=max(resultDictionary, key = resultDictionary.get)

        '''
        If the maxvalue of  resultDictionary  is 0 then we know all our values are also 0. 

        We dont want these to influence our results so have removed them.
        '''
        if resultDictionary[max_key]==0:
            continue

        '''
        Here we are normalizing the resultDictionary . 
        '''
        totalSum = sum(resultDictionary.values())
        resultDictionary = {key:(value/totalSum) for key,value in resultDictionary.items()}

        '''
        I have added bunch of extra outputs. Namely , 
        
        Mapped Interval (which is interval quantity)
        Total Duplication
        Predictive distribution==> Normalized resultDictionary.s

        '''

        #print(file2[interval_quantity])
        print("True Tree: ", file[0]) 
        print("Predicted Tree: ",max_key)
        print("------------------------------------------------------")
        savedValues['True']+=[file[0]]
        savedValues['Total Duplication']+=[total_duplication]
        savedValues['Predicted']+=[max_key]
        savedValues['RF']+=[isequal_value(file[0],max_key)]
        savedValues['Predictive Distribution']+=[resultDictionary]
        savedValues['Mapped Interval']+=[interval_quantity]
        print(resultDictionary)
        pd.DataFrame(savedValues).to_csv('./results.csv', index=True)
    except Exception as e: 
        print(e)

# Plot generation
categories = []
for value in savedValues['RF']: 
    if value not in categories: 
        categories.append(value)
counts = []
for categories_var in categories: 
    counts.append(savedValues['RF'].count(categories_var))
plt.bar(categories, counts, color = 'black')
plt.xlabel('RF Distance')
plt.ylabel('Values')
plt.title('Quantity of each RF')

plt.show()
print("Accuracy: " + str(savedValues['RF'].count(0)/len(savedValues['RF'])))
print("Total trees: " + str(len(savedValues['RF'])))
print("Mean: " + str(statistics.mean(savedValues['RF'])))
print("Standard deviation: " + str(statistics.stdev(savedValues['RF'])))
print("Quantity of 1 RF: " + str(savedValues['RF'].count(1)))
print("Quantity of 2 RF: " + str(savedValues['RF'].count(2)))



