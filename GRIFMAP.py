#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 15:49:23 2024

@author: tanyapawliszak

Description: This is script is used to create a pipeline, where it uses the
terminal command for INFERNAL and use the results from it to check what are
the best possible introns from the given tblout.
"""
import os
import re
import sys
from Bio import SeqIO
import pandas as pd
from collections import defaultdict
from pathlib import Path
from tabulate import tabulate

###############################################################################

## Returns the results of INFERNAL and converts them to a dataframe list (pandas)
## TODO: DO NOT TOUCH
def get_infernal_df(fileResults):
    #we will make a pandas DF
    table_rows = []
    with open(fileResults) as f:
        col_names, positions = f.readline().strip(), f.readline().strip()
        positions = [0] + [i for i, x in enumerate(positions) if x == ' '] + [len(positions)]
        col_names = [x.strip() for x in cut_string_by_pos(col_names, positions)]

        for line in f:
            if line.startswith('#'):
                continue
            line = [x.strip() for x in cut_string_by_pos(line, positions)]
            table_rows.append(line)

    df = pd.DataFrame(table_rows, columns=col_names)

    df['seq from'], df['seq to'] = df['seq from'].apply(int), df['seq to'].apply(int)
    return df

###############################################################################

## Returns the results of HMMER and converts them to a dataframe list (pandas)
## TODO: DO NOT TOUCH
def get_hmmer_df(fileResults):
    #we will make a pandas DF
    table_rows = []
    with open(fileResults) as f:
        col_names, positions = f.readline().strip(), f.readline().strip()
        positions = [0] + [i for i, x in enumerate(positions) if x == ' '] + [len(positions)]
        col_names = [x.strip() for x in cut_string_by_pos(col_names, positions)]

        for line in f:
            if line.startswith('#'):
                continue
            line = [x.strip() for x in cut_string_by_pos(line, positions)]
            table_rows.append(line)

    df = pd.DataFrame(table_rows, columns=col_names)

    df['alifrom'], df['ali to'] = df['alifrom'].apply(int), df['ali to'].apply(int)
    df['hmmfrom'], df['hmm to'] = df['hmmfrom'].apply(int), df['hmm to'].apply(int)

    return df

###############################################################################

## Cuts string to different sections based on the given positions
## TODO: DO NOT TOUCH
def cut_string_by_pos(string, positions):
    out = [string[positions[i]:positions[i+1]] for i in range(0, len(positions) - 1)]
    return out

###############################################################################

## This returns a dataframe list of the best intron hits based on their E-value
def best_intron_output(dfResults):
    
    #storing the best hits in a different dataframe (df) list
    bestdf = pd.DataFrame(columns=dfResults.columns.values)
    
    #iterating through the df list and comparing all of the introns with each other
    for index, x in dfResults.iterrows():
        
        #if bestdf list is empty, add our first row
        if bestdf.empty:
           bestdf.loc[index] = x
        #otherwise, start comparing them
        else:
            overlap = False #checks if there is overlap
            insertRow = False #check if we can insert row
            removeIndex = 0 #keep track which index to remove
            
            #Go through the current df we have to check if there is any overlap
            for yndex, y in bestdf.iterrows():
                #Checking if the two rows we are comparing does overlap
                if x['seq from'] <= y['seq to'] and y['seq from'] <= x['seq to']:
                    overlap = True #there is an overlap!
                    #Now compare the E-values between the 2 rows
                    if float(x['E-value']) < float(y['E-value']):
                        insertRow = True #we can insert the row
                        removeIndex = yndex #which row we need to drop
                        
            #if no overlap found, just insert row
            if not overlap:
                bestdf.loc[index] = x
            #otherwise, add our new row and drop the one that had the lower E-value
            elif overlap and insertRow:
                bestdf.loc[index] = x
                bestdf = bestdf.drop(removeIndex, axis='index')

    return bestdf #return the best intron hits df list

###############################################################################

## This returns the dataframe with the best Twin-trons hits and their updated 
## intron start and end positions
def best_twintrons(dfResults):
        
    #copied the df results we had gotten and sort the new column 'start' in decending order
    besttwintronsdf = dfResults.sort_values('start', ascending=False)
    besttwintronsdf['E-value'] = besttwintronsdf['E-value'].apply(float)
    #print(besttwintronsdf)
    #print(dfResults)
    
    #this is where we will store our best twin-tron results for each possible flag segment
    besttwins = pd.DataFrame(columns=dfResults.columns)
    starttwin = 0
    endtwin = 0
    
    #looping through our sorted besttwintronsdf dataframe
    for xindex, x in besttwintronsdf.iterrows():
        
        #if our temp dataframe is empty and the current index is a possible twin-tron
        #(x['start'] != 0)
        if besttwins.empty and (x['Twin-trons'] == 'Yes'):
            besttwins.loc[xindex] = x #add it at this index
            starttwin = x['start'] #save the start twin-tron position from the first gene
            endtwin = x['end'] #save the end twin-tron position from the second gene
        else: #if not empty
            #check if the current row are under the same segment
            if (starttwin == x['start']) and (endtwin == x['end']):
                besttwins.loc[xindex] = x #add it to the list of best twin-trons results
            else: #if the current row is not part of the same segment
                #if there is only 1 introns in the possible twin-tron segment
                if len(besttwins.index) == 1:   
                    #extend the start and end position of the 1 intron only
                    dfResults.loc[besttwins['position'].iloc[0], 'seq from'] = (starttwin + 1)
                    dfResults.loc[besttwins['position'].iloc[0], 'seq to'] = (endtwin - 1)
                    #empty the list
                    besttwins = besttwins.iloc[0:0]
                    #if the current row is a possible intron
                    if (x['Twin-trons'] == 'Yes'):
                        besttwins.loc[xindex] = x #add it at this index
                        starttwin = x['start']  #save the start twin-tron position from the first gene
                        endtwin = x['end'] #save the end twin-tron position from the second gene
                
                #if there are 2 introns in the possible twin-tron segment
                elif len(besttwins.index) == 2:
                    #sort them in ascending order based on the seq from column
                    besttwins = besttwins.sort_values('seq from')
                    #retrieve the 2 values to calculate the middle of those 2 introns
                    intronendOne = dfResults.loc[besttwins['position'].iloc[0], 'seq to'] #get the 1st intron end position
                    intronstartTwo = dfResults.loc[besttwins['position'].iloc[1], 'seq from'] #get the 2nd intron start position
                    
                    #calculate the midpoint of these 2 introns
                    midPoint = int((intronendOne + intronstartTwo)/2)
                    
                    #extend the start and end position of the first intron (from gene1 end + 1 to the midpoint)
                    dfResults.loc[besttwins['position'].iloc[0], 'seq from'] = (starttwin + 1)
                    dfResults.loc[besttwins['position'].iloc[0], 'seq to'] = midPoint
                    
                    #extend the start and end position of the second intron (from midpoint + 1 to the gene2 start - 1)
                    dfResults.loc[besttwins['position'].iloc[1], 'seq from'] = (midPoint + 1)
                    dfResults.loc[besttwins['position'].iloc[1], 'seq to'] = (endtwin - 1)
                    
                    #empty the list
                    besttwins = besttwins.iloc[0:0]
                    #if the current row is a possible intron
                    if (x['Twin-trons'] == 'Yes'):
                        besttwins.loc[xindex] = x #add it at this index
                        starttwin = x['start'] #save the start twin-tron position from the first gene
                        endtwin = x['end'] #save the end twin-tron position from the second gene
                
                #if there are more than 2 introns in the possible twin-tron segment
                elif len(besttwins.index) > 2:  
                    #first organise them by the lowest E-value (desending order)
                    besttwins = besttwins.sort_values('E-value')
                    besttwins = besttwins.head(2) #keep the best 2 hits
                    besttwins = besttwins.sort_values('seq from') #sort them again by seq from in ascending order
                    
                    #retrieve the 2 values to calculate the middle of those 2 introns
                    intronendOne = dfResults.loc[besttwins['position'].iloc[0], 'seq to'] #get the 1st intron end position
                    intronstartTwo = dfResults.loc[besttwins['position'].iloc[1], 'seq from'] #get the 2nd intron start position
                    
                    #calculate the midpoint of these 2 introns
                    midPoint = int((intronendOne + intronstartTwo)/2)
                    
                    #extend the start and end position of the first intron (from gene1 end + 1 to the midpoint)
                    dfResults.loc[besttwins['position'].iloc[0], 'seq from'] = (starttwin + 1)
                    dfResults.loc[besttwins['position'].iloc[0], 'seq to'] = midPoint
                    
                    #extend the start and end position of the second intron (from midpoint + 1 to the gene2 start - 1)
                    dfResults.loc[besttwins['position'].iloc[1], 'seq from'] = (midPoint + 1)
                    dfResults.loc[besttwins['position'].iloc[1], 'seq to'] = (endtwin - 1)
                    
                    #empty the list
                    besttwins = besttwins.iloc[0:0]
                    #if the current row is a possible intron
                    if (x['Twin-trons'] == 'Yes'):
                        besttwins.loc[xindex] = x #add it at this index
                        starttwin = x['start'] #save the start twin-tron position from the first gene
                        endtwin = x['end'] #save the end twin-tron position from the second gene
                

    return dfResults #return the updated intron dataframe

###############################################################################

## This returns a dataframe list of the best gene hits based on their E-value
def best_gene_output(dfResults):
    
    #storing the best hits in a different dataframe (df) list
    bestdf = pd.DataFrame(columns=dfResults.columns.values)
    
    #iterating through the df list and comparing all of the genes with each other
    for index, x in dfResults.iterrows():
        
        #if bestdf list is empty, add our first row
        if bestdf.empty:
           bestdf.loc[index] = x
        #otherwise, start comparing them
        else:
            overlap = False #checks if there is overlap
            insertRow = False #check if we can insert row
            removeIndex = 0 #keep track which index to remove
            
            #Go through the current df we have to check if there is any overlap
            for yndex, y in bestdf.iterrows():
                #Checking if the two rows we are comparing does overlap
                if x['alifrom'] <= y['ali to'] and y['alifrom'] <= x['ali to']:
                    overlap = True #there is an overlap!
                    #Now comparing the E-values between the 2 rows
                    if float(x['E-value']) < float(y['E-value']):
                        insertRow = True #we can insert the row
                        removeIndex = yndex #which row we need to drop

            #if no overlap found, just insert row
            if not overlap:
                bestdf.loc[index] = x
            #otherwise, add our new row and drop the one that had the lower E-value
            elif overlap and insertRow:
                bestdf.loc[index] = x
                bestdf = bestdf.drop(removeIndex, axis='index')

    return bestdf #return the best gene hits df list

###############################################################################

## This checks if the current Gene Results have paired genes and to remove
## any outliers (such as ones that does not have a pair gene)
def pair_gene(hmmerOrgResults):
    
    #copy the current HMMER results into to a different dataframe to compare the match genes
    pairedGeneDF = hmmerOrgResults.copy()
    pairedGeneDF['hmmfrom'], pairedGeneDF['hmm to'] = pairedGeneDF['hmmfrom'].apply(int), pairedGeneDF['hmm to'].apply(int)
        
    #create a new dataframe to store our updated gene list
    newhmmerpairDF = pd.DataFrame(columns=hmmerOrgResults.columns.values)
    newhmmerpairDF['hmmfrom'], newhmmerpairDF['hmm to'] = newhmmerpairDF['hmmfrom'].apply(int), newhmmerpairDF['hmm to'].apply(int)
    
    offset = 50 #our offset of how much of a difference of nucleotides
    cox1offset = 15
    
    
    jndex = 0
    while jndex < len(hmmerOrgResults):
        pairfound = False
        
        #special case for rnl gene
        if (hmmerOrgResults.at[jndex, '# target name'] == 'rnl'):
            newhmmerpairDF.loc[jndex] = hmmerOrgResults.loc[jndex]
            pairfound = True
            jndex += 1
        
        for yndex, y in pairedGeneDF.iterrows():
            if pairfound:
                break
           
            #special case for nad4L and nad5
            if ((hmmerOrgResults.at[jndex, '# target name'] == 'nad4L') and (y['# target name'] == 'nad5') and (hmmerOrgResults.at[jndex, 'strand'] == y['strand']) and (jndex < yndex) and (yndex-jndex < 6)):
                newhmmerpairDF.loc[jndex] = hmmerOrgResults.loc[jndex]
                newhmmerpairDF.loc[yndex] = y
                jndex = yndex
                pairfound = True
            
            #special case for cob and cytb
            elif ((hmmerOrgResults.at[jndex, '# target name'] == 'cob' or hmmerOrgResults.at[jndex, '# target name'] == 'cytb') and (y['# target name'] == 'cob' or y['# target name'] == 'cytb') and (jndex < yndex)):
                if ((y['hmmfrom'] - offset) <= (hmmerOrgResults.at[jndex, 'hmm to']) <= (y['hmmfrom'] + offset)) and (jndex < yndex) and (yndex-jndex < 6) and (hmmerOrgResults.at[jndex, 'strand'] == y['strand']):
                    newhmmerpairDF.loc[jndex] = hmmerOrgResults.loc[jndex]
                    newhmmerpairDF.loc[jndex, '# target name'] = 'cob'
                    newhmmerpairDF.loc[yndex] = y
                    newhmmerpairDF.loc[jndex, '# target name'] = 'cob'
                    jndex = yndex
                    pairfound = True
            
            #special case for cox1
            elif ((hmmerOrgResults.at[jndex, '# target name'] == 'cox1' and y['# target name'] == 'cox1') and (hmmerOrgResults.at[jndex, 'strand'] == y['strand']) and (jndex < yndex) and (yndex-jndex < 6)):
                if (hmmerOrgResults.at[(yndex + 1), 'hmmfrom'] == y['hmmfrom'] or (y['hmmfrom'] - cox1offset) <= (hmmerOrgResults.at[(yndex + 1), 'hmmfrom']) <= (y['hmmfrom'] + cox1offset)):
                    if float(hmmerOrgResults.at[(yndex + 1), 'E-value']) < float(y['E-value']):
                        newhmmerpairDF.loc[jndex] = hmmerOrgResults.loc[jndex]
                        newhmmerpairDF.loc[(yndex  + 1)] = hmmerOrgResults.loc[(yndex  + 1)]
                        jndex = (yndex + 1)
                        pairfound = True
                    else:
                        newhmmerpairDF.loc[jndex] = hmmerOrgResults.loc[jndex]
                        newhmmerpairDF.loc[yndex] = y
                        jndex = yndex
                        pairfound = True
                else:
                    newhmmerpairDF.loc[jndex] = hmmerOrgResults.loc[jndex]
                    newhmmerpairDF.loc[yndex] = y
                    jndex = yndex
                    pairfound = True
                
            #normal case
            #Checking if x[hmm to] falls between (y[hmmfrom] - offset) and (y[hmmfrom] + offset) and have the same gene name
            elif ((y['hmmfrom'] - offset) <= (hmmerOrgResults.at[jndex, 'hmm to']) <= (y['hmmfrom'] + offset)) and (hmmerOrgResults.at[jndex, '# target name'] == y['# target name']) and (jndex < yndex) and (yndex-jndex < 6) and (hmmerOrgResults.at[jndex, 'strand'] == y['strand']):
                newhmmerpairDF.loc[jndex] = hmmerOrgResults.loc[jndex]
                newhmmerpairDF.loc[yndex] = y
                jndex = yndex
                pairfound = True
            
        if not pairfound:
            jndex += 1
                
    return newhmmerpairDF
    

###############################################################################

## We are combining our results we had received from Infernal and Hmmer into
## one list. We are currently taking the following info from both df set:
## Infernal: #target name, query name, seq from, seq to, E-value, score
## Hmmer: # target name, query name, alifrom, ali to, E-value, score
def combined_list(bestInfernalResults, bestHmmerResults):
    
    #creating our combined dataframe list with the column info we want to keep
    #combinedDF = pd.DataFrame(columns=['# target name', 'query name', 'from', 'to', 'E-value', 'score', 'Possible Twin-tron found', 'Start', 'End'])    
    combinedDF = pd.DataFrame(columns=['# target name', 'query name', 'from', 'to', 'E-value', 'score', 'Possible Twin-tron found'])    

    
    #looping through our dataframe list of bestInfernalResults and take the necessary data 
    for yindex, y in bestInfernalResults.iterrows():
        #combinedDF.loc[len(combinedDF.index)] = [y['#target name'], y['query name'], y['seq from'], y['seq to'], y['E-value'], y['score'], y['Twin-trons'], y['start'], y['end']]
        combinedDF.loc[len(combinedDF.index)] = [y['#target name'], y['query name'], y['seq from'], y['seq to'], y['E-value'], y['score'], y['Twin-trons']]

    #looping through our dataframe list of bestHmmerResults and take the necessary data 
    for xindex, x in bestHmmerResults.iterrows():
        #combinedDF.loc[len(combinedDF.index)] = [x['# target name'], x['query name'], x['alifrom'], x['ali to'], x['E-value'], x['score'], '-', '-', '-'] 
        combinedDF.loc[len(combinedDF.index)] = [x['# target name'], x['query name'], x['alifrom'], x['ali to'], x['E-value'], x['score'], '-'] 

        
    #sort the data based on the 'from' column
    combinedDFsorted = combinedDF.sort_values('from')
    
    combinedDFsorted.reset_index(drop=True, inplace=True)
    
    #print(combinedDFsorted.to_string()) #print out the results
    print(tabulate(combinedDFsorted.values.tolist(), list(combinedDFsorted.columns), tablefmt = 'simple')) #print out the results in table format

###############################################################################
## This would return the best introns hits with their updated lengths
## based on the gene info we had received from hmmer
## TODO: Note - Used 'ali from' and 'ali to' columns since those ones were
## matching closely to the gene bank files we currently had
## TODO: AN INTRON MUST FALL INBETWEEN TWO OF THE SAME GENE FAMILY PAIR!!!
def intron_estimated_sequence_length(infernalResults, hmmerResults):
    
    #if the infernal results are empty, we exit out
    if infernalResults.empty:
        print("No INFERNAL hits!")
        sys.exit()
    #if the hmmer results are empty, we exit out
    if hmmerResults.empty:
        print("No HMMER hits!")
        sys.exit()
    
    averageIntronLength = 2850
    
    #storing our updated infernal hits in a different dataframe list
    infernaldf = pd.DataFrame(columns=infernalResults.columns.values)
    twintrondf = pd.DataFrame(columns=['Twin-trons', 'start', 'end', 'position'])
    twintrondf['start'], twintrondf['end'] = twintrondf['start'].apply(int), twintrondf['end'].apply(int)

           
    #this is our offset to get some of the intron hits that we wanted but were overlapping with other genes
    offset = 15
    
    #iterating our best infernal results we had received
    for xindex, x in infernalResults.iterrows():
        #if both results are not empty, we can continue
        #infernaldf.loc[xindex] = x #adding our intron right away in the new list
        foundRange = False #check if we have found our updated range
        twintrondf.loc[xindex] = 'None'
            
        #creating our temp dataframe list to track 3 genes at a time
        hmmerdf = None
        hmmerdf = pd.DataFrame(columns=hmmerResults.columns.values)
            
        #iterating our best hmmer results
        for yindex, y in hmmerResults.iterrows():
            #checks if we found our updated ranges
            if foundRange:
                break
            #check if our temp hmmerdf is less than 3. If so, add it to the temp hmmerdf
            if len(hmmerdf) < 3:
                hmmerdf.loc[yindex] = y
            #if the temp hmmerdf is of length 3, we can proceed to check where the introns
            #fall between which 2 genes
            if len(hmmerdf) == 3:
                #creating df series to compare the 3 hmmer gene results all at once
                firstPos = None
                firstPos = pd.DataFrame(columns=hmmerResults.columns.values)
                secondPos = None
                secondPos = pd.DataFrame(columns=hmmerResults.columns.values)
                thirdPos = None
                thirdPos = pd.DataFrame(columns=hmmerResults.columns.values)
                    
                #here we can compare the 3 gene results we currentky have
                for windex, w in hmmerdf.iterrows():
                    #store the first gene results
                    if firstPos.empty:
                        firstPos.loc[windex] = w
                        firstPos['alifrom'], firstPos['ali to'] = firstPos['alifrom'].apply(int), firstPos['ali to'].apply(int)
                    #store the second gene results
                    elif secondPos.empty:
                        secondPos.loc[windex] = w
                        secondPos['alifrom'], secondPos['ali to'] = secondPos['alifrom'].apply(int), secondPos['ali to'].apply(int)
                    #store the third gene results
                    elif thirdPos.empty:
                        thirdPos.loc[windex] = w
                        thirdPos['alifrom'], thirdPos['ali to'] = thirdPos['alifrom'].apply(int), thirdPos['ali to'].apply(int)
                        
                    #seqTo = best['seq to'].iloc[0]
                    #first check if any of the df series we created are not empty
                    if not firstPos.empty and not secondPos.empty and not thirdPos.empty:
                        #between first and second gene and the gene name for both does match, add it to the list
                        #and updated their start and stop position
                        if ((firstPos['# target name'].iloc[0] == secondPos['# target name'].iloc[0]) and 
                            (((firstPos['ali to'].iloc[0] - offset) <= x['seq from'] <= (secondPos['alifrom'].iloc[0] + offset)) or ((firstPos['ali to'].iloc[0] - offset) <= x['seq to'] <= (secondPos['alifrom'].iloc[0] + offset)))):
                            
                            #check if the length of the 2 genes are bigger than the averageIntronLength
                            if ((secondPos['alifrom'].iloc[0] - firstPos['ali to'].iloc[0]) > averageIntronLength):
                                infernaldf.loc[xindex] = x
                                twintrondf.loc[xindex, 'Twin-trons'] = 'Yes'
                                twintrondf.loc[xindex, 'start'] = firstPos['ali to'].iloc[0]
                                twintrondf.loc[xindex, 'end'] = secondPos['alifrom'].iloc[0]
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            else: #if not, proceed to extend the intron length
                                infernaldf.loc[xindex] = x
                                infernaldf.loc[xindex, 'seq from'] = (firstPos['ali to'].iloc[0] + 1)
                                infernaldf.loc[xindex, 'seq to'] = (secondPos['alifrom'].iloc[0] - 1)
                                twintrondf.loc[xindex, 'Twin-trons'] = 'No'
                                twintrondf.loc[xindex, 'start'] = 0
                                twintrondf.loc[xindex, 'end'] = 0
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            break
                                
                            #between second and third gene and the gene name for both does match, add it to the list
                            #and updated their start and stop position
                        elif ((secondPos['# target name'].iloc[0] == thirdPos['# target name'].iloc[0]) and
                              (((secondPos['ali to'].iloc[0] - offset) <= x['seq from'] <= (thirdPos['alifrom'].iloc[0] + offset)) or ((secondPos['ali to'].iloc[0] - offset) <= x['seq to'] <= (thirdPos['alifrom'].iloc[0] + offset)))):
                            
                            #check if the length of the 2 genes are bigger than the averageIntronLength
                            if ((thirdPos['alifrom'].iloc[0] - secondPos['ali to'].iloc[0]) > averageIntronLength):
                                infernaldf.loc[xindex] = x
                                twintrondf.loc[xindex, 'Twin-trons'] = 'Yes'
                                twintrondf.loc[xindex, 'start'] = secondPos['ali to'].iloc[0]
                                twintrondf.loc[xindex, 'end'] = thirdPos['alifrom'].iloc[0]
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            else: #if not, proceed to extend the intron length
                                infernaldf.loc[xindex] = x
                                infernaldf.loc[xindex, 'seq from'] = (secondPos['ali to'].iloc[0] + 1)
                                infernaldf.loc[xindex, 'seq to'] = (thirdPos['alifrom'].iloc[0] - 1)
                                twintrondf.loc[xindex, 'Twin-trons'] = 'No'
                                twintrondf.loc[xindex, 'start'] = 0
                                twintrondf.loc[xindex, 'end'] = 0
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            break
                                
                        #between first and third gene and the gene name for both does match and if the second gene name does not match to either of them
                        #add it to the list and updated their start and stop position
                        elif ((firstPos['# target name'].iloc[0] == thirdPos['# target name'].iloc[0]) and
                              (((firstPos['ali to'].iloc[0] - offset) <= x['seq from'] <= (thirdPos['alifrom'].iloc[0] + offset)) or ((firstPos['ali to'].iloc[0] - offset) <= x['seq to'] <= (thirdPos['alifrom'].iloc[0] + offset))) and
                              ((firstPos['# target name'].iloc[0] != secondPos['# target name'].iloc[0]) and (thirdPos['# target name'].iloc[0] != secondPos['# target name'].iloc[0]))):
                            
                            #check if the length of the 2 genes are bigger than the averageIntronLength
                            if ((thirdPos['alifrom'].iloc[0] - firstPos['ali to'].iloc[0]) > averageIntronLength):
                                infernaldf.loc[xindex] = x
                                twintrondf.loc[xindex, 'Twin-trons'] = 'Yes'
                                twintrondf.loc[xindex, 'start'] = firstPos['ali to'].iloc[0]
                                twintrondf.loc[xindex, 'end'] = thirdPos['alifrom'].iloc[0]
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            else: #if not, proceed to extend the intron length
                                infernaldf.loc[xindex] = x
                                infernaldf.loc[xindex, 'seq from'] = (firstPos['ali to'].iloc[0] + 1)
                                infernaldf.loc[xindex, 'seq to'] = (thirdPos['alifrom'].iloc[0] - 1)
                                twintrondf.loc[xindex, 'Twin-trons'] = 'No'
                                twintrondf.loc[xindex, 'start'] = 0
                                twintrondf.loc[xindex, 'end'] = 0
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            break
                        
                        #Special case where we have nad4L and nad5 one after the other and a large gap between them
                        elif ((firstPos['# target name'].iloc[0] == 'nad4L') and (secondPos['# target name'].iloc[0] ==  'nad5') and 
                            (((firstPos['ali to'].iloc[0] - offset) <= x['seq from'] <= (secondPos['alifrom'].iloc[0] + offset)) or ((firstPos['ali to'].iloc[0] - offset) <= x['seq to'] <= (secondPos['alifrom'].iloc[0] + offset)))):
                            #check if the length of the 2 genes are bigger than the averageIntronLength
                            if ((secondPos['alifrom'].iloc[0] - firstPos['ali to'].iloc[0]) > averageIntronLength):
                                infernaldf.loc[xindex] = x
                                twintrondf.loc[xindex, 'Twin-trons'] = 'Yes'
                                twintrondf.loc[xindex, 'start'] = firstPos['ali to'].iloc[0]
                                twintrondf.loc[xindex, 'end'] = (secondPos['alifrom'].iloc[0] - 30)
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            else: #if not, proceed to extend the intron length
                                infernaldf.loc[xindex] = x
                                infernaldf.loc[xindex, 'seq from'] = (firstPos['ali to'].iloc[0] + 1)
                                infernaldf.loc[xindex, 'seq to'] = (secondPos['alifrom'].iloc[0] - 31)
                                twintrondf.loc[xindex, 'Twin-trons'] = 'No'
                                twintrondf.loc[xindex, 'start'] = 0
                                twintrondf.loc[xindex, 'end'] = 0
                                twintrondf.loc[xindex, 'position'] = xindex
                                foundRange = True
                            break

                        #if it does not fall inbetween the current genes, drop the first gene
                        else:
                            hmmerdf = hmmerdf.iloc[1:]
                            #hmmerdf.loc[yindex] = y
                                    
    infernaldf['Twin-trons']= twintrondf['Twin-trons']
    infernaldf['start'] = twintrondf['start']
    infernaldf['end'] = twintrondf['end']
    infernaldf['position'] = twintrondf['position']
    #print(infernaldf)
    updatedInfernaldf = best_twintrons(infernaldf)
    
    bestinfernaldf = best_intron_output(updatedInfernaldf) #filter our INFERNAL results again
    #print(bestinfernaldf)
    bestinfernaldf.reset_index(drop=True, inplace=True) #reset our indexes
    return bestinfernaldf

###############################################################################

## First it gets the INFERNAL results and converts it to a dataframe list.
## Then sorts the introns based on their 'seq from' column and takes the best
## introns hits based on their E-value and returns its as a dataframe list
def best_infernal_results(outputFile):
    df = get_infernal_df(outputFile) #First we get our INFERNAL data in dataframe format

    updatedDF = df.sort_values('seq from') #Sort them by 'seq from' in ascending order

    bestDF = best_intron_output(updatedDF) #Get the best hits
    
    bestDF.reset_index(drop=True, inplace=True) #Reset the dataframe indexes
    
    print(bestDF.to_string(index=False))

    return bestDF #return our filtered INFERNAL best hit list

###############################################################################

## First it gets the HMMER results and converts it to a dataframe list.
## Then sorts the genes based on their 'alifrom' column and takes the best
## gene hits based on their E-value, then removes any outlier genes w/out any 
## pairs and returns its as a dataframe list
def sort_hmmer_results(outputFile):
    df = get_hmmer_df(outputFile) #First we get our HMMER data in dataframe format

    updatedDF = df.sort_values('alifrom') #Sort them by 'alifrom' in ascending order
    
    #print(updatedDF.to_string())
    
    bestDF = best_gene_output(updatedDF) #Get the best hits
    
    bestDF.reset_index(drop=True, inplace=True) #Reset the dataframe indexes
    
    #print(bestDF.to_string())
    
    betterDF = pair_gene(bestDF) #Get only the pair genes
        
    betterDF.reset_index(drop=True, inplace=True) #Reset the dataframe indexes
    
    print(betterDF.to_string(index=False))

    return betterDF #return our filtered HMMER best hit list
    
###############################################################################

## Print our updated INFERNAL results to a tableout in a table format
def convert_to_tableout(infernalResults):
    
    outputTestFile = open(infernalResults['query name'].iloc[0] + "_GRIFMAP.tblout", 'w') #Create a tblout file with the current gene name
    
    infernalResults = infernalResults.iloc[:,:-4]
    #infernalResults.pop(infernalResults.columns[-1])
    outputTestFile.write(tabulate(infernalResults.values.tolist(), list(infernalResults.columns), tablefmt = 'simple')) #Print our table into the file
    
    outputTestFile.close() #close the file

###############################################################################         


if len(sys.argv) != 2:
    print("USAGE: python GRIFMAP.py fastaFile")
    sys.exit()

print("Welcome to GRIFMAP!\n")

print("In GRIFMAP, we are using a combination of both HMMER & Infernal to identify Group I Introns with a few adjustments to get better and accurate results!\n")

print("First we will be running Infernal with our improved models to identify all possible Group I introns from the file " + sys.argv[1])
print("Please note this process may take a while. Have a cup of coffee or tea while you wait! :) \n")

os.system("cmscan --tblout GII_Results.tblout AugmentedModels.cm " + sys.argv[1])

print("\nInfernal has finished scanning for all possible Group I Introns\n")
print("Running next stage of GRIFMAP")
print("Extracting the scan results from Infernal and finding the best identified Group I Introns for certain regions\n")

bestIntronOutput = best_infernal_results("GII_Results.tblout")

print("\nCollected best possible Group I intron hits from Infernal\n")

print("Running next stage of GRIFMAP")

print("We will be running HMMER with our models to identify all possible exons from the file " + sys.argv[1])
print("Please note this process may take a while. Have a cup of coffee or tea while you wait! :) \n")

os.system("nhmmscan -E 15 --tblout GIH_Results.tblout geneDatabasewoCountModel3 " + sys.argv[1])

print("\nHMMER has finished scanning for all possible exons\n")
print("Running next stage of GRIFMAP")
print("\nExtracting the HMMER results and finding the best identified exons for certain regions\n")

sortedHmmerOutput = sort_hmmer_results("GIH_Results.tblout")

print("\nSorting both Infernal and HMMER results\n")

updatedInferalResults = intron_estimated_sequence_length(bestIntronOutput, sortedHmmerOutput)

print("Printing GRIFMAP results:\n")

combined_list(updatedInferalResults, sortedHmmerOutput)

convert_to_tableout(updatedInferalResults)
