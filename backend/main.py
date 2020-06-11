import os

from flask import Flask, jsonify, request
from flask_cors import CORS
app = Flask(__name__)
CORS(app)


import pandas as pd
import numpy as np
'''
import pymysql
from DBUtils.PooledDB import PooledDB
dbPool = PooledDB(creator   = pymysql,
                host      = "localhost",
                user      = "root",
                password  = "root",
                database  = "edn",
                autocommit    = True,
                charset       = "utf8mb4",
                cursorclass   = pymysql.cursors.DictCursor,
                blocking      = False,
                maxconnections = 10)
'''

import sqlite3

#Used to convert DB entris from list to dictionary: column-value
def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def _connect():
    path_db = "./ic.db"
    conn = sqlite3.connect(path_db)
    conn.text_factory = str
    conn.row_factory = dict_factory
    return conn

def _execute(query, params=None):
    _conn = _connect()
    cursor = _conn.cursor()
    if (params):
        cursor.execute(query, params)
    else:
        cursor.execute(query)
    _conn.commit()
    return cursor

@app.route('/api/tissues/list', methods=["GET"])
def get_tissues():
    query="""SELECT
    Name AS name,
    BTOid AS BTOID,
    isCancerTissue AS isCancer 
    FROM Tissue
    ORDER BY name COLLATE NOCASE ASC;"""

    tissues=_execute(query).fetchall()
    return jsonify(tissues)

def _uniqueQuery(params):
    #Select genes expressed in the selected tissues
    query="""
    SELECT 
    e.TissueName AS tissue,
    e.ProteinUniProtAccNum,
    e.ExprLevelQual,
    e.ExprLevel,
    p.GeneSymbol AS geneSymbol,
    p.Name AS geneName,
    p.InBETSE AS inBETSE,
    s.SpecificityScore AS specificityScore
    FROM `expression` e
    LEFT JOIN `Protein` p
    ON (e.ProteinUniProtAccNum = p.UniProtAccNum)
    LEFT JOIN `Specificity` s
    ON (s.TissueName=e.TissueName AND s.UniProtAccNum=e.ProteinUniProtAccNum)
    WHERE e.ProteinUniProtAccNum IN (
        SELECT DISTINCT
        e.ProteinUniProtAccNum
        FROM `expression` e
        WHERE ((e.TissueName=? OR e.TissueName=?) AND (e.ExprLevel > 0 OR e.ExprLevelQual != 'Not detected'))
    );
    """
    resp=_execute(query, params).fetchall()
    resp=pd.DataFrame(resp)

    res=pd.DataFrame()
    #for each gene, check if uniquely expressed in the given tissue
    for geneUniprot in resp['ProteinUniProtAccNum'].unique():
        gene_data=resp[resp['ProteinUniProtAccNum']==geneUniprot]
        #drop all categorical expressions containing 'not expressed'
        gene_data.drop(gene_data[(gene_data['ExprLevelQual']=="Not detected")].index, inplace=True)

        #drop all quantitative expression equals 0
        gene_data.drop(gene_data[(gene_data['ExprLevel']==0) & (gene_data['ExprLevelQual']=="")].index, inplace=True)

        #if only records remaning are about the selected tissues, save it
        if (not gene_data.empty and (len(gene_data)==len(gene_data[gene_data['tissue'].isin(list(params))]))):
            res=res.append(gene_data)

    data=[]
    for _, row in res.iterrows():

        data.append({
            'tissue': row['tissue'],
            'exprType': 'quantitative' if row['ExprLevelQual'] == "" else 'qualitative',
            'exprLevel': row['ExprLevel'] if row['ExprLevelQual'] == "" else row['ExprLevelQual'],
            'geneSymbol': row['geneSymbol'],
            'geneName': row['geneName'],
            'specificityScore': row['specificityScore']
        })
    return data

def _comprehensiveQuery(params):
    query="""SELECT
    e.TissueName AS tissue,
    CASE e.ExprLevelQual WHEN '' THEN e.ExprLevel ELSE e.ExprLevelQual END AS exprLevel,
    CASE e.ExprLevelQual WHEN '' THEN 'quantitative' ELSE 'qualitative' END AS exprType,
    p.GeneSymbol AS geneSymbol,
    p.Name AS geneName,
    p.InBETSE AS inBETSE,
    s.SpecificityScore AS specificityScore
    FROM `expression` e
    LEFT JOIN `Protein` p
    ON (e.ProteinUniProtAccNum = p.UniProtAccNum)
    LEFT JOIN `Specificity` s
    ON (s.TissueName=e.TissueName AND s.UniProtAccNum=e.ProteinUniProtAccNum)
    WHERE (e.TissueName=? OR e.TissueName=?);
    """
    return _execute(query, params).fetchall()

@app.route('/api/tissues/expressions', methods=["POST"])
def get_ion_channel_expression_data():
    content = request.json

    params=(content['tissueA'], content['tissueB'])

    if (content['type']=="comprehensive"):
        expressions=_comprehensiveQuery(params)
    else:
        expressions=_uniqueQuery(params)

    if (len(expressions)==0):
        return jsonify({
            "genes": [],
            "genesInfo": {},
            "tissueA": {},
            "tissueB": {}
        })
    #Get unique genes from results
    data=pd.DataFrame(expressions)
    genes=list(data['geneSymbol'].unique())

    #Build genes info
    genesInfo={}
    for geneSymbol in genes:
        records=data[(data['geneSymbol']==geneSymbol)]
        genesInfo[geneSymbol]={
            "inBETSE": False if records.iloc[0]['inBETSE']=='N' else True,
            "geneName": records.iloc[0]['geneName']
        }

    #build TissuesInfo
    tissueData={
        content['tissueA']:{},
        content['tissueB']:{}
    }
    for geneSymbol in genes:
        geneRecords=data[(data['geneSymbol']==geneSymbol)]
        for tissueName in [content['tissueA'], content['tissueB']]:
            tissueRecords=geneRecords[geneRecords['tissue']==tissueName]

            if (tissueName=="--"):
                tissueData[tissueName][geneSymbol]={
                    'specificityScore': None,
                    'QtExpression_mean': None, 
                    'QtExpression_std': None, 
                    'catExpression': None}
            else:
                specificityScore=tissueRecords['specificityScore'].astype(float).mean()
                specificityScore=None if np.isnan(specificityScore) else specificityScore

                #Get quantitative expressions
                quantitativeExpressions=tissueRecords[tissueRecords['exprType']=='quantitative']
                if (len(quantitativeExpressions)==0):
                    QtExpression_mean=None
                    QtExpression_std=None
                else:
                    QtExpression_mean=quantitativeExpressions['exprLevel'].astype(float).mean()
                    QtExpression_std=quantitativeExpressions['exprLevel'].astype(float).std()
                    QtExpression_std=None if np.isnan(QtExpression_std) else QtExpression_std

                #Get qualitative expressions
                qualitativeExpressions=tissueRecords[tissueRecords['exprType']=='qualitative']
                if (len(qualitativeExpressions) > 0):
                    grouped=qualitativeExpressions.groupby('exprLevel').count().reset_index()

                    catExpression=[0, 0, 0, 0]
                    for _, row in grouped.iterrows():
                        exprLevel=row['exprLevel']
                        qt=row['tissue']
                        if (exprLevel=='Not detected'):
                            catExpression[0]=qt
                        elif (exprLevel=='Low'):
                            catExpression[1]=qt
                        elif (exprLevel=='Medium'):
                            catExpression[2]=qt
                        else:
                            catExpression[3]=qt
                else:
                    catExpression=None

                tissueData[tissueName][geneSymbol]={
                    'specificityScore': specificityScore,
                    'QtExpression_mean': QtExpression_mean, 
                    'QtExpression_std': QtExpression_std, 
                    'catExpression': catExpression}

    return jsonify({
        "genes": genes,
        "genesInfo": genesInfo,
        "tissueA": tissueData[content['tissueA']],
        "tissueB": tissueData[content['tissueB']]
    })

@app.route('/api/interactions', methods=["POST"])
def get_protein_compound_interactions():
    content = request.json
    geneSymbol=content['geneSymbol']

    query="""SELECT
    p.UniProtAccNum AS uniprotACC,
    i.ActionType AS actionType,
    i.AssayValue AS assayValue,
    i.AssayUnits AS assayUnits,
    i.AssayStandardType AS assayType,
    c.ChemblId AS chemblID,
    c.Name AS compoundName,
    c.CompoundType AS compoundType
    FROM `Protein` p
    LEFT JOIN `interaction` i
    ON (i.TargetUniProtAccNum=p.UniProtAccNum)
    LEFT JOIN `Compound` c
    ON (c.id=i.compoundID)
    WHERE p.GeneSymbol=?
    """
    params=(geneSymbol, )

    resp=_execute(query, params).fetchall()

    resp=pd.DataFrame(resp)

    #prepare data
    resp['assayValue_max']=resp['assayValue'] #create new columns to simplify group min and max
    resp.loc[resp['actionType']=="", 'actionType']="UNKNOWN" #set not available action type as unknown

    compoundIDs=list(resp['chemblID'].unique())

    data={}
    if (compoundIDs[0]==None): 
        compoundIDs=[]

    for comp_id in compoundIDs:
        compoundData=resp[resp['chemblID']==comp_id]

        #group assay by action type
        compoundActionTypes=compoundData.groupby('actionType').count().reset_index()[['actionType', 'assayValue']]

        compoundActionTypes=compoundActionTypes.rename({
            'actionType': 'actionType',
            'assayValue': 'count'
        }, axis='columns')
        compoundActionTypes.loc[:, 'tot']=compoundActionTypes['count'].sum()

        #group assay measurement values by type of assay measurement
        compoundValues=compoundData.groupby('assayType').agg({
            'assayValue': 'min',
            'assayValue_max': 'max',
            'assayUnits': 'first'
        }).reset_index()

        compoundValues= compoundValues.rename({
            'assayType': 'assayType',
            'assayValue': 'min',
            'assayValue_max': 'max',
            'assayUnits': 'units'
        }, axis='columns')

        data[comp_id]= {
                "name": compoundData.iloc[0]['compoundName'] if compoundData.iloc[0]['compoundName'] != "" else compoundData.iloc[0]['chemblID'],
                "chemblID": compoundData.iloc[0]['chemblID'],
                "compoundType": compoundData.iloc[0]['compoundType'],
                "assayValues": compoundValues.to_dict('records'),
                "compoundActionTypes": compoundActionTypes.to_dict('records'),
                "url_ref": "https://www.ebi.ac.uk/chembl/compound_report_card/"+comp_id
            }

    return jsonify({
        "success": True,
        "data": {
            "compoundIDs": compoundIDs,
            "compoundData": data
        }
    })

if __name__ == '__main__':
    debug=False if os.environ.get('PORT') else True
    app.run(debug=debug, port=int(os.environ.get('PORT', 5000)))
    #get_ion_channel_expression_data()