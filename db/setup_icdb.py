import csv

import go_obo_parser

import icdb

import pandas as pd
import numpy as np

# Path for DB
PATH_DB = 'ic.db'

# Constants to be stored in the DB
ACTION_TYPE_ACTIVATOR = 'ACTIVATOR'
ACTION_TYPE_AGONIST = 'AGONIST'
ACTION_TYPE_ALLOSTERIC_ANTAGONIST = 'ALLOSTERIC ANTAGONIST'
ACTION_TYPE_ANTAGONIST = 'ANTAGONIST'
ACTION_TYPE_BLOCKER = 'BLOCKER'
ACTION_TYPE_INHIBITOR = 'INHIBITOR'
ACTION_TYPE_MODULATOR = 'MODULATOR'
ACTION_TYPE_NEGATIVE_ALLOSTERIC_MODULATOR = 'NEGATIVE ALLOSTERIC MODULATOR'
ACTION_TYPE_OPENER = 'OPENER'
ACTION_TYPE_PARTIAL_AGONIST = 'PARTIAL AGONIST'
ACTION_TYPE_POSITIVE_ALLOSTERIC_MODULATOR = 'POSITIVE ALLOSTERIC MODULATOR'
ACTION_TYPE_POSITIVE_MODULATOR = 'POSITIVE MODULATOR'
ACTION_TYPE_STABILISER = 'STABILISER'
EXPR_ASSAY_IMMUNO = 'immuno'
EXPR_ASSAY_MICROARRAY = 'microarray'
EXPR_ASSAY_RNASEQ = 'RNA-seq'
ASSAY_TYPE_ADMET = 'A'
ASSAY_TYPE_BINDING = 'B'
ASSAY_TYPE_FUNCTIONAL = 'F'
ASSAY_TYPE_TOXICITY = 'T'
ASSAY_TYPE_UNASSIGNED = 'U'
COMPOUND_TYPE_COMPOUND = 'compound'
COMPOUND_TYPE_DRUG = 'drug'
DATASET_BIOGPS = 'U133AGNF1B'
DATASET_HPA = 'hpa'
EXPR_H = 'High'
EXPR_L = 'Low'
EXPR_M = 'Medium'
EXPR_ND = 'Not detected'
EXTDB_BIOGPS = 'biogps'
EXTDB_BRENDA = 'brenda'
EXTDB_CHANNELPEDIA = 'channelpedia'
EXTDB_CHEMBL = 'chembl'
EXTDB_DRUGBANK = 'drugbank'
EXTDB_HMDB = 'hmdb'
EXTDB_HPA = 'hpa'
EXTDB_ZINC15 = 'zinc15'
NO = 'N'
YES = 'Y'
ACTION_TYPE_CHOICES = [
    ACTION_TYPE_ACTIVATOR,
    ACTION_TYPE_AGONIST,
    ACTION_TYPE_ALLOSTERIC_ANTAGONIST,
    ACTION_TYPE_ANTAGONIST,
    ACTION_TYPE_BLOCKER,
    ACTION_TYPE_INHIBITOR,
    ACTION_TYPE_MODULATOR,
    ACTION_TYPE_NEGATIVE_ALLOSTERIC_MODULATOR,
    ACTION_TYPE_OPENER,
    ACTION_TYPE_PARTIAL_AGONIST,
    ACTION_TYPE_POSITIVE_ALLOSTERIC_MODULATOR,
    ACTION_TYPE_POSITIVE_MODULATOR,
    ACTION_TYPE_STABILISER,
    ]
ASSAY_TYPE_CHOICES = [
    ASSAY_TYPE_ADMET,
    ASSAY_TYPE_BINDING,
    ASSAY_TYPE_FUNCTIONAL,
    ASSAY_TYPE_TOXICITY,
    ASSAY_TYPE_UNASSIGNED,
    ]

# Paths for data files
PATH_CANCER_TISSUES = 'data/bto/cancers.txt'
PATH_BIOGPS_GCRMA = 'data/biogps/U133AGNF1B.gcrma.avg.csv'
PATH_BIOGPS_GNF1H_ANNOT = 'data/biogps/gnf1h.annot2007.tsv'
PATH_BIOGPS_TRANSLATION_TABLE = 'data/biogps/biogps_translation_table.csv'
PATH_BTO_BRENDA_TISSUE_OBO = 'data/bto/BrendaTissueOBO.obo'
PATH_BTO_TISSUES = 'data/bto/tissues.txt'
PATH_CHANNEL_CLASSES = 'data/channel-classes/channel-classes.csv'
PATH_CHANNELPEDIA_INFO = 'data/channel-classes/channelpedia_info.csv'
PATH_CHEMBL_ASSAYS_COMPOUND = 'data/chembl-assays/output_compound-assays_v1.1.dat' # Not used
PATH_CHEMBL_ASSAYS_DRUG = 'data/chembl-assays/output_drug-assays_v1.1.dat' # Not used
PATH_CHEMBL_HUMAN_ASSAYS_COMPOUND = 'data/chembl-assays/output_human_compounds-assays_v1.1.dat'
PATH_CHEMBL_HUMAN_ASSAYS_DRUG = 'data/chembl-assays/output_human_drug-assays_v1.1.dat'
PATH_GO_ION_CHANNELS = 'data/channel-classes/master_ion_channel_list.csv'
PATH_GO_QUICKGO = 'data/go/QuickGO-ion-channel-COMBINED-human.dat'
PATH_GO_AMIGO = 'data/go/go_0006811_taxon_9606.dat'
PATH_HPA_NORMAL = 'data/hpa/normal_tissue.tsv'
PATH_HPA_PATHOLOGY = 'data/hpa/pathology.tsv'
PATH_HPA_RNA_TISSUE = 'data/hpa/rna_tissue.tsv'
PATH_HPA_TRANSLATION_TABLE = 'data/hpa/hpa_tissue_translation_table.csv'
PATH_IN_BETSE = 'data/channel-classes/in_betse.csv'
PATH_TARGET_COMPOUND = 'data/target-compound/target-compound.csv' # Not used
PATH_TC_CHEMBL_COMPOUND = 'data/target-compound/output-query_ChEMBL-uniprot-compound.dat' # Not used
PATH_TC_CHEMBL_DRUG = 'data/target-compound/output-query_ChEMBL-uniprot-drug.dat' # Not used
PATH_TC_DRUG_INFO_DRUGBANK = 'data/target-compound/output-organized_drug-info-Drugbank.dat' # Not used
PATH_TC_DRUG_INFO_HMDB = 'data/target-compound/output-organized_drug-info-HMDB.dat' # Not used
PATH_TC_TARGET_INFO_DRUGBANK = 'data/target-compound/output-organized_target-info-Drugbank.dat' # Not used
PATH_TC_TARGET_INFO_HMDB = 'data/target-compound/output-organized_target-info-HMDB.dat' # Not used
PATH_TC_TTD = 'data/target-compound/output-organized-TTD.dat' # Not used
PATH_UNIPROT_ENSEMBL_IDMAPPING = 'data/uniprot/HUMAN_9606_idmapping_Ensembl.dat' # Not used
PATH_UNIPROT_GENBANK = 'data/uniprot/GenBankUniProt.txt'
PATH_UNIPROT_HUMAN_IDMAPPING = 'data/uniprot/HUMAN_9606_idmapping.dat'


# DB setup routines
def _setup_external_db_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_external_db_table()
    db.add_external_db(EXTDB_BIOGPS, 'http://biogps.org')
    db.add_external_db(EXTDB_CHEMBL, 'https://www.ebi.ac.uk/chembl/')
    db.add_external_db(EXTDB_BRENDA, 'http://www.brenda-enzymes.org')
    db.add_external_db(EXTDB_DRUGBANK, 'https://www.drugbank.ca')
    db.add_external_db(EXTDB_ZINC15, 'http://zinc15.docking.org')
    db.add_external_db(EXTDB_HMDB, 'http://www.hmdb.ca')
    db.add_external_db(EXTDB_CHANNELPEDIA, 'http://channelpedia.epfl.ch')

def _setup_tissue_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_tissue_table()
    tissue_name_list = []
    cancer_name_list = []
    with open(PATH_BTO_TISSUES, 'rU') as tissue_file:
        for line in tissue_file:
            tissue_name_list.append(line.strip())
    with open(PATH_CANCER_TISSUES, 'rU') as cancer_file:
        for line in cancer_file:
            cancer_name_list.append(line.strip())
    bto_parser = go_obo_parser.parseGOOBO(PATH_BTO_BRENDA_TISSUE_OBO)
    for bto_record in bto_parser:
        tissue_name = bto_record['name']
        if tissue_name in tissue_name_list:
            bto_id = bto_record['id']
            if tissue_name in cancer_name_list:
                is_cancer = 'Y'
            else:
                is_cancer = 'N'
            db.add_tissue(tissue_name, bto_id, '', is_cancer)

def input_biogps_tissue_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_BIOGPS_TRANSLATION_TABLE, 'rU') as translation_file:
        translation_reader = csv.reader(translation_file)
        next(translation_reader)
        for row in translation_reader:
            biogps_tissue = row[0]
            bto_tissue = row[1]
            db.add_db_tissue(EXTDB_BIOGPS, bto_tissue, biogps_tissue)


def input_hpa_tissue_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_HPA_TRANSLATION_TABLE, 'rU') as translation_file:
        translation_reader = csv.reader(translation_file)
        next(translation_reader)
        for row in translation_reader:
            hpa_tissue = row[0]
            bto_tissue = row[1]
            if bto_tissue != '':
                db.add_db_tissue(EXTDB_HPA, bto_tissue, hpa_tissue)

def input_ion_channel_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_GO_ION_CHANNELS, 'rU') as prot_file:
        prot_reader = csv.reader(prot_file)
        next(prot_reader)
        for row in prot_reader:
            upac = row[0]
            gene_symbol = row[1]
            name = row[2]
            ions = row[3]
            gating = row[4]

            if (upac == ""):
                print("Record with empty UPAC:", row)
                continue

            if gene_symbol == '-':
                gene_symbol = ''
                
            if not db.exists_protein(upac):
                db.add_protein(
                    upac, gene_symbol=gene_symbol, name=name, ions=ions,
                    gating=gating)
            else:
                old_protein_record = db.lookup_protein(upac)
                if old_protein_record is not None:
                    old_gene_symbol = old_protein_record['GeneSymbol']
                    if old_gene_symbol == '':
                        db.update_protein_gene_symbol(upac, gene_symbol)
                    elif old_gene_symbol != gene_symbol:
                        db.update_protein_gene_symbol(upac, gene_symbol)
                        print ('Gene symbol mismatch upac=%s old_gene_symbol=%s gene_symbol=%s' \
                            % (upac, old_gene_symbol, gene_symbol))
                    old_name = old_protein_record['Name']
                    if old_name == '':
                        db.update_protein_name(upac, name)
                    db.update_protein_ions(upac, ions)
                    db.update_protein_gating(upac, gating)

def input_chembl_compound_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    proteinsUPAC = db.get_proteins_upac()
    
    #First add in batch compounds
    compoundsData = db.get_compounds()
    compoundsCHEMBLIDs = [c[5] for c in compoundsData]

    compoundsRecords = []
    with open(PATH_CHEMBL_HUMAN_ASSAYS_COMPOUND, 'rU') as compound_file:
        for line in compound_file:
            line_split = line.split('\t')
            target_upac = line_split[0].strip()

            compound_type=''
            smiles=''
            inchi=''
            compound_name = line_split[5].strip()
            compound_chembl_id = line_split[4].strip()
            iupac_name = line_split[6].strip()
            approval_status = ''
            first_approval_year = ''

            if (target_upac in proteinsUPAC and compound_chembl_id not in compoundsCHEMBLIDs):
                compoundsRecords.append([
                    COMPOUND_TYPE_COMPOUND, smiles, inchi, 
                    compound_name, compound_chembl_id, iupac_name,
                    approval_status, first_approval_year, EXTDB_CHEMBL])
                compoundsCHEMBLIDs.append(compound_chembl_id)

    db.add_multiple_compound(compoundsRecords)

    #Then add interactions in batch
    compoundsData = db.get_compounds()
    compoundsCHEMBLID = {}
    for compound in compoundsData:
        ID = compound[0]
        chemblID = compound[5]
        compoundsCHEMBLID[chemblID] = ID
    
    interactionRecords = []
    with open(PATH_CHEMBL_HUMAN_ASSAYS_COMPOUND, 'rU') as compound_file:
        for line in compound_file:
            line_split = line.split('\t')
            target_upac = line_split[0].strip()
            target_chembl_id = line_split[1].strip()
            target_type = line_split[2].strip() # Ignored
            target_name = line_split[3].strip()
            compound_chembl_id = line_split[4].strip()
            compound_name = line_split[5].strip()
            iupac_name = line_split[6].strip()
            assay_chembl_id = line_split[7].strip()
            assay_standard_type = line_split[8].strip()
            assay_standard_relation = line_split[9].strip() # Ignored
            assay_value = line_split[10].strip()
            assay_units = line_split[11].strip()
            assay_type = line_split[12].strip()
            assay_description = line_split[13].strip() # Ignored

            if assay_type not in ASSAY_TYPE_CHOICES:
                print ('Unexpected assay_type:', assay_type)

            action_type = ""

            if (target_upac not in proteinsUPAC):
                #print("Unexpected protein target UPAC:", target_upac)
                continue

            if (compound_chembl_id not in compoundsCHEMBLID):
                print("Unexpected compound CHEMBLID:", compound_chembl_id)
                continue

            compound_id = compoundsCHEMBLID[compound_chembl_id]

            #Use UPAC and compoundID to create interaction record
            interactionRecords.append([
                target_upac, compound_id, action_type, assay_value,
                assay_units, assay_standard_type, assay_type, 
                assay_chembl_id, EXTDB_CHEMBL
            ])

    db.add_multiple_interaction(interactionRecords)

def input_chembl_drug_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    proteinsUPAC = db.get_proteins_upac()
    #First add in batch compounds
    compoundsData = db.get_compounds()
    compoundsCHEMBLIDs = [c[5] for c in compoundsData]

    compoundsRecords = []
    with open(PATH_CHEMBL_HUMAN_ASSAYS_DRUG, 'rU') as drug_file:
        for line in drug_file:
            line_split = line.split('\t')
            target_upac = line_split[0].strip()

            compound_type=''
            smiles=''
            inchi=''
            compound_name = line_split[8].strip()
            compound_chembl_id = line_split[7].strip()
            dosed_compound_name = line_split[6].strip()
            approval_status = ''
            first_approval_year = ''

            if (target_upac in proteinsUPAC and compound_chembl_id not in compoundsCHEMBLIDs):
                compoundsRecords.append([
                    COMPOUND_TYPE_DRUG, smiles, inchi, 
                    compound_name, compound_chembl_id, dosed_compound_name,
                    approval_status, first_approval_year, EXTDB_CHEMBL])
                compoundsCHEMBLIDs.append(compound_chembl_id)

    db.add_multiple_compound(compoundsRecords)

    #Then add interactions in batch
    compoundsData = db.get_compounds()
    compoundsCHEMBLID = {}
    for compound in compoundsData:
        ID = compound[0]
        chemblID = compound[5]
        compoundsCHEMBLID[chemblID] = ID
    
    interactionRecords = []
    with open(PATH_CHEMBL_HUMAN_ASSAYS_DRUG, 'rU') as drug_file:
        for line in drug_file:
            line_split = line.split('\t')
            target_upac = line_split[0].strip()
            target_chembl_id = line_split[1].strip()
            target_type = line_split[2].strip() # Ignored
            target_name = line_split[3].strip()
            drug_mechanism = line_split[4].strip()
            dosed_compound_chembl_id = line_split[5].strip()
            dosed_compound_name = line_split[6].strip()
            active_compound_chembl_id = line_split[7].strip()
            active_compound_name = line_split[8].strip()
            assay_chembl_id = line_split[9].strip()
            assay_standard_type = line_split[10].strip()
            assay_standard_relation = line_split[11].strip() # Ignored
            assay_value = line_split[12].strip()
            assay_units = line_split[13].strip()
            assay_type = line_split[14].strip()
            assay_description = line_split[15].strip() # Ignored

            if drug_mechanism not in ACTION_TYPE_CHOICES:
                print ('Unexpected action_type:', drug_mechanism)
                continue

            if (target_upac not in proteinsUPAC):
                #print("Unexpected protein target UPAC:", target_upac)
                continue

            if (active_compound_chembl_id not in compoundsCHEMBLID):
                print("Unexpected compound CHEMBLID:", compound_chembl_id)
                continue

            compound_id = compoundsCHEMBLID[active_compound_chembl_id]

            #Use UPAC and compoundID to create interaction record
            interactionRecords.append([
                target_upac, compound_id, drug_mechanism, assay_value,
                assay_units, assay_standard_type, assay_type, 
                assay_chembl_id, EXTDB_CHEMBL
            ])

    db.add_multiple_interaction(interactionRecords)

def check_interaction_data():
    '''
    This function if used to check for the correctness of interaction measurement units
    SOme measurement infact have huge values compared to the other measurements. Set them as outliers
    '''
    db = icdb.IonChannelDatabase(PATH_DB)
    interactions=db.get_interactions()
    data=pd.DataFrame(interactions, columns=['interactionID', 'TargetUniprot', 'AssayValue', 'AssayType', 'CompoundChemblID'])
    data['AssayValue'].replace('', np.nan, inplace=True)
    data.dropna(subset=['AssayValue'], inplace=True)

    uniprots=list(data['TargetUniprot'].unique())
    assayTypes=list(data['AssayType'].unique())

    interactionsIDs = [] #list of interactionIDs to set zscore
    zscoresList = []
    for uniprot in uniprots:
        for assayType in assayTypes:
            values=data[(data['TargetUniprot']==uniprot) & (data['AssayType']==assayType)]
            if (len(values) > 0):
                compounds=values['CompoundChemblID'].unique()
                for compound in compounds:
                    compoundValues=values[values['CompoundChemblID']==compound]
                    if (len(compoundValues) > 2):
                        #compute zscore
                        std=compoundValues['AssayValue'].std()
                        zscores=(compoundValues['AssayValue'] - compoundValues['AssayValue'].mean()) / (std if std > 0 else 1)
                        compoundValues = compoundValues.assign(zscore=zscores.abs().values)

                        #save records
                        interactionsIDs.extend(list(compoundValues['interactionID']))
                        zscoresList.extend(list(compoundValues['zscore']))

    records = []
    for interactionID, zscore in zip(interactionsIDs, zscoresList):
        records.append([zscore, interactionID])

    db.update_multiple_interaction_zscores(records)

def input_in_betse_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    with open(PATH_IN_BETSE, 'rU') as in_betse_file:
        in_betse_reader = csv.reader(in_betse_file)
        next(in_betse_reader)
        for row in in_betse_reader:
            gene_symbol = row[0].strip()
            in_betse = row[1].strip()
            upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_symbol)
            for upac in upac_list:
                db.update_protein_in_betse(upac, in_betse)

def setup_go_term_table():
    '''
    Setup GO Terms table.
    Parse AmiGO file 

    Extract GO terms associated with Ion Channels UPAC notation
    TODO:
    -Fill table columns: name, qualifier, aspect
    -Parse and extract data from QuickGO as well
    '''
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_go_term_table()
    with open(PATH_GO_AMIGO, 'rU') as go_file:
        go_reader = csv.reader(go_file, delimiter='\t')
        next(go_reader) #skipt headers
        inserted = []
        for row in go_reader:
            bioentity = row[0]
            db_name, upac = bioentity.split(':')
            bioentity_name = row[1]
            qualifier = row[2]
            goTerm = row[3]
            annotation_extension_json = row[4]
            assigned_by = row[5]
            taxon = row[6]
            evidence_type = row[7]
            evidence_with = row[8]
            panther_family = row[9]
            type = row[10]
            bioentity_isoform = row[11]
            reference = row[12]
            date = row[13]
            if db_name == 'UniProtKB' and db.exists_protein(upac) and (upac, goTerm) not in inserted:
                db.add_go_term(upac, goTerm)
                inserted.append((upac, goTerm))

def _geneSymbolToUPACs():
    db = icdb.IonChannelDatabase(PATH_DB)
    proteinsData = db.get_proteins()
    geneSymbolToUPACs = {}
    for protein in proteinsData:
        upac = protein[0]
        geneSymbol = protein[1]
        if (geneSymbol not in geneSymbolToUPACs):
            geneSymbolToUPACs[geneSymbol] = [upac]
        else:
            geneSymbolToUPACs[geneSymbol].append(upac)
    return geneSymbolToUPACs

def _DBTissueNameToTissueName():
    db = icdb.IonChannelDatabase(PATH_DB)
    DBTissuesData = db.get_dbtissues()
    DBTissueNameToTissueName = {}
    for DBTissue in DBTissuesData:
        DBname = DBTissue[1]
        tissueName = DBTissue[2]
        DBTissueName = DBTissue[3]
        name = "{}-{}".format(DBname, DBTissueName)
        if (name in DBTissueNameToTissueName):
            print("Duplicate tissue name", name)
        DBTissueNameToTissueName[name] = tissueName
    return DBTissueNameToTissueName

def get_biogps_probeset_id_to_gene_symbol_dict():
    annot_dict = {}
    with open(PATH_BIOGPS_GNF1H_ANNOT, 'rU') as chip_annot_file:
        chip_annot_reader = csv.reader(chip_annot_file, delimiter='\t')
        next(chip_annot_reader)
        for row in chip_annot_reader:
            probeset_id = row[0]
            gene_symbol = row[6]
            if gene_symbol not in ['', 'obsoleted by Celera']:
                annot_dict[probeset_id] = gene_symbol
        chip_annot_file.close()
        return annot_dict
        
def input_biogps_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)
    annot_dict = get_biogps_probeset_id_to_gene_symbol_dict()

    geneSymbolToUPACs = _geneSymbolToUPACs()
    DBTissueNameToTissueName = _DBTissueNameToTissueName()

    expressionData = []
    with open(PATH_BIOGPS_GCRMA, 'rU') as microarray_file:
        microarray_reader = csv.reader(microarray_file)
        header = next(microarray_reader)
        tissue_list = header[1:]
        for row in microarray_reader:
            probeset_id = row[0]
            if probeset_id not in annot_dict:
                #print("ProbsetID not in annot_dict", probeset_id)
                continue

            gene_symbol = annot_dict[probeset_id].strip()
            
            if (gene_symbol not in geneSymbolToUPACs):
                #print("Gene {} not in DB".format(gene_symbol))
                continue
            upac_list = geneSymbolToUPACs[gene_symbol]
            for i in range(0, len(tissue_list)):
                DBTissueName = tissue_list[i]
                expr_level = float(row[i + 1])
                expr_level_qual = '' #quantitative values. So set Categoricals ones
                expr_units = ''

                key = "{}-{}".format(EXTDB_BIOGPS, DBTissueName)
                if (key not in DBTissueNameToTissueName):
                    print(key, "not in DBTissueNameToTissueName")
                    continue

                tissue_name = DBTissueNameToTissueName[key]
                for upac in upac_list:
                    expressionData.append([
                        tissue_name, upac, expr_level, expr_level_qual, expr_units,
                        EXPR_ASSAY_MICROARRAY, DATASET_BIOGPS, EXTDB_BIOGPS
                    ])

    db.add_multiple_expression(expressionData)

def input_hpa_rna_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)

    geneSymbolToUPACs = _geneSymbolToUPACs()
    DBTissueNameToTissueName = _DBTissueNameToTissueName()

    expressionData = []
    with open(PATH_HPA_RNA_TISSUE, 'rU') as hpa_file:
        next(hpa_file)
        for line in hpa_file:
            row = line.strip().split('\t')
            if len(row) >= 5 and row[0].strip() != '':
                gene_symbol = row[1]
                DBTissueName = row[2]
                key = "{}-{}".format(EXTDB_HPA, DBTissueName)
                if (key not in DBTissueNameToTissueName):
                    print(key, "not in DBTissueNameToTissueName")
                    continue
                tissue_name = DBTissueNameToTissueName[key]

                expr_level = float(row[3])
                expr_units = row[4]
                expr_level_qual = ''

                if (gene_symbol not in geneSymbolToUPACs):
                    #print("Gene {} not in DB".format(gene_symbol))
                    continue
                upac_list = geneSymbolToUPACs[gene_symbol]

                for upac in upac_list:
                    expressionData.append([
                        tissue_name, upac, expr_level, expr_level_qual, expr_units,
                        EXPR_ASSAY_RNASEQ, DATASET_HPA, EXTDB_HPA
                    ])

    db.add_multiple_expression(expressionData)

def input_hpa_cancer_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)

    geneSymbolToUPACs = _geneSymbolToUPACs()
    DBTissueNameToTissueName = _DBTissueNameToTissueName()

    expressionData = []
    with open(PATH_HPA_PATHOLOGY, 'rU') as hpa_file:
        next(hpa_file)
        expr_level_qual_choices = [EXPR_H, EXPR_M, EXPR_L, EXPR_ND]
        for line in hpa_file:
            row = line.strip().split('\t')
            if len(row) >= 7:
                gene_symbol = row[1]
                cancer = row[2]
                key = "{}-{}".format(EXTDB_HPA, cancer)
                if (key not in DBTissueNameToTissueName):
                    #print(key, "not in DBTissueNameToTissueName")
                    continue
                tissue_name = DBTissueNameToTissueName[key]

                if (gene_symbol not in geneSymbolToUPACs):
                    #print("Gene {} not in DB".format(gene_symbol))
                    continue
                upac_list = geneSymbolToUPACs[gene_symbol]
                for upac in upac_list:
                    if '' in row[3:7]:
                        #empty expression record.
                        continue
                    patient_count_list = [int(i) for i in row[3:7]]
                    patient_count_max_index = np.argmax(patient_count_list)
                    expr_level_qual = expr_level_qual_choices[patient_count_max_index]

                    expr_level = 0
                    expr_units = ''
                    expressionData.append([
                        tissue_name, upac, expr_level, expr_level_qual, expr_units,
                        EXPR_ASSAY_IMMUNO, DATASET_HPA, EXTDB_HPA
                    ])

    db.add_multiple_expression(expressionData)

def input_hpa_normal_expression_data():
    db = icdb.IonChannelDatabase(PATH_DB)

    geneSymbolToUPACs = _geneSymbolToUPACs()
    DBTissueNameToTissueName = _DBTissueNameToTissueName()

    expressionData = []
    with open(PATH_HPA_NORMAL, 'rU') as hpa_file:
        next(hpa_file)
        for line in hpa_file:
            row = line.strip().split('\t')
            if len(row) >= 6:
                gene_symbol = row[1]
                tissue = row[2]
                key = "{}-{}".format(EXTDB_HPA, tissue)
                if (key not in DBTissueNameToTissueName):
                    #print(key, "not in DBTissueNameToTissueName")
                    continue
                tissue_name = DBTissueNameToTissueName[key]

                expr_level_qual = row[4]
                expr_level = 0
                expr_units = ''

                if expr_level_qual not in [EXPR_H, EXPR_M, EXPR_L, EXPR_ND]:
                    print("expr_level_qual not in [EXPR_H, EXPR_M, EXPR_L, EXPR_ND]")
                    continue

                if (gene_symbol not in geneSymbolToUPACs):
                    #print("Gene {} not in DB".format(gene_symbol))
                    continue
                upac_list = geneSymbolToUPACs[gene_symbol]
                for upac in upac_list:
                    expressionData.append([
                        tissue_name, upac, expr_level, expr_level_qual, expr_units,
                        EXPR_ASSAY_IMMUNO, DATASET_HPA, EXTDB_HPA
                    ])

    db.add_multiple_expression(expressionData)

# NOT USED
def setup_genbank_uniprot_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_genbank_uniprot_table()
    with open(PATH_UNIPROT_GENBANK, 'rU') as gu_file:
        next(gu_file)
        for line in gu_file:
            line_split = line.split('\t')
            if len(line_split) == 2:
                gbac = line_split[0].strip()
                upac = line_split[1].strip()
                if upac != '' and not db.exists_genbank_uniprot(gbac, upac):
                    db.add_genbank_uniprot(gbac, upac)

# NOT USED
def get_ensembl_to_uniprot_dict():
    ensembl_to_uniprot = {}
    with open(PATH_UNIPROT_ENSEMBL_IDMAPPING, 'rU') as map_file:
        for line in map_file:
            row = line.strip().split('\t')
            upac = row[0]
            ensg = row[2]
            ensembl_to_uniprot[ensg] = upac
        map_file.close()
        return ensembl_to_uniprot


def get_channelpedia_dict():
    '''
    Helper function for setup_channel_class_tables()
    '''
    channelpedia_dict = {}
    with open(PATH_CHANNELPEDIA_INFO, 'rU') as channelpedia_file:
        channelpedia_reader = csv.reader(channelpedia_file)
        next(channelpedia_reader)
        for row in channelpedia_reader:
            subclass = row[1]
            intro_text = row[2]
            url = row[3]
            channelpedia_dict[subclass] = (intro_text, url)
    return channelpedia_dict


def setup_channel_class_tables():
    '''
    Setup ChannelPedia tables: store info relative to ion channels with info from channelpedia
    1) Build helper dictionary {channelIdentifier: (description, url)}
    2) Parse Proteins-Channel class data to populate tables:
        -superclass, class and subclass
    3) Update proteins table column: IonChannelSubClass with channel subclass
    '''
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_channel_super_class_table()
    db.create_channel_class_table()
    db.create_channel_sub_class_table()
    channelpedia_dict = get_channelpedia_dict()
    with open(PATH_CHANNEL_CLASSES, 'rU') as class_file:
        class_reader = csv.reader(class_file)
        next(class_reader)
        for row in class_reader:
            channel_superclass = row[0]
            channel_class = row[1]
            channel_subclass = row[2]
            subfamily = row[3]
            gene_symbol = row[4]
            if not db.exists_channel_super_class(channel_superclass):
                db.add_channel_super_class(channel_superclass)
            if not db.exists_channel_class(channel_class):
                db.add_channel_class(channel_class, channel_superclass)
            if not db.exists_channel_sub_class(channel_subclass):
                if channel_subclass in channelpedia_dict:
                    channelpedia_text, channelpedia_url = channelpedia_dict[
                        channel_subclass]
                else:
                    channelpedia_text = ''
                    channelpedia_url = ''
                db.add_channel_sub_class(
                    channel_subclass, channel_class,
                    channelpedia_text=channelpedia_text,
                    channelpedia_url=channelpedia_url, subfamily=subfamily)
            upac_list = db.get_uniprot_accnums_by_gene_symbol(gene_symbol)
            for upac in upac_list:
                db.update_protein_sub_class(upac, channel_subclass)

def setup_specificity_table_with_jp_method():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_specificity_table()
    upac_list = db.get_protein_uniprot_accnums()
    tissue_name_list = db.get_tissue_names()
    specificity_score_threshold = 5.
    specificity_score_dict = {}
    expr_table = {}

    #for each protein
    for upac in upac_list:
        specificity_score = 0.
        tissue_count = 0
        #for each protein-tissue pair
        for tissue_name in tissue_name_list:
            #get expression level HPA DB
            expr_level_list = \
                db.get_expression_level_by_uniprot_accnum_tissue_dataset(
                    upac, tissue_name, EXTDB_HPA)
            if expr_level_list:
                #If available save it. Increase score only if above thres
                expr_level_rna = expr_level_list[0]
                expr_table[(tissue_name, upac)] = expr_level_rna
                if expr_level_rna > specificity_score_threshold:
                    specificity_score += 1.
                tissue_count += 1
        if tissue_count > 0:
            specificity_score = 1. - specificity_score / float(tissue_count)
        else:
            specificity_score = None
        specificity_score_dict[upac] = specificity_score
    for tissue_name, upac in expr_table:
        expr_level_rna = expr_table[(tissue_name, upac)]
        specificity_score_for_tissue = \
            specificity_score_dict[upac] * expr_level_rna
        if specificity_score_for_tissue != None:
            db.add_specificity(tissue_name, upac, specificity_score_for_tissue)


def vacuum_db():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.vacuum_db()


def dump_db():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.dump_db()


def print_db_stats():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.print_db_stats()

#FUNCTIONS CALLED DURING SETUP. IN FOLLOWING ORDER
def setup_external_db_table():
    _setup_external_db_table()

def setup_tissue_table():
    _setup_tissue_table()

def setup_db_tissue_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_db_tissue_table()
    print("input biogps tissue data..")
    input_biogps_tissue_data()
    print("input hpa tissue data..")
    input_hpa_tissue_data()

def setup_protein_compound_interaction_tables():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_protein_table()
    db.create_compound_table()
    db.create_interaction_table()
    print("input ion channel data..")
    input_ion_channel_data()
    print("input chembl compound data..")
    input_chembl_compound_data()
    print("input chembl drug data..")
    input_chembl_drug_data()
    print("check interaction data..")
    check_interaction_data()

def setup_expression_table():
    db = icdb.IonChannelDatabase(PATH_DB)
    db.create_expression_table()
    print("input BioGPS expression data..")
    input_biogps_expression_data()
    print("input HPA RNA expression data..")
    input_hpa_rna_expression_data()
    print("input HPA cancer expression data..")
    input_hpa_cancer_expression_data()
    print("input HPA normal expression data..")
    input_hpa_normal_expression_data()

if __name__=='__main__':
    check_interaction_data()