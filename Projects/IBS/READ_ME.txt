Diagnosis RomeIII questionnaire

# Version 1: 21 July 2016 Ettje Tigchelaar
# Update on 25 march 2020 by Trishla Sinha: translating, editing scripts as per new requirements 
# Update on 30 june 2020 by Trishla Sinha: editing scripts as per new requirements for all LL participants 

# accompanying document: RomeIII_functional bowel disorders questionnaire.pdf
# execution of the code below will results in 29 additional variables:
  # - 26 scoring variables representing the 26 questions of the questionnaire
  # - variable IBS stating yes or no for IBS
  # - variable Subtype specifying the IBS subtype C/D/M/U
  # - variable FGID listing all functional gastrointestinal disorders from this questionnaire: IBS, IBS strict, functional constipation, functional diarrhea, functional bloating

# name of dataframe: RO 
# If containing no sex information: Add a column sex to this data R0$sex with labels  male and female
# Variable names questions: ROME1 to ROME26 


## score questions based on numeric values in LL catalogue ( https://catalogue.lifelines.nl/menu/main/protocolviewer) 
