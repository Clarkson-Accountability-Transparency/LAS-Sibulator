import argparse
import csv
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)

def convert_csv(inFilename, outFilename): 
    logger.info(f'converting {inFilename} to {outFilename}')
    with open(inFilename, 'r') as csvF, open(outFilename, 'w', newline='') as csvOut:
        csvReader = csv.reader(csvF)
        csvWriter = csv.writer(csvOut)
        headers = ['Sample']
        rowData = []
        readHeaders = next(csvReader)
        for row in csvReader:
            if len(rowData) == 0:
                rowData.append(row[0]) #set the Sample
            rowData.extend( (row[2], row[3]) )
            if row[1] == 'VWA':
                #KLUDGE to fix the naming discrepancy
                headers.append('vWA')
                headers.append('vWA')
            else:
                headers.append(row[1])
                headers.append(row[1])
        csvWriter.writerow(headers)
        csvWriter.writerow(rowData)
    
if __name__ == '__main__':
    argParser =argparse.ArgumentParser()
    argParser.add_argument('inputcsv', help='File to convert')
    args = argParser.parse_args()

    convert_csv(args.inputcsv, f'{args.inputcsv}_sibs.csv')
