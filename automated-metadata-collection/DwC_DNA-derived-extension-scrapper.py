import requests
from bs4 import BeautifulSoup
import pandas as pd

url = 'http://rs.gbif.org/extension/gbif/1.0/dna_derived_data_2021-07-05.xml'

response = requests.get(url)
soup = BeautifulSoup(response.content, 'xml')

namespace = {'dc': 'http://purl.org/dc/terms/'}

elements = soup.find_all('property')
data = []

for element in elements:
    name = element.get('name', 'N/A')
    qualName = element.get('qualName', 'N/A')
    uri = qualName.split('/')[-1]
    group = element.get('group', 'N/A')
    description = element.get('dc:description', 'N/A')
    required = element.get('required', 'N/A')
    data.append([name, uri, group, required, description])

columns = ['field_name', 'uri', 'group', 'required', 'description']

df = pd.DataFrame(data, columns=columns)

df.to_csv('DwC_DNA-derived-extension.csv', index=False)