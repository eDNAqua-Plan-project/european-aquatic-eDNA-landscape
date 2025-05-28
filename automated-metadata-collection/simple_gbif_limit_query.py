import requests
import json

url = "https://api.gbif.org/v1/occurrence/search?gbifRegion=EUROPE&habitat=aquatic&basis_of_record=MATERIAL_SAMPLE&extension=dna-derived-data&limit=300&fields=extensions,gbifID"

response = requests.get(url)
data = response.json()

dna_extension_key = 'http://rs.gbif.org/terms/1.0/DNADerivedData'
matching_occurrences = []

for result in data['results']:
    extensions = result.get('extensions', {})
    dna_data = extensions.get(dna_extension_key, [])

    if dna_data:  # If DNA-derived data is present, add it to matching_occurrences
        matching_occurrences.append({
            'gbifID': result.get('gbifID'),
            'DNADerivedData': dna_data
        })

output_data = {
    'count': len(matching_occurrences),
    'occurrences': matching_occurrences
}

# Write the output data to a JSON file
with open('dna_derived_data.json', 'w') as file:
    json.dump(output_data, file, indent=2)

print(f"Data saved to dna_derived_data.json with {len(matching_occurrences)} occurrences.")
