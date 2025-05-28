import requests
import json
import time
from requests.exceptions import RequestException, JSONDecodeError

def fetch_data_with_limit(base_url, max_records, max_retries=3, retry_delay=5):
    all_results = []
    offset = 0
    limit = 300  #arbitrary

    while len(all_results) < max_records:
        url = f"{base_url}&offset={offset}&limit={limit}" # Add offset and limit to the URL to control pagination

        for attempt in range(max_retries):
            try:
                response = requests.get(url, timeout=30)
                response.raise_for_status()
                data = response.json()
                break  # If successful, break out of the retry loop
            except (RequestException, JSONDecodeError) as e:
                if attempt < max_retries - 1:
                    print(f"Error occurred: {e}. Retrying in {retry_delay} seconds...")
                    time.sleep(retry_delay)
                else:
                    print(f"Failed to fetch data after {max_retries} attempts. Last error: {e}")
                    return all_results

        results = data.get('results', []) # Get the 'results' list from the response data and default to an empty list
        all_results.extend(results[:max_records - len(all_results)]) # Append the results to the all_results list and limit the total to max_records

        if len(results) < limit or offset + limit >= data.get('count', 0): # checks if there are fewer results than the limit (no more pages with data) or if the next batch would exceed the total number of available records
            break

        offset += limit # Increment the offset to get the next page of results
        print(f"Fetched {len(all_results)} records so far...")

    return all_results[:max_records]

def count_entries_with_extra_keys(data):
    count = 0
    standard_keys = {'http://rs.gbif.org/terms/dna_sequence', 'https://w3id.org/mixs/0000044'}

    for entry in data['occurrences']:
        for dna_data in entry['DNADerivedData']:
            if set(dna_data.keys()) > standard_keys: # set comparison (>) checks if the left set is
                # a superset of the right set, meaning all elements of B are in A, and A has at
                # least one additional element not in B.
                count += 1
                break

    return count

base_url = "https://api.gbif.org/v1/occurrence/search?gbifRegion=EUROPE&habitat=aquatic&basis_of_record=MATERIAL_SAMPLE&extension=dna-derived-data&fields=extensions,gbifID"

max_records = 100000
all_data = fetch_data_with_limit(base_url, max_records)

dna_extension_key = 'http://rs.gbif.org/terms/1.0/DNADerivedData'

# Filter out occurrences with DNA-derived data to write to the output file
matching_occurrences = []

for result in all_data:
    extensions = result.get('extensions', {})
    dna_data = extensions.get(dna_extension_key, [])
    if dna_data:
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

# Count entries with extra keys
extra_keys_count = count_entries_with_extra_keys(output_data)
print(f"Number of entries with additional keys: {extra_keys_count}")