import requests

base_url = "https://api.gbif.org/v1/occurrence"
response = requests.get(f"{base_url}&offset=0&limit=1")  # Request a small amount of data
data = response.json()

total_count = data.get('count', 0)  # Get the total number of records
print(total_count)

