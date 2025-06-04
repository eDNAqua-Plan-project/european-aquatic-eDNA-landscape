from Bio import Entrez
import csv
import time

Entrez.email = "camila.babo@cibio.up.pt"

def get_pmid_from_doi(doi):
    """
    Search PubMed for DOI to get PMID
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=doi, field="doi")
        record = Entrez.read(handle)
        handle.close()
        pmid = record["IdList"][0] if record["IdList"] else None
        return pmid
    except Exception as e:
        print(f"Error fetching PMID for DOI {doi}: {e}")
        return None

def get_metadata(pmid):
    """
    Fetch metadata using PubMed EFetch
    """
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]

        # Get year
        pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        year = pub_date.get("Year", "N/A")

        # Get keywords
        keywords = records["PubmedArticle"][0]["MedlineCitation"].get("KeywordList", [])
        keywords_flat = [kw for sublist in keywords for kw in sublist] if keywords else []
        keywords_str = "; ".join(keywords_flat) if keywords_flat else "N/A"

        return year, keywords_str

    except Exception as e:
        print(f"Error fetching metadata for PMID {pmid}: {e}")
        return "N/A", "N/A"

def main():
    input_file = "LLM_DOIs.txt"
    output_file = "output.csv"

    with open(input_file, "r") as f:
        dois = [line.strip() for line in f if line.strip()]

    total = len(dois)

    with open(output_file, "w", newline='', encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["DOI", "Year", "Keywords"])

        for i, doi in enumerate(dois, start=1):
            pmid = get_pmid_from_doi(doi)
            if not pmid:
                writer.writerow([doi, "N/A", "N/A"])
                print(f"[{i}/{total}] DOI not found: {doi}")
                continue

            year, keywords = get_metadata(pmid)
            writer.writerow([doi, year, keywords])
            print(f"[{i}/{total}] Retrieved: {doi} (PMID: {pmid})")
            time.sleep(1)  # Be polite to NCBI's servers

    print(f"\n Done! Mined metadata for {total} DOIs. Output saved to '{output_file}'.")


if __name__ == "__main__":
    main()
