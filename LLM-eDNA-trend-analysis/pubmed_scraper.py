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

        article_data = records["PubmedArticle"][0]
        article = article_data["MedlineCitation"]["Article"]

        # Get year
        pub_date = article.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        year = pub_date.get("Year", "N/A")

        # Get keywords
        keywords = article_data["MedlineCitation"].get("KeywordList", [])
        keywords_flat = [kw for sublist in keywords for kw in sublist] if keywords else []
        keywords_str = "; ".join(keywords_flat) if keywords_flat else "N/A"

        # Get MeSH terms
        mesh_headings = article_data["MedlineCitation"].get("MeshHeadingList", [])
        mesh_terms = []
        major_topics = []
        subheadings = []

        for mesh in mesh_headings:
            descriptor = mesh["DescriptorName"]
            descriptor_name = descriptor.title()
            mesh_terms.append(descriptor_name)

            if descriptor.attributes.get("MajorTopicYN") == "Y":
                major_topics.append(descriptor_name)

            for sub in mesh.get("QualifierName", []):
                sub_name = sub.title()
                subheadings.append(sub_name)

        mesh_terms_str = "; ".join(mesh_terms) if mesh_terms else "N/A"
        major_topics_str = "; ".join(major_topics) if major_topics else "N/A"
        subheadings_str = "; ".join(subheadings) if subheadings else "N/A"

        return year, keywords_str, mesh_terms_str, major_topics_str, subheadings_str

    except Exception as e:
        print(f"Error fetching metadata for PMID {pmid}: {e}")
        return "N/A", "N/A", "N/A", "N/A", "N/A"

def main():
    input_file = "LLM-eDNA-trend-analysis/LLM_DOIs.txt"
    output_file = "LLM-eDNA-trend-analysis/output-mesh-terms.csvrt  "

    with open(input_file, "r") as f:
        dois = [line.strip() for line in f if line.strip()]

    total = len(dois)

    with open(output_file, "w", newline='', encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["DOI", "Year", "Keywords", "MeSH Terms", "Major Topics", "Subheadings"])

        for i, doi in enumerate(dois, start=1):
            pmid = get_pmid_from_doi(doi)
            if not pmid:
                writer.writerow([doi, "N/A", "N/A", "N/A", "N/A", "N/A"])
                print(f"[{i}/{total}] DOI not found: {doi}")
                continue

            year, keywords, mesh_terms, major_topics, subheadings = get_metadata(pmid)
            writer.writerow([doi, year, keywords, mesh_terms, major_topics, subheadings])
            print(f"[{i}/{total}] Retrieved: {doi} (PMID: {pmid})")
            time.sleep(1)  # Be polite to NCBI's servers

    print(f"\n Mined metadata for {total} DOIs. Output saved to '{output_file}'.")

if __name__ == "__main__":
    main()
