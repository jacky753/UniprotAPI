"""
import requests

def get_interpro_accessions_from_uniprot(uniprot_id: str):
    
    #UniProt の REST API を使って、指定された UniProt ID に紐づく InterPro エントリーの accession を取得。
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    resp = requests.get(url)
    if resp.status_code != 200:
        raise RuntimeError(f"UniProt API error: status {resp.status_code} for ID {uniprot_id}")
    data = resp.json()
    interpro_accessions = set()
    for feature in data.get("features", []):
        if feature.get("type") == "Domain" and "interpro" in feature:
            interpro_accessions.add(feature["interpro"]["id"])
    return list(interpro_accessions)
"""

import requests

def get_interpro_accessions_from_uniprot(uniprot_id: str):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    resp = requests.get(url)
    if resp.status_code != 200:
        raise RuntimeError(f"UniProt API error: status {resp.status_code} for ID {uniprot_id}")
    data = resp.json()

    interpro_ids = set()
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") == "InterPro" and "id" in ref:
            interpro_ids.add(ref["id"])

    return list(interpro_ids)



def get_interpro_entry(interpro_id: str):
    """
    InterPro API を使って、指定された InterPro エントリー ID の詳細を取得。
    """
    url = f"https://www.ebi.ac.uk/interpro/api/entry/InterPro/{interpro_id}"
    resp = requests.get(url, params={"format": "json"})
    if resp.status_code != 200:
        raise RuntimeError(f"InterPro API error: status {resp.status_code} for entry {interpro_id}")
    return resp.json()

def fetch_domain_structure_for_uniprot(uniprot_id: str):
    """
    UniProt ID -> 関連する InterPro エントリー -> 各エントリーの詳細を取得。
    """
    interpro_ids = get_interpro_accessions_from_uniprot(uniprot_id)
    results = {}
    for ipr in interpro_ids:
        entry = get_interpro_entry(ipr)
        results[ipr] = entry
    return results

if __name__ == "__main__":
    uniprot_id = "P69905"  # 例
    try:
        interpro_entries = fetch_domain_structure_for_uniprot(uniprot_id)
        print(f"UniProt ID: {uniprot_id}")
        for ipr, entry in interpro_entries.items():
            print("ipr: ")
            print(ipr)
            print("entry: ")
            print(entry)

            print(f"InterPro ID: {ipr}")
            # entry の中から、domain の名称・説明・関連ドメイン構造等を表示
            print("  Name:", entry.get("metadata", {}).get("name"))
            print("  Description:", entry.get("metadata", {}).get("short_description"))
            # 必要に応じて「構造」など他のフィールドも
    except Exception as e:
        print("Error:", e)
