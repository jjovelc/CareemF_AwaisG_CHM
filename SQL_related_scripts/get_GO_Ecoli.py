import requests

def get_uniprot_ids_for_ecoli(gene_names):
    url = "https://rest.uniprot.org/idmapping/run"
    params = {
        'from': 'Gene_Name',
        'to': 'UniProtKB_AC-ID',
        'taxonId': '83333'  # Taxonomy ID for Escherichia coli (strain K12)
    }
    data = {
        'ids': ' '.join(gene_names)
    }
    
    # Submit the ID mapping job
    response = requests.post(url, params=params, data=data)
    
    if response.status_code != 200:
        print(f"Error: Received status code {response.status_code}")
        print(response.text)  # Print the response text for debugging
        return []
    
    job_id = response.json()['jobId']
    
    # Check the status of the job
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status_response = requests.get(status_url)
        status = status_response.json()
        if status['jobStatus'] == 'FINISHED':
            break
    
    # Get the results
    result_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}"
    result_response = requests.get(result_url)
    
    lines = result_response.text.strip().split('\n')
    if len(lines) < 2:
        print("Error: No results returned from UniProt")
        return []
    
    uniprot_ids = []
    for line in lines[1:]:
        parts = line.split('\t')
        if len(parts) < 2:
            print(f"Error: Unexpected line format: {line}")
        else:
            uniprot_ids.append(parts[1])
    
    return uniprot_ids

def get_go_annotations(uniprot_ids):
    url = "https://www.uniprot.org/uniprot/"
    annotations = {}
    for uniprot_id in uniprot_ids:
        response = requests.get(f"{url}{uniprot_id}.txt")
        if response.status_code != 200:
            print(f"Error: Failed to retrieve data for {uniprot_id}, status code {response.status_code}")
            continue
        data = response.text
        go_terms = [line for line in data.split('\n') if line.startswith('DR   GO;')]
        annotations[uniprot_id] = go_terms
    return annotations

# Example usage
gene_names = ['clpC', 'cusC']
uniprot_ids = get_uniprot_ids_for_ecoli(gene_names)
print(uniprot_ids)  # Print the UniProt IDs to verify if they were retrieved correctly

if uniprot_ids:
    go_annotations = get_go_annotations(uniprot_ids)

    for uniprot_id, go_terms in go_annotations.items():
        print(f"UniProt ID: {uniprot_id}")
        for term in go_terms:
            print(term)