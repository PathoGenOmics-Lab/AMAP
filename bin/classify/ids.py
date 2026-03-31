"""Gene ID cleaning and formatting utilities."""


def clean_id(gene_id):
    """
    Cleans a gene ID based on the 4-part rule.
    If the ID split by '_' has 4 elements, it removes the last one.
    Example: 'L8_G26689_Rv0609A_703830-704057' becomes 'L8_G26689_Rv0609A'.
    """
    gene_id = str(gene_id)
    parts = gene_id.split('_')
    if len(parts) == 4:
        return '_'.join(parts[:-1])
    return gene_id


def format_cluster_id(gene_id, base_name):
    """
    Formats the final gene ID.
    - Single part (e.g., 'Rv1234') → 'Sample1_Rv1234'
    - 3+ parts → keep first 3
    - Otherwise → return as-is
    """
    parts = str(gene_id).split('_')
    if len(parts) == 1:
        return f"{base_name}_{parts[0]}"
    elif len(parts) >= 3:
        return '_'.join(parts[:3])
    return gene_id
