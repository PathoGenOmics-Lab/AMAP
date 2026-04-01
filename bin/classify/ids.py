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
    - 2 parts where second looks like embedded coords (e.g., 'Rv0001_1-1524')
      → 'Sample1_Rv0001' (strip the coordinates)
    - 3+ parts → keep first 3
    - Otherwise → return as-is
    """
    import re
    parts = str(gene_id).split('_')
    if len(parts) == 1:
        return f"{base_name}_{parts[0]}"
    elif len(parts) == 2:
        # Check if second part is embedded coordinates (digits-digits)
        if re.match(r'^\d+-\d+$', parts[1]):
            return f"{base_name}_{parts[0]}"
        return gene_id
    elif len(parts) >= 3:
        return '_'.join(parts[:3])
    return gene_id
