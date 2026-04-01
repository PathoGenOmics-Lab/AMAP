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
    Formats the final gene ID with consistent sample prefix.
    
    Rules:
    - Strip embedded coordinates (trailing _DIGITS-DIGITS)
    - Ensure result always has the sample base_name prefix
    - Preserve multi-part gene names (e.g., PE_PGRS1)
    
    Examples:
        ('Rv0001', 'L8_G26689')           → 'L8_G26689_Rv0001'
        ('Rv0001_1-1524', 'L8_G26689')    → 'L8_G26689_Rv0001'
        ('PE_PGRS1', 'L8_G26689')         → 'L8_G26689_PE_PGRS1'
        ('L8_G26689_Rv0001', 'L8_G26689') → 'L8_G26689_Rv0001'
        ('L8_G26689_00042', 'L8_G26689')  → 'L8_G26689_00042'
    """
    import re
    gene_id = str(gene_id)

    # 1. Strip trailing embedded coordinates (_DIGITS-DIGITS)
    gene_id = re.sub(r'_\d+-\d+$', '', gene_id)

    # 2. If already starts with base_name, return as-is
    if gene_id.startswith(base_name):
        return gene_id

    # 3. If it has 3+ underscore parts and the first two look like a sample prefix
    #    (e.g., L8_G26689_Rv0001), keep first 3
    parts = gene_id.split('_')
    if len(parts) >= 3:
        return '_'.join(parts[:3])

    # 4. Otherwise prepend base_name (handles 'Rv0001', 'PE_PGRS1', etc.)
    return f"{base_name}_{gene_id}"
