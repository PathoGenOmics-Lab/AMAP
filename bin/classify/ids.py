"""
Gene ID cleaning and formatting utilities.

Handles the AMAP naming convention where gene IDs may include:
- Sample prefix: L4_G12345_
- Gene name: Rv0001, PE_PGRS1, dnaA, etc.
- Embedded coordinates: _703830-704057

These utilities normalize IDs for consistent downstream processing.
"""

import re


def clean_id(gene_id: str) -> str:
    """
    Strip trailing embedded coordinates from a protein/gene ID.
    
    Reference proteins have format GENE_ID_START-END (e.g., Rv0609A_703830-704057).
    When the protein is from a typed sample, the full ID is:
    L8_G26689_Rv0609A_703830-704057 (sample prefix + gene + coords).
    
    This function removes the trailing _DIGITS-DIGITS suffix if present.
    It's safe for IDs without coordinates (e.g., L8_G26689_00042 → unchanged).
    
    Old behavior: stripped the last underscore part only when exactly 4 parts
    existed. This was brittle and failed for PE_PGRS genes (5+ parts).
    """
    gene_id = str(gene_id)
    # Only strip trailing _DIGITS-DIGITS (embedded coordinates)
    return re.sub(r'_\d+-\d+$', '', gene_id)


def format_cluster_id(gene_id: str, base_name: str) -> str:
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
