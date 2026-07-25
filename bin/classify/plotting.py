"""Visualization: interactive overlap plots."""

import plotly.graph_objects as go


def generate_overlap_plot(overlap_groups, output_path):
    """Generate an interactive HTML plot of gene overlaps."""
    fig = go.Figure()
    y_pos = 0

    for key, group in overlap_groups.items():
        prot_id = group['Prot Id']
        locus_tags = group['Locus Tags']

        # Bakta CDS (blue)
        case4_x, case4_y, case4_text = [], [], []
        for idx, (start, end) in enumerate(group['case4']):
            tag = locus_tags[idx] if idx < len(locus_tags) else "Unknown"
            case4_x.extend([start, end, None])
            case4_y.extend([y_pos, y_pos, None])
            case4_text.extend([tag, tag, None])

        fig.add_trace(go.Scatter(
            x=case4_x, y=case4_y, mode='lines+markers',
            line=dict(color='blue', width=6), marker=dict(size=10),
            name=f"Bakta: {locus_tags[0] if locus_tags else ''}",
            hoverinfo='text', text=case4_text))

        # Miniprot alignment (red)
        mp_x, mp_y, mp_text = [], [], []
        for start, end in group['miniprot_df']:
            mp_x.extend([start, end, None])
            mp_y.extend([y_pos, y_pos, None])
            mp_text.extend([prot_id, prot_id, None])

        fig.add_trace(go.Scatter(
            x=mp_x, y=mp_y, mode='lines+markers',
            line=dict(color='red', width=10), marker=dict(size=14),
            name=f"Miniprot: {prot_id}", hoverinfo='text', text=mp_text, opacity=0.7))

        y_pos += 1

    fig.update_layout(
        title='Visualization of Overlapping Gene Annotations',
        xaxis_title='Genomic Coordinate', yaxis_title='Overlap Group',
        showlegend=False)
    fig.write_html(output_path)
