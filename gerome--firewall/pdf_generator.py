import io
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors

def generate_clinical_pdf(analysis, report, infection_site="Systemic / General"):
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(
        buffer,
        pagesize=letter,
        rightMargin=36,
        leftMargin=36,
        topMargin=36,
        bottomMargin=36
    )
    story = []
    styles = getSampleStyleSheet()

    # Document Typography Styles
    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=15,
        leading=18,
        textColor=colors.HexColor('#0F172A')
    )
    subtitle_style = ParagraphStyle(
        'SubTitle',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8,
        leading=11,
        textColor=colors.HexColor('#475569')
    )
    h2_style = ParagraphStyle(
        'H2',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=10,
        leading=13,
        textColor=colors.HexColor('#1E293B'),
        spaceBefore=8,
        spaceAfter=4
    )
    body_style = ParagraphStyle(
        'Body',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=8,
        leading=10,
        textColor=colors.HexColor('#334155')
    )
    header_cell_style = ParagraphStyle(
        'HeaderCell',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=8,
        leading=10,
        textColor=colors.white
    )

    # Title & Metadata
    story.append(Paragraph("GENOME FIREWALL — CLINICAL DIAGNOSTIC REPORT", title_style))
    story.append(Paragraph(f"ANTIMICROBIAL RESISTANCE INFERENCE & PK/PD ASSESSMENT // SITE: {infection_site.upper()}", subtitle_style))
    story.append(Spacer(1, 8))
    story.append(HRFlowable(width="100%", thickness=1, color=colors.HexColor('#0F172A'), spaceBefore=0, spaceAfter=8))

    meta_data = [
        [
            Paragraph(f"<b>Accession ID:</b> {analysis.get('accession_id', 'N/A')}", body_style),
            Paragraph(f"<b>Organism:</b> {analysis.get('species', 'N/A')}", body_style)
        ],
        [
            Paragraph(f"<b>Assembly Length:</b> {analysis.get('total_bp', 0):,} bp", body_style),
            Paragraph(f"<b>Sequence Type:</b> {analysis.get('assembly_type', 'N/A')}", body_style)
        ]
    ]
    t_meta = Table(meta_data, colWidths=[270, 270])
    t_meta.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#F8FAFC')),
        ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
        ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#E2E8F0')),
        ('PADDING', (0, 0), (-1, -1), 5),
    ]))
    story.append(t_meta)
    story.append(Spacer(1, 10))

    # Resistance Determinants Table
    story.append(Paragraph("DETECTED RESISTANCE DETERMINANTS", h2_style))
    genes = analysis.get("detected_genes", [])
    if genes:
        gene_table_data = [[
            Paragraph("Gene Marker", header_cell_style),
            Paragraph("Target Class", header_cell_style),
            Paragraph("Identity %", header_cell_style),
            Paragraph("Coverage %", header_cell_style),
            Paragraph("Database", header_cell_style)
        ]]
        for g in genes:
            gene_table_data.append([
                Paragraph(str(g.get("name", "")), body_style),
                Paragraph(str(g.get("drug_class", "")), body_style),
                Paragraph(f"{g.get('identity_pct', 0):.1f}%", body_style),
                Paragraph(f"{g.get('coverage_pct', 0):.1f}%", body_style),
                Paragraph(str(g.get("ref_db", "")), body_style)
            ])
        t_gene = Table(gene_table_data, colWidths=[110, 180, 80, 80, 90])
        t_gene.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1E293B')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#E2E8F0')),
            ('PADDING', (0, 0), (-1, -1), 4),
        ]))
        story.append(t_gene)
    else:
        story.append(Paragraph("No primary resistance determinants identified.", body_style))

    story.append(Spacer(1, 10))

    # Susceptibility Matrix Table
    story.append(Paragraph("ANTIMICROBIAL SUSCEPTIBILITY & PK/PD EVALUATION MATRIX", h2_style))
    matrix_data = [[
        Paragraph("Antibiotic", header_cell_style),
        Paragraph("Prediction", header_cell_style),
        Paragraph("Primary Marker", header_cell_style),
        Paragraph("PK/PD Site Context", header_cell_style)
    ]]
    for r in report:
        matrix_data.append([
            Paragraph(str(r.get("Antibiotic", "")), body_style),
            Paragraph(str(r.get("Prediction", "")), body_style),
            Paragraph(str(r.get("Primary Determinant", "None")), body_style),
            Paragraph(str(r.get("PKPD_Notes", "Standard Spectrum")), body_style)
        ])
    t_matrix = Table(matrix_data, colWidths=[110, 100, 130, 200])
    t_matrix.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0F172A')),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
        ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#CBD5E1')),
        ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#E2E8F0')),
        ('PADDING', (0, 0), (-1, -1), 4),
    ]))
    story.append(t_matrix)

    story.append(Spacer(1, 14))
    story.append(HRFlowable(width="100%", thickness=0.5, color=colors.HexColor('#94A3B8'), spaceBefore=4, spaceAfter=6))
    story.append(Paragraph("CONFIDENTIAL CLINICAL DECISION SUPPORT SYSTEM // NOT FOR STANDALONE DIAGNOSTIC USE", subtitle_style))

    doc.build(story)
    buffer.seek(0)
    return buffer.getvalue()