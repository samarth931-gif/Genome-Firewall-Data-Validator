# predictor.py
import json

SPECTRUM_RULES = {
    "Methicillin": {"target_gram": "Gram-Positive", "reason": "Methicillin is intrinsically ineffective against Gram-negative species (outer membrane barrier)."},
    "Vancomycin": {"target_gram": "Gram-Positive", "reason": "Large glycopeptide molecules cannot cross Gram-negative outer membrane channels."}
}

DYNAMIC_MECHANISMS = {
    ("Amoxicillin", "blaKPC-2"): "Hydrolyzed by KPC-2 serine carbapenemase.",
    ("Amoxicillin", "blaNDM-1"): "Hydrolyzed by NDM-1 metallo-beta-lactamase.",
    ("Amoxicillin", "blaCTX-M-15"): "Hydrolyzed by CTX-M-15 extended-spectrum beta-lactamase (ESBL).",
    ("Ceftriaxone", "blaCTX-M-15"): "Hydrolyzed by CTX-M-15 extended-spectrum beta-lactamase (ESBL).",
    ("Ceftriaxone", "blaNDM-1"): "Hydrolyzed by NDM-1 metallo-beta-lactamase.",
    ("Imipenem", "blaKPC-2"): "Carbapenem ring hydrolyzed by KPC-2 serine carbapenemase.",
    ("Imipenem", "blaNDM-1"): "Zinc-dependent carbapenem ring cleavage by NDM-1 metallo-beta-lactamase.",
    ("Gentamicin", "aac(6')-Ib"): "Covalently modified by aac(6')-Ib acetyltransferase.",
    ("Gentamicin", "armA"): "16S rRNA G1405 methylation by ArmA prevents Gentamicin ribosomal binding.",
    ("Ciprofloxacin", "gyrA_S83L"): "Ser83Leu mutation in DNA Gyrase subunit A disrupts quinolone binding pocket.",
    ("Colistin", "mcr-1"): "Phosphoethanolamine transferase modifies Lipid A, neutralizing LPS negative charge."
}

RESCUE_CANDIDATES = {
    "blaKPC-2": [
        {"drug": "Ceftazidime-Avibactam", "class": "Cephalosporin + BLI", "mechanism": "Avibactam neutralizes KPC-2 serine carbapenemase, protecting Ceftazidime."},
        {"drug": "Meropenem-Vaborbactam", "class": "Carbapenem + BLI", "mechanism": "Vaborbactam boronic acid inhibitor neutralizes Class A KPC enzymes."}
    ],
    "blaNDM-1": [
        {"drug": "Cefiderocol", "class": "Siderophore Cephalosporin", "mechanism": "Siderophore transport bypasses outer membrane; stable against NDM metallo-beta-lactamases."},
        {"drug": "Aztreonam + Avibactam", "class": "Monobactam + BLI", "mechanism": "Aztreonam resists NDM hydrolysis; Avibactam neutralizes co-expressed ESBL/AmpC enzymes."}
    ],
    "blaCTX-M-15": [
        {"drug": "Ertapenem", "class": "Carbapenem", "mechanism": "Stable against CTX-M-15 ESBL hydrolysis."}
    ],
    "armA": [
        {"drug": "Tigecycline", "class": "Glycylcycline", "mechanism": "Binds 30S ribosome at a distinct site, bypassing 16S rRNA G1405 methylation."}
    ],
    "aac(6')-Ib": [
        {"drug": "Plazomicin", "class": "Next-Gen Aminoglycoside", "mechanism": "Engineered side chain blocks acetyltransferase [aac(6')-Ib] binding."}
    ],
    "gyrA_S83L": [
        {"drug": "Trimethoprim-Sulfamethoxazole", "class": "Folate Inhibitor", "mechanism": "Inhibits folic acid biosynthesis; independent of DNA Gyrase mutations."},
        {"drug": "Fosfomycin", "class": "Epoxide", "mechanism": "Inhibits early peptidoglycan MurA enzyme; unaffected by GyrA mutations."}
    ],
    "mcr-1": [
        {"drug": "Tigecycline", "class": "Glycylcycline", "mechanism": "Inhibits protein synthesis; completely bypasses Lipid A membrane modifications."}
    ]
}

PANEL = ["Amoxicillin", "Methicillin", "Imipenem", "Ceftriaxone", "Gentamicin", "Ciprofloxacin", "Vancomycin", "Colistin"]

def generate_susceptibility_report(analysis_result):
    gram_type = analysis_result["gram_type"]
    detected = analysis_result["detected_genes"]
    detected_keys = {g["gene_key"] for g in detected}
    
    primary_failed_drugs = set()
    gene_drug_map = {}
    
    for item in detected:
        for drug in item["affected_drugs"]:
            primary_failed_drugs.add(drug)
            gene_drug_map[drug] = item

    report = []
    for drug in PANEL:
        # Intrinsic Spectrum Constraints
        if drug in SPECTRUM_RULES and gram_type == "Gram-Negative":
            report.append({
                "Antibiotic": drug,
                "Prediction": "Intrinsic Non-Spectrum ⚠️",
                "Confidence": "High (Biological Spectrum)",
                "Primary Determinant": "Gram-Negative Outer Membrane Barrier",
                "Mechanism Explanation": SPECTRUM_RULES[drug]["reason"],
                "Alternatives": [],
                "Alternatives_JSON": "[]"
            })
            continue

        # Genomic Resistance Hits
        if drug in gene_drug_map:
            hit = gene_drug_map[drug]
            gene_key = hit["gene_key"]
            
            explanation = DYNAMIC_MECHANISMS.get((drug, gene_key), hit["mechanism"])
            
            # Cross-resistance filter
            raw_candidates = RESCUE_CANDIDATES.get(gene_key, [])
            valid_alternatives = []
            
            for cand in raw_candidates:
                cand_name = cand["drug"]
                
                if cand_name in primary_failed_drugs:
                    continue
                if cand_name == "Trimethoprim-Sulfamethoxazole" and "dfrA14" in detected_keys:
                    continue
                if cand_name == "Fosfomycin" and "fosA" in detected_keys:
                    continue
                if cand_name == "Plazomicin" and "armA" in detected_keys:
                    continue
                if cand_name == "Colistin" and "mcr-1" in detected_keys:
                    continue
                
                valid_alternatives.append(cand)

            report.append({
                "Antibiotic": drug,
                "Prediction": "Predicted Resistant ❌",
                "Confidence": f"High ({hit['identity_pct']:.1f}% Ident / {hit['coverage_pct']:.1f}% Cov)",
                "Primary Determinant": hit["name"],
                "Mechanism Explanation": explanation,
                "Alternatives": valid_alternatives,
                "Alternatives_JSON": json.dumps(valid_alternatives)  # Valid JSON String for CSV Export
            })
        else:
            report.append({
                "Antibiotic": drug,
                "Prediction": "Predicted Susceptible ✅",
                "Confidence": "Moderate (No Known ARG Detected)",
                "Primary Determinant": "None Detected",
                "Mechanism Explanation": "No cataloged resistance genes identified in assembly. Susceptibility assumed pending lab AST.",
                "Alternatives": [],
                "Alternatives_JSON": "[]"
            })

    return report