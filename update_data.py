#!/usr/bin/env python3
"""
GEO CNS Tumor (中枢神经系统肿瘤) 数据增量更新脚本
"""

import os
import json
import time
import re
import requests
from datetime import datetime, timedelta
from Bio import Entrez

NCBI_EMAIL = os.environ.get('NCBI_EMAIL', '')
NCBI_API_KEY = os.environ.get('NCBI_API_KEY', '')
MINIMAX_API_KEY = os.environ.get('MINIMAX_API_KEY', '')

SEARCH_CONFIG = {
    "keywords": [
        # 胶质瘤
        "glioma", "glioblastoma", "GBM", "astrocytoma",
        "oligodendroglioma", "diffuse glioma", "low grade glioma",
        "high grade glioma", "IDH mutant glioma",
        # 脑转移瘤（通用术语）
        "brain metastasis", "brain metastases",
        "cerebral metastasis", "cerebral metastases",
        "leptomeningeal metastasis", "leptomeningeal metastases",
        "intracranial metastasis", "intracranial metastases",
        "secondary brain tumor", "secondary brain tumour",
        "metastatic brain tumor", "metastatic brain tumour",
        "brain colonization",
        # 脑转移瘤（按常见原发瘤来源）
        "lung cancer brain metastas", "NSCLC brain metastas",
        "breast cancer brain metastas", "melanoma brain metastas",
        "renal cell carcinoma brain", "colorectal cancer brain",
        "HER2 brain metastas",
        # 脑膜瘤
        "meningioma",
        # 垂体瘤
        "pituitary adenoma", "pituitary tumor", "pituitary tumour",
        "pituitary neuroendocrine tumor",
        # 听神经瘤
        "vestibular schwannoma", "acoustic neuroma",
        # 其他中枢神经系统肿瘤
        "medulloblastoma", "ependymoma", "craniopharyngioma",
        "CNS lymphoma", "brain tumor", "brain tumour",
        "central nervous system tumor", "central nervous system tumour",
        "intracranial tumor", "intracranial tumour",
        "DIPG", "diffuse intrinsic pontine glioma",
        "choroid plexus tumor", "neurocytoma",
    ],
    "organisms": ["Homo sapiens", "Mus musculus"],
    "data_types": [
        "Expression profiling by high throughput sequencing",
        "Expression profiling by array",
        "Methylation profiling by array",
        "Methylation profiling by high throughput sequencing",
        "Genome binding/occupancy profiling by high throughput sequencing",
        "Non-coding RNA profiling by high throughput sequencing",
        "Genome variation profiling by high throughput sequencing",
    ]
}

DATA_FILE = "data/geo_data.json"

# CNS-tumor-relevance terms: datasets must mention at least one of these in title/summary
CNS_TERMS = [
    # 胶质瘤
    "glioma", "glioblastoma", "gbm", "astrocytoma", "oligodendroglioma",
    "glial tumor", "glial tumour", "gliogenesis",
    # 脑转移瘤（通用）
    "brain metastas", "cerebral metastas", "intracranial metastas",
    "leptomeningeal metastas", "secondary brain tumor", "secondary brain tumour",
    "metastatic brain", "brain coloniz", "brm ",
    # 脑转移瘤（原发瘤语境下的脑转移）
    "brain metastatic niche", "metastatic niche brain",
    # 脑膜瘤
    "meningioma",
    # 垂体瘤
    "pituitary",
    # 听神经瘤
    "schwannoma", "acoustic neuroma",
    # 其他 CNS 肿瘤
    "medulloblastoma", "ependymoma", "craniopharyngioma",
    "cns lymphoma", "cns tumor", "cns tumour",
    "intracranial tumor", "intracranial tumour", "intracranial neoplasm",
    "brain cancer", "brain tumor", "brain tumour",
    "neuro-oncol", "neurooncol",
    "dipg", "pontine glioma", "choroid plexus", "neurocytoma",
]

# Datasets primarily about these topics (without CNS relevance) are excluded
EXCLUDE_TERMS = ["skin fibrosis", "dermal fibrosis", "cardiac fibrosis"]

# ── 肿瘤类型标准化映射 ──────────────────────────────────────────────
# 按优先级从高到低排列：先匹配更具体的子类型，再匹配宽泛类别
# 每条规则: (标准名称, [匹配关键词列表])
# ── 组学数据类型标准化映射 ─────────────────────────────────────────
# 按优先级从高到低：先匹配更具体的技术，再匹配宽泛类别
DATA_TYPE_RULES = [
    # ── 单细胞/空间组学（优先匹配）──
    ("Spatial Transcriptomics",  ["spatial transcriptom", "visium", "slide-seq", "merfish",
                                  "seqfish", "stereo-seq", "10x visium", "spatial rna",
                                  "st-seq", "hdst", "dbit-seq"]),
    ("scRNA-seq",                ["single-cell rna", "scrna-seq", "scrna seq", "single cell rna",
                                  "10x genomics", "10x chromium", "drop-seq", "dropseq",
                                  "smart-seq", "smartseq", "cel-seq", "celseq",
                                  "single-cell transcriptom", "single cell transcriptom",
                                  "indrops", "sci-rna"]),
    ("scATAC-seq",               ["scatac-seq", "scatac seq", "single-cell atac",
                                  "single cell atac", "snATAC", "snatac"]),
    ("snRNA-seq",                ["snrna-seq", "snrna seq", "single-nucleus rna",
                                  "single nucleus rna", "single-nuclei rna"]),
    ("Multi-omics (single-cell)", ["multiome", "multi-ome", "cite-seq", "citeseq",
                                   "tea-seq", "share-seq", "paired-seq",
                                   "single-cell multi", "single cell multi"]),

    # ── 表观遗传组学 ──
    ("CUT&Tag/CUT&RUN",         ["cut&tag", "cut&run", "cuttag", "cutrun",
                                  "cut and tag", "cut and run", "cleavage under targets"]),
    ("ChIP-seq",                 ["chip-seq", "chipseq", "chip seq",
                                  "chromatin immunoprecipitation sequencing"]),
    ("ATAC-seq",                 ["atac-seq", "atac seq", "atacseq"]),
    ("Methylation Array",        ["methylation array", "methylation profiling by array",
                                  "450k", "850k", "epic array", "infinium",
                                  "methylation bead", "hm450", "hm27"]),
    ("Bisulfite-seq",            ["bisulfite", "wgbs", "rrbs", "methylation profiling by high throughput",
                                  "em-seq", "methylc-seq"]),
    ("Hi-C",                     ["hi-c", "hic ", "hi c ", "chromosome conformation",
                                  "3c-seq", "4c-seq", "capture-c"]),

    # ── 非编码 RNA ──
    ("miRNA/lncRNA",             ["mirna", "microrna", "lncrna", "long non-coding",
                                  "long noncoding", "circrna", "circular rna",
                                  "small rna", "non-coding rna profiling"]),

    # ── 基因组变异 ──
    ("WES/WGS",                  ["whole exome", "whole genome sequencing", "wes ", "wgs ",
                                  "exome sequencing", "genome variation profiling",
                                  "targeted sequencing", "panel sequencing"]),

    # ── 蛋白组学 ──
    ("Proteomics",               ["proteom", "mass spectrometry", "tmt labeling",
                                  "itraq", "silac", "phosphoproteom"]),

    # ── 转录组（兜底）──
    ("Microarray",               ["expression profiling by array", "microarray",
                                  "affymetrix", "agilent", "illumina beadchip",
                                  "gene expression array"]),
    ("bulk RNA-seq",             ["rna-seq", "rna seq", "rnaseq", "mrna-seq",
                                  "transcriptome sequencing",
                                  "expression profiling by high throughput sequencing"]),
]


def classify_data_type(title, summary, overall_design, geo_data_type=""):
    """根据标题、摘要、实验设计和GEO原始数据类型，识别组学数据类型。"""
    combined = (title + " " + summary + " " + overall_design + " " + geo_data_type).lower()
    for std_name, keywords in DATA_TYPE_RULES:
        if any(kw in combined for kw in keywords):
            return std_name
    return "Other"


TUMOR_TYPE_RULES = [
    # ── 胶质瘤子类型（先匹配具体型再匹配宽泛型）──
    ("Glioma - GBM",          ["glioblastoma", "gbm "]),
    ("Glioma - Astrocytoma",  ["astrocytoma"]),
    ("Glioma - Oligodendroglioma", ["oligodendroglioma"]),
    ("Glioma - DIPG",         ["dipg", "diffuse intrinsic pontine glioma"]),
    ("Glioma",                ["glioma", "glial tumor", "glial tumour"]),

    # ── 脑转移瘤（先按原发瘤分，再归通用类）──
    ("Brain Metastasis - Lung",       ["lung cancer brain", "nsclc brain", "lung adenocarcinoma brain",
                                       "small cell lung cancer brain", "sclc brain"]),
    ("Brain Metastasis - Breast",     ["breast cancer brain", "her2 brain", "triple negative brain"]),
    ("Brain Metastasis - Melanoma",   ["melanoma brain metastas", "melanoma brain coloniz",
                                       "melanoma cerebral", "melanoma leptomeningeal"]),
    ("Brain Metastasis - Renal",      ["renal cell carcinoma brain", "renal brain metastas",
                                       "kidney cancer brain"]),
    ("Brain Metastasis - Colorectal", ["colorectal cancer brain", "colon cancer brain"]),
    ("Brain Metastasis",              ["brain metastas", "cerebral metastas", "intracranial metastas",
                                       "leptomeningeal metastas", "secondary brain tumor",
                                       "secondary brain tumour", "metastatic brain", "brain coloniz",
                                       "brm "]),

    # ── 脑膜瘤 ──
    ("Meningioma",            ["meningioma"]),

    # ── 垂体瘤 ──
    ("Pituitary Tumor",       ["pituitary adenoma", "pituitary tumor", "pituitary tumour",
                               "pituitary neuroendocrine", "pituitary carcinoma",
                               "pituitary blastoma"]),

    # ── 听神经瘤 ──
    ("Vestibular Schwannoma", ["vestibular schwannoma", "acoustic neuroma", "acoustic schwannoma"]),

    # ── 其他 CNS 肿瘤 ──
    ("Medulloblastoma",       ["medulloblastoma"]),
    ("Ependymoma",            ["ependymoma"]),
    ("Craniopharyngioma",     ["craniopharyngioma"]),
    ("CNS Lymphoma",          ["cns lymphoma", "primary central nervous system lymphoma", "pcnsl"]),
    ("Choroid Plexus Tumor",  ["choroid plexus"]),
    ("Neurocytoma",           ["neurocytoma"]),

    # ── 兜底 ──
    ("Other CNS Tumor",       ["brain tumor", "brain tumour", "brain cancer",
                               "central nervous system tumor", "central nervous system tumour",
                               "intracranial tumor", "intracranial tumour",
                               "intracranial neoplasm"]),
]


def classify_tumor_type(title, summary):
    """根据标题和摘要将数据集归类到标准肿瘤类型。
    返回 (主类型, 子类型) 元组，例如 ("Glioma", "Glioma - GBM")。
    """
    combined = (title + " " + summary).lower()
    for std_name, keywords in TUMOR_TYPE_RULES:
        if any(kw in combined for kw in keywords):
            # 主类型取 " - " 之前的部分
            main_type = std_name.split(" - ")[0] if " - " in std_name else std_name
            return main_type, std_name
    return "Other CNS Tumor", "Other CNS Tumor"


def is_cns_relevant(record):
    """Check if a dataset is relevant to CNS tumors."""
    title = record.get("title", "").lower()
    summary = record.get("summary", "").lower()
    combined = title + " " + summary
    has_cns = any(term in combined for term in CNS_TERMS)
    if has_cns:
        return True
    # Exclude datasets that are clearly about non-CNS topics
    has_exclude = any(term in combined for term in EXCLUDE_TERMS)
    if has_exclude:
        return False
    # No CNS terms and no exclude terms — still reject (be strict for relevance)
    return False


def setup_entrez():
    Entrez.email = NCBI_EMAIL
    if NCBI_API_KEY:
        Entrez.api_key = NCBI_API_KEY


def build_query():
    keyword_query = " OR ".join([f'"{kw}"' for kw in SEARCH_CONFIG["keywords"]])
    org_query = " OR ".join([f'"{org}"[Organism]' for org in SEARCH_CONFIG["organisms"]])
    type_query = " OR ".join([f'"{t}"[DataSet Type]' for t in SEARCH_CONFIG["data_types"]])
    return f"({keyword_query}) AND ({org_query}) AND ({type_query})"


def search_geo(max_retries=3, is_initial=False):
    """搜索 GEO 数据库（带重试机制）
    is_initial=True 时不限制日期范围，用于首次全量抓取。
    """
    query = build_query()
    end_date = datetime.now().strftime("%Y/%m/%d")

    search_params = {
        "db": "gds", "term": query, "retmax": 10000, "usehistory": "y",
    }

    if is_initial:
        print("首次运行：不限制日期范围，抓取全部历史数据")
    else:
        start_date = (datetime.now() - timedelta(days=30)).strftime("%Y/%m/%d")
        search_params["mindate"] = start_date
        search_params["maxdate"] = end_date
        search_params["datetype"] = "pdat"
        print(f"日期范围: {start_date} - {end_date}")

    print(f"搜索查询: {query[:100]}...")

    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(**search_params)
            results = Entrez.read(handle)
            handle.close()
            return results.get("IdList", [])
        except Exception as e:
            print(f"搜索失败 (尝试 {attempt + 1}/{max_retries}): {e}")
            if attempt < max_retries - 1:
                time.sleep(10)
            else:
                print("所有重试都失败了")
                return []


def fetch_summaries(id_list, max_retries=3):
    """获取数据集摘要（带重试机制）"""
    if not id_list:
        return []

    for attempt in range(max_retries):
        try:
            handle = Entrez.esummary(db="gds", id=",".join(id_list))
            records = Entrez.read(handle)
            handle.close()
            return records
        except Exception as e:
            print(f"获取摘要失败 (尝试 {attempt + 1}/{max_retries}): {e}")
            if attempt < max_retries - 1:
                time.sleep(10)
            else:
                return []


def clean_pubmed_ids(pubmed_str):
    if not pubmed_str:
        return ""
    numbers = re.findall(r'IntegerElement\((\d+)', str(pubmed_str))
    if numbers:
        return "; ".join(numbers)
    numbers = re.findall(r'\d+', str(pubmed_str))
    if numbers:
        return "; ".join(numbers)
    return str(pubmed_str)


def generate_ai_summary(title, summary, data_type):
    if not MINIMAX_API_KEY:
        return ""

    prompt = f"""请用中文为以下GEO数据集生成一个精炼的科研摘要（80-120字）：

标题: {title}
数据类型: {data_type}
研究摘要: {summary[:800]}

请直接输出中文摘要："""

    try:
        response = requests.post(
            'https://api.minimaxi.com/v1/chat/completions',
            headers={
                "Authorization": f'Bearer {MINIMAX_API_KEY}',
                "Content-Type": "application/json"
            },
            json={
                "model": "MiniMax-M2.1",
                "messages": [{"role": "user", "content": prompt}],
                "max_tokens": 1500,
                "temperature": 0.7
            },
            timeout=60
        )
        if response.status_code == 200:
            content = response.json()["choices"][0]["message"]["content"]
            return re.sub(r'<think>.*?</think>', '', content, flags=re.DOTALL).strip()
    except Exception as e:
        print(f"AI 摘要生成失败: {e}")
    return ""


def fetch_geo_soft(accession):
    """获取GEO SOFT格式的详细信息（Country, Lab, Institute等）"""
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=full"
    try:
        response = requests.get(url, timeout=30)
        if response.status_code != 200:
            return {}

        info = {
            "overall_design": "",
            "contributors": [],
            "lab": "",
            "institute": "",
            "country": "",
            "series_type": [],
        }

        for line in response.text.split('\n'):
            line = line.strip()
            if line.startswith('!Series_overall_design'):
                info["overall_design"] = line.split('=', 1)[1].strip()
            elif line.startswith('!Series_type'):
                info["series_type"].append(line.split('=', 1)[1].strip())
            elif line.startswith('!Series_contributor'):
                contributor = line.split('=', 1)[1].strip()
                parts = contributor.split(',')
                if len(parts) >= 2:
                    name = f"{parts[-1]} {parts[0]}".strip()
                    if name and name not in info["contributors"]:
                        info["contributors"].append(name)
            elif line.startswith('!Series_contact_laboratory'):
                info["lab"] = line.split('=', 1)[1].strip()
            elif line.startswith('!Series_contact_institute'):
                info["institute"] = line.split('=', 1)[1].strip()
            elif line.startswith('!Series_contact_country'):
                info["country"] = line.split('=', 1)[1].strip()

        return info
    except Exception as e:
        print(f"    获取SOFT信息失败: {e}")
        return {}


def parse_record(record):
    accession = record.get("Accession", "")
    if not accession.startswith("GSE"):
        return None

    pubmed_ids = record.get("PubMedIds", [])
    pubmed_str = clean_pubmed_ids("; ".join(str(p) for p in pubmed_ids) if pubmed_ids else "")

    title = record.get("title", "")
    summary = record.get("summary", "")

    # 获取详细SOFT信息
    soft_info = fetch_geo_soft(accession)
    time.sleep(0.3)

    # 智能识别数据类型
    geo_data_type = "; ".join(soft_info.get("series_type", []))
    overall_design = soft_info.get("overall_design", "")
    data_type = classify_data_type(title, summary, overall_design, geo_data_type)

    # 标准化肿瘤类型
    tumor_main, tumor_sub = classify_tumor_type(title, summary)

    ai_summary = generate_ai_summary(title, summary, data_type)
    if ai_summary:
        time.sleep(1)

    return {
        "Accession": accession,
        "Title": title,
        "Organism": record.get("taxon", ""),
        "Tumor_Type": tumor_main,
        "Tumor_Subtype": tumor_sub,
        "Data_Type": data_type,
        "Sample_Count": record.get("n_samples", 0),
        "Platform": record.get("GPL", ""),
        "Country": soft_info.get("country", ""),
        "Lab": soft_info.get("lab", ""),
        "Institute": soft_info.get("institute", ""),
        "Contributors": "; ".join(soft_info.get("contributors", [])),
        "PubMed_IDs": pubmed_str,
        "Supplementary_Size": "N/A",
        "Summary": summary,
        "Overall_Design": soft_info.get("overall_design", ""),
        "AI_Summary_CN": ai_summary,
        "AI_Summary": ai_summary,
        "GEO_Link": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
        "Submission_Date": record.get("PDAT", ""),
    }


def main():
    print(f"开始更新 CNS Tumor 数据 - {datetime.now()}")

    if not NCBI_EMAIL:
        print("错误: 未设置 NCBI_EMAIL")
        return

    setup_entrez()

    if os.path.exists(DATA_FILE):
        with open(DATA_FILE, 'r', encoding='utf-8') as f:
            existing_data = json.load(f)
    else:
        existing_data = []

    existing_accessions = {d["Accession"] for d in existing_data}
    print(f"现有数据集: {len(existing_data)}")

    # 数据为空时进行全量抓取，否则增量更新（最近30天）
    is_initial = len(existing_data) == 0
    id_list = search_geo(is_initial=is_initial)
    print(f"搜索到: {len(id_list)} 条记录")

    if not id_list:
        print("没有新数据")
        return

    summaries = fetch_summaries(id_list)

    new_count = 0
    skipped_irrelevant = 0
    for record in summaries:
        accession = record.get("Accession", "")
        if accession in existing_accessions or not accession.startswith("GSE"):
            continue

        if not is_cns_relevant(record):
            skipped_irrelevant += 1
            print(f"  跳过 (非CNS肿瘤相关): {accession}")
            continue

        parsed = parse_record(record)
        if parsed:
            existing_data.insert(0, parsed)
            existing_accessions.add(accession)
            new_count += 1
            print(f"  新增: {accession}")

    if skipped_irrelevant:
        print(f"跳过 {skipped_irrelevant} 条非CNS肿瘤相关数据集")

    if new_count > 0:
        with open(DATA_FILE, 'w', encoding='utf-8') as f:
            json.dump(existing_data, f, ensure_ascii=False, indent=2)
        print(f"完成! 新增 {new_count} 条，总计 {len(existing_data)} 条")
    else:
        print("没有新数据需要添加")


if __name__ == "__main__":
    main()
