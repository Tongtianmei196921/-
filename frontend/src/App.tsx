import { useState, useCallback, useEffect, useRef } from 'react'
import { useDropzone } from 'react-dropzone'
import { toPng, toSvg } from 'html-to-image'
import { jsPDF } from 'jspdf'
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
  ResponsiveContainer, Cell, ScatterChart, Scatter, ZAxis,
} from 'recharts'

// ─── Types ───────────────────────────────────────────────────────────────────

type Page = 'tool' | 'tutorial'
type Step = 'upload' | 'config' | 'running' | 'results'
type InputMode = 'file' | 'geo'
type Locale = 'zh' | 'en'

interface CheckpointStatus {
  all_ready: boolean
  checkpoints: Record<string, boolean>
  checkpoint_dir?: string
}

interface CompoundAnnotation {
  display_name?: string | null
  full_broad_id?: string | null
  clinical_phase?: string | null
  moa?: string | null
  target?: string | null
  disease_area?: string | null
  indication?: string | null
  vendor_name?: string | null
  smiles?: string | null
  inchikey?: string | null
  trade_name?: string | null
  pubchem_title?: string | null
  chemical_name?: string | null
  molecular_formula?: string | null
  molecular_weight?: number | null
  canonical_smiles?: string | null
  isomeric_smiles?: string | null
  synonyms?: string[] | null
  pubchem_cid?: number | null
  pubchem_url?: string | null
  structure_image?: string | null
  source?: string | null
}

interface DirectionEvidence {
  label?: string | null
  score?: number | null
  source?: string | null
  reason?: string | null
  query_up_genes?: number | null
  query_down_genes?: number | null
}

interface CompoundResult {
  compound: string
  rank: number
  prob: number
  logit: number
  annotation?: CompoundAnnotation | null
  direction?: DirectionEvidence | null
}

interface PredictionResponse {
  success: boolean
  data: { samples: string[]; results: Record<string, CompoundResult[]> }
  meta: { n_samples: number; n_compounds: number; n_genes: number; filename: string; [key: string]: unknown }
}

interface PreparationSummary {
  filename: string
  status: 'ready' | 'needs_configuration'
  mode: string
  original_shape: number[]
  prepared_shape: number[] | null
  sample_names: string[]
  top_genes: Array<{ gene: string; score: number }>
  notes: string[]
  candidate_groups: Array<{ column: string; values: string[]; n_unique: number; score: number }>
  candidate_sample_ids?: Array<{ column: string; n_unique: number; score: number }>
  detected_data_type?: string | null
  used_sample_id_column?: string | null
}

interface PrepOptions {
  group_column?: string
  group1_value?: string
  group2_value?: string
  sample_id_column?: string
}

type SortKey = 'rank' | 'prob' | 'logit' | 'compound'

interface GeoPreviewResponse {
  accession: string
  n_samples: number
  n_probe_rows: number
  n_genes: number
  platform_id?: string | null
  organism?: string | null
  symbol_source: string
  expression_source?: string
  ortholog_mapping?: {
    source: string
    input_genes: number
    mapped_input_genes: number
    unique_human_symbols: number
    orthology_types: string[]
  } | null
  used_log2: boolean
  detected_grouping?: {
    group_column: string
    group1_value: string
    group2_value: string
    group1_count: number
    group2_count: number
    mode: string
  } | null
  candidate_columns: Array<{ column: string; n_unique: number; values: string[]; score: number }>
}

interface GeoRunConfig {
  accession: string
  control_keyword?: string
  case_keyword?: string
  group_column?: string
  group1_value?: string
  group2_value?: string
}

// ─── API ─────────────────────────────────────────────────────────────────────

const API_BASE_URL = (import.meta.env.VITE_API_BASE_URL || '').replace(/\/+$/, '')

function apiUrl(path: string): string {
  return `${API_BASE_URL}${path}`
}

function normalizeUserError(message: unknown): string {
  const text = String(message || '').trim()
  if (!text) return '未知错误'
  if (text.includes('unknown url type') && text.includes('NONE')) {
    return '该 GEO 记录包含 NONE 这类空的补充文件链接，系统已判定它不是可下载表达矩阵。请更换 GEO 编号，或下载原始数据整理成 h5ad / 表达矩阵后上传。'
  }
  if (
    text.includes('does not expose an expression matrix') ||
    text.includes('no per-sample supplementary expression files')
  ) {
    return '该 GEO 暂时不能自动解析：series matrix 中没有可用表达矩阵，也没有逐样本表达补充文件。请换一个包含表达矩阵的 GEO，或把数据整理成 h5ad / CSV 后上传。'
  }
  return text
}

async function fetchCheckpoints(): Promise<CheckpointStatus> {
  const res = await fetch(apiUrl('/api/checkpoints'))
  if (!res.ok) throw new Error(`HTTP ${res.status}`)
  return res.json()
}

function generateSampleCSV(): string {
  // LCG pseudo-random, seed=42 — deterministic, no backend needed
  let s = 42
  const rand = () => {
    s = Math.imul(s, 1664525) + 1013904223 | 0
    return (s >>> 0) / 0xffffffff * 4 - 2   // range [-2, 2]
  }
  const genes = Array.from({ length: 978 }, (_, i) => `GENE_${String(i).padStart(4, '0')}`)
  const samples = ['Sample_A', 'Sample_B', 'Sample_C']
  const rows = [
    ['sample_id', ...genes].join(','),
    ...samples.map(name => [name, ...genes.map(() => rand().toFixed(4))].join(',')),
  ]
  return rows.join('\n')
}

function downloadSampleCSV() {
  const csv = generateSampleCSV()
  const blob = new Blob([csv], { type: 'text/csv' })
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url; a.download = 'sample_vscores.csv'; a.click()
  URL.revokeObjectURL(url)
}

async function runPrediction(file: File, nTop: number, prep?: PrepOptions): Promise<PredictionResponse> {
  const form = new FormData()
  form.append('file', file)
  form.append('n_top', String(nTop))
  if (prep?.group_column) form.append('group_column', prep.group_column)
  if (prep?.group1_value) form.append('group1_value', prep.group1_value)
  if (prep?.group2_value) form.append('group2_value', prep.group2_value)
  if (prep?.sample_id_column) form.append('sample_id_column', prep.sample_id_column)
  const res = await fetch(apiUrl('/api/predict'), { method: 'POST', body: form })
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: '未知错误' }))
    throw new Error(err.detail || `HTTP ${res.status}`)
  }
  return res.json()
}

async function previewGeo(accession: string): Promise<GeoPreviewResponse> {
  const res = await fetch(apiUrl(`/api/geo/preview?accession=${encodeURIComponent(accession)}`))
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: 'Unknown GEO error' }))
    throw new Error(normalizeUserError(err.detail || `HTTP ${res.status}`))
  }
  return res.json()
}

async function runGeoPrediction(config: GeoRunConfig, nTop: number): Promise<PredictionResponse> {
  const res = await fetch(apiUrl('/api/geo/predict'), {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ ...config, n_top: nTop }),
  })
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: 'Unknown GEO error' }))
    throw new Error(normalizeUserError(err.detail || `HTTP ${res.status}`))
  }
  return res.json()
}

async function prepareFileInput(file: File, prep?: PrepOptions): Promise<PreparationSummary> {
  const form = new FormData()
  form.append('file', file)
  if (prep?.group_column) form.append('group_column', prep.group_column)
  if (prep?.group1_value) form.append('group1_value', prep.group1_value)
  if (prep?.group2_value) form.append('group2_value', prep.group2_value)
  if (prep?.sample_id_column) form.append('sample_id_column', prep.sample_id_column)
  const res = await fetch(apiUrl('/api/prepare-input'), { method: 'POST', body: form })
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: 'Unknown preparation error' }))
    throw new Error(err.detail || `HTTP ${res.status}`)
  }
  const data = await res.json()
  return data.preparation
}

async function downloadWordReport(response: PredictionResponse, sample: string, locale: Locale) {
  const res = await fetch(apiUrl('/api/report/docx'), {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ response, sample, locale }),
  })
  if (!res.ok) {
    const err = await res.json().catch(() => ({ detail: 'Word report generation failed' }))
    throw new Error(err.detail || `HTTP ${res.status}`)
  }
  const blob = await res.blob()
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url
  a.download = `drugreflector_${sample.replace(/[^\w.-]+/g, '_')}_report.docx`
  a.click()
  URL.revokeObjectURL(url)
}

function downloadCSV(data: CompoundResult[], filename: string) {
  const header = [
    'compound',
    'display_name',
    'vendor_name',
    'trade_name',
    'chemical_name',
    'pubchem_title',
    'molecular_formula',
    'molecular_weight',
    'direction_label',
    'direction_score',
    'direction_source',
    'direction_reason',
    'rank',
    'prob',
    'logit',
    'clinical_phase',
    'moa',
    'target',
    'pubchem_cid',
    'smiles',
  ].join(',')
  const rows = data.map(r =>
      [
        r.compound,
        csvEscape(getCompoundDisplayName(r)),
        csvEscape(r.annotation?.vendor_name),
        csvEscape(r.annotation?.trade_name),
        csvEscape(r.annotation?.chemical_name),
        csvEscape(r.annotation?.pubchem_title),
        csvEscape(r.annotation?.molecular_formula),
        r.annotation?.molecular_weight ?? '',
        csvEscape(r.direction?.label),
        r.direction?.score ?? '',
        csvEscape(r.direction?.source),
        csvEscape(r.direction?.reason),
        r.rank,
        r.prob.toFixed(6),
        r.logit.toFixed(6),
      csvEscape(r.annotation?.clinical_phase),
      csvEscape(r.annotation?.moa),
      csvEscape(r.annotation?.target),
      r.annotation?.pubchem_cid ?? '',
      csvEscape(r.annotation?.smiles),
    ].join(',')
  )
  const blob = new Blob([[header, ...rows].join('\n')], { type: 'text/csv' })
  const url = URL.createObjectURL(blob)
  const a = document.createElement('a')
  a.href = url; a.download = filename; a.click()
  URL.revokeObjectURL(url)
}

function csvEscape(value: unknown) {
  if (value === null || value === undefined) return ''
  const text = String(value)
  if (!/[",\n]/.test(text)) return text
  return `"${text.replace(/"/g, '""')}"`
}

function getCompoundDisplayName(result: CompoundResult) {
  return result.annotation?.display_name || result.compound
}

function getCompoundShortLabel(result: CompoundResult, maxLength = 28) {
  const name = getCompoundDisplayName(result)
  return name.length > maxLength ? `${name.slice(0, maxLength - 1)}…` : name
}

function getCompoundSecondaryLabel(result: CompoundResult) {
  return result.annotation?.display_name ? result.compound : null
}

function getDirectionClass(direction?: DirectionEvidence | null) {
  const label = direction?.label?.toLowerCase()
  if (label === 'reverse') return 'reverse'
  if (label === 'mimic') return 'mimic'
  return 'no-evidence'
}

function getDirectionDisplayLabel(direction: DirectionEvidence | null | undefined, locale: Locale) {
  const label = direction?.label?.trim().toLowerCase()
  if (label === 'reverse') return locale === 'zh' ? '更像逆转' : 'Reverse'
  if (label === 'mimic') return locale === 'zh' ? '更像增强' : 'Mimic'
  return locale === 'zh' ? '暂无客观证据' : 'No objective evidence'
}

function getDirectionSummary(direction: DirectionEvidence | null | undefined, locale: Locale) {
  const label = direction?.label?.trim().toLowerCase()
  if (label === 'reverse') {
    return locale === 'zh'
      ? '提示该药物更像逆转当前输入 signature，但这仍是签名层面的数据证据，不等于已证实上调或下调某个单基因。'
      : 'Suggests this compound may reverse the input signature, but this remains signature-level evidence rather than proof of up- or down-regulation of any single gene.'
  }
  if (label === 'mimic') {
    return locale === 'zh'
      ? '提示该药物更像增强或模拟当前输入 signature，但这仍是签名层面的数据证据，不等于已证实上调或下调某个单基因。'
      : 'Suggests this compound may mimic the input signature, but this remains signature-level evidence rather than proof of up- or down-regulation of any single gene.'
  }
  return locale === 'zh'
    ? '当前部署缺少可客观判定方向性的外部签名连通性证据，因此这里不能负责任地判断“更像逆转”还是“更像增强”。'
    : 'This deployment lacks objective external signed-connectivity evidence, so it cannot responsibly classify the compound as reverse or mimic.'
}

function getAnnotationCompleteness(result: CompoundResult) {
  let score = 0
  if (result.annotation?.display_name) score += 1
  if (result.annotation?.target || result.annotation?.moa) score += 1
  if (result.annotation?.pubchem_cid || result.annotation?.structure_image) score += 1

  if (score >= 3) return { label: 'High', className: 'high' }
  if (score === 2) return { label: 'Medium', className: 'medium' }
  return { label: 'Low', className: 'low' }
}

function summarizeTokens(values: Array<string | null | undefined>, limit = 8) {
  const counts = new Map<string, number>()
  values.forEach(value => {
    if (!value) return
    value
      .split('|')
      .map(token => token.trim())
      .filter(Boolean)
      .forEach(token => {
        counts.set(token, (counts.get(token) ?? 0) + 1)
      })
  })

  return [...counts.entries()]
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
    .slice(0, limit)
}

function formatRatio(value?: number | null) {
  if (value === null || value === undefined || Number.isNaN(value)) return 'Unavailable'
  return `${(value * 100).toFixed(1)}%`
}

function getDisplaySampleLabel(sample: string, locale: Locale) {
  const arrowParts = sample.split('->').map(part => part.trim())
  if (arrowParts.length !== 2) return sample

  const left = arrowParts[0]
  const right = arrowParts[1]
  const leftParts = left.split(':').map(part => part.trim()).filter(Boolean)
  const rightParts = right.split(':').map(part => part.trim()).filter(Boolean)

  const accession = leftParts[0]?.match(/^GSE\d+$/i) ? leftParts[0] : null
  const controlLabel = leftParts[leftParts.length - 1] || left
  const caseLabel = rightParts[rightParts.length - 1] || right
  const groupLabel = leftParts.length >= 2 ? leftParts[leftParts.length - 2] : null

  if (locale === 'zh') {
    const prefix = accession ? `${accession} | ` : ''
    const groupText = groupLabel ? `${groupLabel} 分组：` : ''
    return `${prefix}${groupText}${caseLabel} vs ${controlLabel}`
  }

  const prefix = accession ? `${accession} | ` : ''
  const groupText = groupLabel ? `${groupLabel}: ` : ''
  return `${prefix}${groupText}${caseLabel} vs ${controlLabel}`
}

function getOrganismLabel(value: string | null | undefined, locale: Locale) {
  if (!value) return locale === 'zh' ? '未知' : 'unknown'
  const normalized = value.toLowerCase()
  if (normalized === 'mus musculus') return locale === 'zh' ? '小鼠' : 'Mus musculus'
  if (normalized === 'homo sapiens') return locale === 'zh' ? '人' : 'Homo sapiens'
  return value
}

function getExpressionSourceLabel(value: string | null | undefined, locale: Locale) {
  if (!value) return locale === 'zh' ? '未知' : 'unknown'
  if (locale === 'zh' && value === 'supplementary:supplementary_file_2') return 'GEO 补充文件（supplementary_file_2）'
  return value
}

function getGeneMappingLabel(value: string, locale: Locale) {
  if (locale !== 'zh') return value
  if (value === 'mouse_to_human_orthologs') return '小鼠基因已客观映射为人同源基因'
  if (value === 'ensembl_ids') return 'Ensembl 基因 ID'
  if (value === 'gene_symbols') return '人类基因符号'
  return value
}

async function copyText(text: string) {
  try {
    await navigator.clipboard.writeText(text)
  } catch {
    const input = document.createElement('textarea')
    input.value = text
    document.body.appendChild(input)
    input.select()
    document.execCommand('copy')
    input.remove()
  }
}

function fmtBytes(n: number) {
  if (n < 1024) return `${n} B`
  if (n < 1024 ** 2) return `${(n / 1024).toFixed(1)} KB`
  return `${(n / 1024 ** 2).toFixed(1)} MB`
}

// ─── Tooltip helper ──────────────────────────────────────────────────────────

function Tip({ text }: { text: string }) {
  return (
    <span className="tip">
      <em className="tip-icon">?</em>
      <span className="tip-box">{text}</span>
    </span>
  )
}

// ─── Step indicator ───────────────────────────────────────────────────────────

const STEP_LABELS: Record<Locale, { id: Step; label: string }[]> = {
  zh: [
    { id: 'upload', label: '上传数据' },
    { id: 'config', label: '配置参数' },
    { id: 'running', label: '运行中' },
    { id: 'results', label: '预测结果' },
  ],
  en: [
    { id: 'upload', label: 'Upload' },
    { id: 'config', label: 'Configure' },
    { id: 'running', label: 'Running' },
    { id: 'results', label: 'Results' },
  ],
}

function StepIndicator({ current, locale }: { current: Step; locale: Locale }) {
  const steps = STEP_LABELS[locale]
  const idx = steps.findIndex(s => s.id === current)
  return (
    <div className="steps">
      {steps.map((s, i) => {
        const state = i < idx ? 'done' : i === idx ? 'active' : 'idle'
        return (
          <div key={s.id} className="step-item">
            {i > 0 && <div className={`step-connector ${i <= idx ? 'done' : ''}`} />}
            <div className={`step-circle ${state}`}>{i < idx ? '✓' : i + 1}</div>
            <span className={`step-label ${state}`}>{s.label}</span>
          </div>
        )
      })}
    </div>
  )
}

// ─── Tutorial page ────────────────────────────────────────────────────────────

function TutorialPage({ onStart }: { onStart: () => void }) {
  return (
    <>
      {/* Hero */}
      <div className="tut-hero">
        <div className="tut-hero-tag">使用教程</div>
        <div className="tut-hero-title">5 分钟上手 DrugReflector</div>
        <div className="tut-hero-desc">
          DrugReflector 是一个<strong>虚拟药物筛选</strong>工具，
          通过分析基因表达特征，预测哪些化合物最有可能逆转您的疾病或细胞状态——
          无需编程基础。
        </div>
        <button className="btn btn-primary" onClick={onStart}>立即开始使用 →</button>
      </div>

      {/* Workflow */}
      <div className="card tut-section">
        <div className="tut-section-title">🗺 整体工作流程</div>
        <div className="workflow">
          <div className="wf-step">
            <div className="wf-icon c1">🧬</div>
            <div className="wf-title">① 准备数据</div>
            <div className="wf-desc">计算基因表达差异（v-score）或直接使用 CSV</div>
          </div>
          <div className="wf-arrow">→</div>
          <div className="wf-step">
            <div className="wf-icon c2">☁️</div>
            <div className="wf-title">② 上传文件</div>
            <div className="wf-desc">支持 H5AD / CSV / TSV 格式</div>
          </div>
          <div className="wf-arrow">→</div>
          <div className="wf-step">
            <div className="wf-icon c3">⚡</div>
            <div className="wf-title">③ 运行预测</div>
            <div className="wf-desc">三个神经网络集成模型，通常 10–60 秒</div>
          </div>
          <div className="wf-arrow">→</div>
          <div className="wf-step">
            <div className="wf-icon c4">📊</div>
            <div className="wf-title">④ 解读结果</div>
            <div className="wf-desc">查看候选化合物排名并下载 CSV</div>
          </div>
        </div>
      </div>

      {/* Use cases */}
      <div className="card tut-section">
        <div className="tut-section-title">🎯 典型应用场景</div>
        <div className="use-cases">
          <div className="use-case">
            <span className="use-case-icon">🔬</span>
            <div>
              <div className="use-case-title">疾病 vs 正常组织</div>
              <div className="use-case-desc">寻找能将疾病细胞表达谱转变回正常状态的化合物</div>
            </div>
          </div>
          <div className="use-case">
            <span className="use-case-icon">🔄</span>
            <div>
              <div className="use-case-title">细胞状态转变</div>
              <div className="use-case-desc">例如造血干细胞 → 淋巴细胞的分化诱导</div>
            </div>
          </div>
          <div className="use-case">
            <span className="use-case-icon">💊</span>
            <div>
              <div className="use-case-title">药物重定向</div>
              <div className="use-case-desc">为已知靶点寻找新适应症或协同用药</div>
            </div>
          </div>
        </div>
      </div>

      {/* Data format */}
      <div className="card tut-section">
        <div className="tut-section-title">📁 数据格式要求</div>
        <p style={{ fontSize: '.88rem', color: 'var(--mid)', marginBottom: 14, lineHeight: 1.65 }}>
          支持两种输入格式：<strong>H5AD</strong>（AnnData 对象，包含 v-score）或
          <strong> CSV/TSV</strong>（表格文件）。
        </p>
        <table className="fmt-table">
          <thead>
            <tr>
              <th>格式</th><th>行（Row）</th><th>列（Column）</th><th>说明</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td><code>.h5ad</code></td>
              <td>样本 / 细胞转变</td>
              <td>基因（978 个 landmark 基因）</td>
              <td>推荐，保留更多元数据</td>
            </tr>
            <tr>
              <td><code>.csv / .tsv</code></td>
              <td>样本名（第一列）</td>
              <td>基因名（HGNC 格式）</td>
              <td>最简单，用 Excel 也能制作</td>
            </tr>
          </tbody>
        </table>

        <div className="warn-box" style={{ marginTop: 14 }}>
          ⚠️ <strong>关键：</strong>基因名必须使用 <strong>HGNC 格式</strong>（如 <code>TP53</code>、<code>EGFR</code>、<code>CDKN1A</code>）。
          DrugReflector 会自动将小写、带后缀的基因名转换为标准格式。
        </div>

        <div style={{ marginTop: 16 }}>
          <div style={{ fontSize: '.84rem', fontWeight: 600, marginBottom: 8 }}>✅ CSV 示例（前 3 行）</div>
          <div className="code-block">
            <pre>{`sample_id,TP53,EGFR,CDKN1A,MYC,...
Sample_A,1.23,-0.87,0.45,2.10,...
Sample_B,-0.56,1.34,-1.20,0.78,...`}</pre>
          </div>
        </div>

        <div style={{ marginTop: 4 }}>
          <button
            onClick={downloadSampleCSV}
            style={{ background: 'none', border: 'none', cursor: 'pointer', color: 'var(--blue)', fontWeight: 500, fontFamily: 'inherit', fontSize: '.84rem', padding: 0 }}
          >
            ⬇ 下载示例 CSV 文件
          </button>
          <span style={{ fontSize: '.8rem', color: 'var(--mid)', marginLeft: 8 }}>（3 个样本 × 978 个基因，随机生成）</span>
        </div>
      </div>

      {/* V-score */}
      <div className="card tut-section">
        <div className="tut-section-title">🧮 如何计算 v-score？</div>
        <p style={{ fontSize: '.88rem', color: 'var(--mid)', marginBottom: 12, lineHeight: 1.65 }}>
          v-score 是一种衡量两组细胞之间基因表达差异的指标，类似于 fold change，
          但对异质性更鲁棒。如果您已有 v-score 数据，可以直接上传。
          如果没有，可以用以下 Python 代码计算：
        </p>
        <div className="code-block">
          <pre>
<span className="kw">import</span> drugreflector <span className="kw">as</span> dr{'\n'}
<span className="kw">import</span> scanpy <span className="kw">as</span> sc{'\n\n'}
<span className="cm"># 加载您的单细胞数据（AnnData 格式）</span>{'\n'}
<span className="va">adata</span> = sc.read_h5ad(<span className="st">"your_data.h5ad"</span>){'\n\n'}
<span className="cm"># 计算两组细胞之间的 v-score</span>{'\n'}
<span className="va">vscores</span> = dr.<span className="fn">compute_vscores_adata</span>({'\n'}
{'    '}<span className="va">adata</span>,{'\n'}
{'    '}group_col=<span className="st">"cell_type"</span>,       <span className="cm"># 分组列名</span>{'\n'}
{'    '}group1_value=<span className="st">"正常细胞"</span>,    <span className="cm"># 对照组</span>{'\n'}
{'    '}group2_value=<span className="st">"疾病细胞"</span>,    <span className="cm"># 实验组</span>{'\n'}
){'\n\n'}
<span className="cm"># vscores 是一个 pandas Series，可以直接传给模型</span>{'\n'}
<span className="va">model</span> = dr.<span className="fn">DrugReflector</span>(checkpoint_paths=[...]){'\n'}
<span className="va">predictions</span> = <span className="va">model</span>.<span className="fn">predict</span>(<span className="va">vscores</span>, n_top=<span className="st">50</span>)
          </pre>
        </div>
        <div className="info-box">
          💡 如果您是<strong>初学者</strong>，建议先
          <button onClick={downloadSampleCSV} style={{ background: 'none', border: 'none', cursor: 'pointer', color: 'var(--blue)', fontWeight: 500, padding: '0 2px', fontFamily: 'inherit', fontSize: 'inherit' }}>下载示例 CSV</button>
          跑通整个流程，再替换成您自己的数据。
        </div>
      </div>

      {/* Results interpretation */}
      <div className="card tut-section">
        <div className="tut-section-title">📊 如何解读预测结果？</div>
        <div className="result-example">
          <div className="result-card">
            <div className="result-card-label">概率（Prob）</div>
            <div className="result-card-val">0.9432</div>
            <div className="result-card-desc">越接近 1 越好。该化合物逆转您表达特征的置信度。</div>
          </div>
          <div className="result-card">
            <div className="result-card-label">排名（Rank）</div>
            <div className="result-card-val">#1</div>
            <div className="result-card-desc">所有 ~2400 个化合物按概率排名，第 1 名为最优候选。</div>
          </div>
          <div className="result-card">
            <div className="result-card-label">Logit 值</div>
            <div className="result-card-val">2.84</div>
            <div className="result-card-desc">模型原始输出，经 sigmoid 转换即为概率，用于高级分析。</div>
          </div>
        </div>
        <div className="info-box" style={{ marginTop: 0 }}>
          🔬 <strong>下一步建议：</strong>
          优先对<strong>排名前 10–20 名</strong>的化合物进行湿实验验证。
          建议交叉比对 DrugBank 或 ChEMBL 数据库，了解化合物的已知作用机制。
        </div>
      </div>

      {/* FAQ */}
      <div className="card tut-section">
        <div className="tut-section-title">❓ 常见问题</div>
        <div className="faq-list">
          {[
            {
              q: '我的基因名用的是 Ensembl ID（如 ENSG00000141510），可以吗？',
              a: '可以。DrugReflector 会自动将 Ensembl ID 转换为 HGNC 格式。系统也支持小写、带 Affymetrix 后缀（_at）、带版本号（.v2）的基因名，会在运行时自动清洗。',
            },
            {
              q: '需要包含全部 978 个 landmark 基因吗？',
              a: '不强制要求，但覆盖率越高预测越准。模型会自动使用您数据中与 978 个 landmark 基因重叠的部分。您可以在结果页查看基因覆盖率报告。',
            },
            {
              q: '模型 checkpoint 文件在哪里下载？',
              a: '请从 Zenodo 下载（DOI: 10.5281/zenodo.16912444）。安装 zenodo-get 后运行：pip install zenodo-get，然后：zenodo_get --output-dir checkpoints 16912444',
            },
            {
              q: '预测结果大约需要多长时间？',
              a: '使用 CPU 通常 10–60 秒，取决于样本数量。如果您有 NVIDIA GPU，速度可提升 5–10 倍。',
            },
            {
              q: '可以同时预测多个样本吗？',
              a: '可以。在 CSV 文件中每行放一个样本，或上传包含多个 obs 的 H5AD 文件。结果页会提供下拉菜单切换不同样本的结果。',
            },
            {
              q: '预测结果概率都很低（接近 0），是数据有问题吗？',
              a: '概率是相对值，低概率不代表结果无效。关键是看相对排名——排名前列的化合物仍然是最佳候选。同时请检查：基因名是否为 HGNC 格式，数据是否已经过归一化处理。',
            },
          ].map(({ q, a }) => (
            <details key={q} className="faq-item">
              <summary>{q}</summary>
              <div className="faq-body">{a}</div>
            </details>
          ))}
        </div>
      </div>

      {/* CTA */}
      <div style={{ textAlign: 'center', padding: '20px 0 40px' }}>
        <div style={{ fontSize: '.88rem', color: 'var(--mid)', marginBottom: 16 }}>准备好了吗？</div>
        <button className="btn btn-primary" onClick={onStart} style={{ fontSize: '1rem', padding: '.7rem 2rem' }}>
          开始使用 DrugReflector →
        </button>
      </div>
    </>
  )
}

// ─── Upload step ──────────────────────────────────────────────────────────────

function UploadStep({
  onReady,
  onGeoReady,
  checkpoints,
  backendOnline,
  onTutorial,
  onRefreshCheckpoints,
  locale,
}: {
  onReady: (file: File) => void
  onGeoReady: (config: GeoRunConfig, preview: GeoPreviewResponse | null) => void
  checkpoints: CheckpointStatus | null
  backendOnline: boolean | null
  onTutorial: () => void
  onRefreshCheckpoints: () => void
  locale: Locale
}) {
  const isZh = locale === 'zh'
  const [mode, setMode] = useState<InputMode>('file')
  const [file, setFile] = useState<File | null>(null)
  const [refreshing, setRefreshing] = useState(false)
  const [geoAccession, setGeoAccession] = useState('')
  const [geoControl, setGeoControl] = useState('')
  const [geoCase, setGeoCase] = useState('')
  const [geoPreviewData, setGeoPreviewData] = useState<GeoPreviewResponse | null>(null)
  const [geoBusy, setGeoBusy] = useState(false)
  const [geoError, setGeoError] = useState<string | null>(null)
  const onDrop = useCallback((accepted: File[]) => { if (accepted[0]) setFile(accepted[0]) }, [])
  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: {
      'application/octet-stream': ['.h5ad'],
      'text/csv': ['.csv'],
      'text/tab-separated-values': ['.tsv'],
    },
    maxFiles: 1,
  })

  async function handleRefresh() {
    setRefreshing(true)
    onRefreshCheckpoints()
    await new Promise(r => setTimeout(r, 800))
    setRefreshing(false)
  }

  async function handleGeoPreview() {
    if (!geoAccession.trim()) return
    setGeoBusy(true)
    setGeoError(null)
    try {
      const data = await previewGeo(geoAccession.trim())
      setGeoPreviewData(data)
    } catch (e) {
      setGeoPreviewData(null)
      setGeoError((e as Error).message)
    } finally {
      setGeoBusy(false)
    }
  }

  function handleContinue() {
    if (mode === 'file' && file) {
      onReady(file)
      return
    }
    if (mode === 'geo' && geoAccession.trim()) {
      onGeoReady(
        {
          accession: geoAccession.trim().toUpperCase(),
          control_keyword: geoControl.trim() || undefined,
          case_keyword: geoCase.trim() || undefined,
        },
        geoPreviewData,
      )
    }
  }

  const cpNames = ['model_fold_0.pt', 'model_fold_1.pt', 'model_fold_2.pt']
  const allReady = checkpoints?.all_ready ?? false
  const checkpointDir = checkpoints?.checkpoint_dir || 'checkpoints'
  const downloadCmd = `zenodo_get --output-dir "${checkpointDir}" 16912444`
  const geoPreviewReadyForPrediction =
    !!geoPreviewData &&
    (
      !geoPreviewData.organism ||
      geoPreviewData.organism.toLowerCase() === 'homo sapiens' ||
      geoPreviewData.symbol_source === 'mouse_to_human_orthologs'
    ) &&
    geoPreviewData.symbol_source !== 'ensembl_ids'
  const canContinueFile = mode === 'file' && !!file && allReady && backendOnline !== false
  const canContinueGeo =
    mode === 'geo' &&
    !!geoAccession.trim() &&
    geoPreviewReadyForPrediction &&
    allReady &&
    backendOnline !== false
  const inputName = file?.name ?? ''
  const inputMeta = file ? fmtBytes(file.size) : ''
  const orthologRate =
    geoPreviewData?.ortholog_mapping && geoPreviewData.ortholog_mapping.input_genes > 0
      ? geoPreviewData.ortholog_mapping.mapped_input_genes / geoPreviewData.ortholog_mapping.input_genes
      : null

  return (
    <>
      <div className="card">
        <div className="card-title">{isZh ? '数据来源' : 'Data Source'}</div>
        <div className="mode-switch">
          <button
            type="button"
            className={`mode-chip ${mode === 'file' ? 'active' : ''}`}
            onClick={() => setMode('file')}
          >
            {isZh ? '上传文件' : 'Upload file'}
          </button>
          <button
            type="button"
            className={`mode-chip ${mode === 'geo' ? 'active' : ''}`}
            onClick={() => setMode('geo')}
          >
            {isZh ? 'GEO 编号' : 'GEO accession'}
          </button>
        </div>
        <div className="card-desc" style={{ marginBottom: 0 }}>
          {isZh
            ? '可以上传已经准备好的 signature 文件，或让系统自动下载公开 GEO 数据并为你构建 signature。'
            : 'Choose a prepared signature file, or let the app fetch a public GEO series and build the signature for you.'}
        </div>
      </div>

      {mode === 'geo' && (
        <div className="card">
          <div className="card-title">
            {isZh ? 'GEO 编号' : 'GEO accession'}
            <Tip text={isZh
              ? '例如 GSE6631。后端会下载 GEO 数据，在可能时把 probe 映射到基因名，自动识别两组样本，计算 signature，并运行 DrugReflector。'
              : 'Example: GSE6631. The backend downloads the GEO series matrix, maps probes to gene symbols when possible, auto-detects two groups, computes a signature, and runs DrugReflector.'} />
          </div>
          <div className="card-desc">
            {isZh
              ? '当 GEO 元数据不够清晰时，可以填写辅助关键词。必要时分别填一个对照组关键词和一个病例组关键词。'
              : 'Optional keywords help when GEO metadata is ambiguous. Use one control keyword and one case keyword if needed.'}
          </div>

          <div className="geo-grid">
            <div>
              <label className="field-label" htmlFor="geo-accession">{isZh ? 'GEO 编号' : 'GEO accession'}</label>
              <input
                id="geo-accession"
                type="text"
                placeholder="GSE6631"
                value={geoAccession}
                onChange={e => setGeoAccession(e.target.value.toUpperCase())}
              />
            </div>
            <div>
              <label className="field-label" htmlFor="geo-control">{isZh ? '对照组关键词' : 'Control keyword'}</label>
              <input
                id="geo-control"
                type="text"
                placeholder="normal / control"
                value={geoControl}
                onChange={e => setGeoControl(e.target.value)}
              />
            </div>
            <div>
              <label className="field-label" htmlFor="geo-case">{isZh ? '病例组关键词' : 'Case keyword'}</label>
              <input
                id="geo-case"
                type="text"
                placeholder="tumor / disease"
                value={geoCase}
                onChange={e => setGeoCase(e.target.value)}
              />
            </div>
          </div>

          <div className="nav-actions" style={{ marginTop: 16 }}>
            <button
              type="button"
              className="btn btn-secondary"
              disabled={!geoAccession.trim() || geoBusy}
              onClick={handleGeoPreview}
            >
              {geoBusy ? (isZh ? '正在检查 GEO…' : 'Inspecting GEO...') : (isZh ? '检查 GEO' : 'Inspect GEO')}
            </button>
          </div>

          {geoError && (
            <div className="error-box" style={{ marginTop: 14, marginBottom: 0 }}>
              <div className="error-title">{isZh ? 'GEO 获取失败' : 'GEO fetch failed'}</div>
              <div className="error-msg">{geoError}</div>
            </div>
          )}

          {geoPreviewData && (
            <div className="info-box">
              {geoPreviewData.organism?.toLowerCase() === 'mus musculus' && geoPreviewData.symbol_source === 'mouse_to_human_orthologs' && (
                <div className="info-box" style={{ marginBottom: 12 }}>
                  <div><strong>{isZh ? '当前数据物种：' : 'Current dataset species: '}</strong>{isZh ? '小鼠' : 'Mouse'}</div>
                  <div><strong>{isZh ? '预测所用物种体系：' : 'Prediction gene space: '}</strong>{isZh ? '人' : 'Human'}</div>
                  <div><strong>{isZh ? '同源映射状态：' : 'Ortholog mapping: '}</strong>{isZh ? '已自动完成小鼠→人同源基因映射' : 'Mouse genes were automatically mapped to human orthologs'}</div>
                  {orthologRate !== null && (
                    <div><strong>{isZh ? '映射成功率：' : 'Mapping success rate: '}</strong>{formatRatio(orthologRate)}</div>
                  )}
                </div>
              )}
              <div><strong>{geoPreviewData.accession}</strong> · {geoPreviewData.n_samples} {isZh ? '个样本' : 'samples'} · {geoPreviewData.n_genes} {isZh ? '个映射后基因' : 'mapped genes'}</div>
              <div>{isZh ? '平台：' : 'Platform: '}<code>{geoPreviewData.platform_id || (isZh ? '未知' : 'unknown')}</code></div>
              {geoPreviewData.organism && <div>{isZh ? '物种：' : 'Organism: '}<code>{getOrganismLabel(geoPreviewData.organism, locale)}</code></div>}
              {geoPreviewData.expression_source && <div>{isZh ? '表达数据来源：' : 'Expression source: '}<code>{getExpressionSourceLabel(geoPreviewData.expression_source, locale)}</code></div>}
              <div>{isZh ? '基因映射方式：' : 'Gene mapping: '}<code>{getGeneMappingLabel(geoPreviewData.symbol_source, locale)}</code></div>
              {geoPreviewData.ortholog_mapping && (
                <div>
                  {isZh
                    ? <>同源映射：<code>{geoPreviewData.ortholog_mapping.mapped_input_genes}</code> / <code>{geoPreviewData.ortholog_mapping.input_genes}</code> 个小鼠基因成功映射到 <code>{geoPreviewData.ortholog_mapping.unique_human_symbols}</code> 个人类基因符号</>
                    : <>Ortholog mapping: <code>{geoPreviewData.ortholog_mapping.mapped_input_genes}</code> / <code>{geoPreviewData.ortholog_mapping.input_genes}</code> mouse genes mapped to <code>{geoPreviewData.ortholog_mapping.unique_human_symbols}</code> human symbols</>}
                </div>
              )}
              <div>{isZh ? '是否进行了 log2 转换：' : 'Applied log2 transform: '}<code>{String(geoPreviewData.used_log2)}</code></div>
              {geoPreviewData.organism && geoPreviewData.organism.toLowerCase() !== 'homo sapiens' && geoPreviewData.symbol_source !== 'mouse_to_human_orthologs' && (
                <div className="warn-box" style={{ marginTop: 10 }}>
                  {isZh
                    ? '这个 GEO 数据集不是人类数据。当前自动 GEO 预测仅直接支持人类数据，除非基因已经客观映射为人类基因符号。'
                    : 'This GEO series is not human. Automatic DrugReflector GEO prediction currently only supports human datasets unless genes have already been mapped to human symbols objectively.'}
                </div>
              )}
              {geoPreviewData.symbol_source === 'mouse_to_human_orthologs' && (
                <div className="info-box" style={{ marginTop: 10 }}>
                  {isZh
                    ? '在预测前，系统已经依据官方 MGI 报告，把小鼠 Ensembl 基因 ID 客观映射成一对一的人类同源基因符号。'
                    : 'Mouse Ensembl gene IDs were objectively mapped to human one-to-one ortholog symbols using official MGI reports before prediction.'}
                </div>
              )}
              {geoPreviewData.symbol_source === 'ensembl_ids' && (
                <div className="warn-box" style={{ marginTop: 10 }}>
                  {isZh
                    ? '导入矩阵使用的是 Ensembl 风格基因 ID。当前自动 GEO 预测在解析后仍期望使用人类基因符号。'
                    : 'The imported matrix uses Ensembl-style gene IDs. Automatic GEO prediction currently expects human gene symbols after parsing.'}
                </div>
              )}
              {geoPreviewData.detected_grouping && (
                <div style={{ marginTop: 8 }}>
                  {isZh ? '自动分组：' : 'Auto grouping: '}<code>{geoPreviewData.detected_grouping.group_column}</code> :
                  <code> {geoPreviewData.detected_grouping.group1_value}</code> {isZh ? '对比' : 'vs'} <code>{geoPreviewData.detected_grouping.group2_value}</code>
                </div>
              )}
              {geoPreviewData.candidate_columns.length > 0 && (
                <div style={{ marginTop: 8 }}>
                  <div style={{ marginBottom: 6 }}>{isZh ? '候选元数据列' : 'Candidate metadata columns'}</div>
                  <div className="geo-candidates">
                    {geoPreviewData.candidate_columns.slice(0, 4).map(item => (
                      <span key={item.column} className="geo-candidate">
                        {item.column} ({item.n_unique})
                      </span>
                    ))}
                  </div>
                </div>
              )}
            </div>
          )}
        </div>
      )}

      {/* Upload card */}
      {mode === 'file' && <div className="card">
        <div className="card-title">
          📂 上传基因表达数据
          <Tip text="支持 H5AD（AnnData）或 CSV/TSV 格式。基因名需为 HGNC 格式，如 TP53、EGFR。" />
        </div>
        <div className="card-desc">
          拖放或点击上传包含 v-score 的文件。不确定格式？
          <button
            className="btn btn-sm"
            style={{ background: 'none', color: 'var(--blue)', padding: '0 4px', fontWeight: 500 }}
            onClick={onTutorial}
          >
            查看使用教程 →
          </button>
        </div>

        <div {...getRootProps()} className={`dropzone ${isDragActive ? 'active' : ''}`}>
          <input {...getInputProps()} />
          <div className="dropzone-icon">{isDragActive ? '📥' : '☁️'}</div>
          <div className="dropzone-text">
            {isDragActive
              ? <strong>松开即可上传</strong>
              : <><strong>点击选择文件</strong>，或拖放至此</>}
          </div>
          <div className="dropzone-hint">.h5ad · .csv · .tsv · 最大 200 MB</div>
        </div>

        <div className="info-box" style={{ marginTop: 14 }}>
          {isZh ? (
            <>
              <div><strong>单细胞上传支持：</strong>支持上传您自己的单细胞 <code>.h5ad</code> 文件。</div>
              <div>如果 <code>adata.obs</code> 里同时有“样本编号列”和“两组分组列”，系统会先按样本做伪批量，再客观计算 signature。</div>
              <div>如果缺少这两类元数据，页面会明确提示需要补什么，不会瞎编结果。</div>
              <div><strong>GEO 单细胞：</strong>目前优先支持上传整理后的 <code>.h5ad</code>。公开 GEO 的单细胞数据格式差异太大，只有在能客观解析出表达矩阵和分组信息时才会继续。</div>
            </>
          ) : (
            <>
              <div><strong>Single-cell uploads:</strong> your own single-cell <code>.h5ad</code> files are supported.</div>
              <div>If <code>adata.obs</code> contains both a sample identifier column and a two-group column, the app will pseudobulk cells by sample before computing the signature objectively.</div>
              <div>If that metadata is missing, the app will stop and tell you what is missing instead of inventing a result.</div>
              <div><strong>Single-cell GEO:</strong> uploaded <code>.h5ad</code> is the most reliable path. Public GEO single-cell formats vary widely, so the app only continues when an expression matrix and grouping can be parsed objectively.</div>
            </>
          )}
        </div>

        {file && (
          <div className="file-pill">
            <span>📄</span>
            <span className="file-pill-name">{inputName}</span>
            <span className="file-pill-size">{inputMeta}</span>
            <button className="file-pill-remove" onClick={() => setFile(null)}>×</button>
          </div>
        )}

        {!file && (
          <div className="info-box">
            💡 首次使用？
            <button
              onClick={downloadSampleCSV}
              style={{ background: 'none', border: 'none', cursor: 'pointer', color: 'var(--blue)', fontWeight: 500, padding: '0 2px', fontFamily: 'inherit', fontSize: 'inherit' }}
            >
              下载示例 CSV 文件
            </button>
            ，用随机数据跑通完整流程后，再替换成您自己的数据。
          </div>
        )}
      </div>}

      {/* Checkpoints card */}
      <div className="card">
        <div style={{ display: 'flex', alignItems: 'center', gap: 8, marginBottom: 5 }}>
          <div className="card-title" style={{ marginBottom: 0 }}>
            🗂 模型 Checkpoint 文件
            <Tip text="DrugReflector 使用 3 个神经网络集成模型，需要提前下载对应的权重文件放在 checkpoints/ 目录下。" />
          </div>
          <div style={{ flex: 1 }} />
          <button
            className="btn btn-sm btn-secondary"
            onClick={handleRefresh}
            disabled={refreshing}
            style={{ flexShrink: 0 }}
          >
            {refreshing ? '检测中…' : '🔄 刷新状态'}
          </button>
        </div>
        <div className="card-desc">
          需要将 3 个模型文件放在程序目录下的 <code>checkpoints/</code> 文件夹中。
        </div>
        {backendOnline === true && (
          <div style={{ fontSize: '.78rem', color: 'var(--mid)', marginTop: 8 }}>
            当前检测目录：<code>{checkpointDir}</code>
          </div>
        )}

        {/* Backend offline warning */}
        {backendOnline === false && (
          <div className="error-box" style={{ marginBottom: 14 }}>
            <div className="error-title">⚠️ 无法连接到后端服务</div>
            <div className="error-msg">
              请在终端中启动后端服务后，点击「刷新状态」：<br />
              <code style={{ background: 'rgba(0,0,0,.06)', borderRadius: 4, padding: '2px 6px', fontFamily: 'monospace', fontSize: '.82rem' }}>
                uvicorn server:app --reload --port 8000
              </code>
            </div>
          </div>
        )}

        {/* Checkpoint list — only show when backend is reachable */}
        {backendOnline !== false && (
          <div className="cp-list">
            {backendOnline === null
              ? <div style={{ fontSize: '.84rem', color: 'var(--mid)', padding: '8px 0' }}>正在检测…</div>
              : cpNames.map(name => {
                  const ok = checkpoints?.checkpoints[name] ?? false
                  return (
                    <div key={name} className="cp-item">
                      <div className={`cp-dot ${ok ? 'ok' : 'missing'}`} />
                      <span className="cp-name">{name}</span>
                      <span className={`cp-status ${ok ? 'ok' : 'missing'}`}>
                        {ok ? '✓ 已就绪' : '未找到'}
                      </span>
                    </div>
                  )
                })
            }
          </div>
        )}

        {backendOnline === true && !allReady && (
          <div className="info-box">
            <strong>下载模型文件（在终端中执行）：</strong><br />
            <code>pip install zenodo-get</code><br />
            <code>{downloadCmd}</code><br />
            <div style={{ display: 'flex', gap: 10, flexWrap: 'wrap', marginTop: 10 }}>
              <a
                href="https://doi.org/10.5281/zenodo.16912444"
                target="_blank"
                rel="noreferrer"
                style={{ color: 'var(--blue)', fontWeight: 600, textDecoration: 'none' }}
              >
                打开 Zenodo 下载页
              </a>
              <button
                type="button"
                onClick={() => copyText(`pip install zenodo-get\n${downloadCmd}`)}
                style={{ background: 'none', border: 'none', cursor: 'pointer', color: 'var(--blue)', fontWeight: 600, padding: 0, fontFamily: 'inherit', fontSize: '.84rem' }}
              >
                复制下载命令
              </button>
            </div>
            <span style={{ marginTop: 6, display: 'block', fontSize: '.8rem' }}>
              ✅ 下载完成后点击上方「刷新状态」按钮即可识别，无需刷新页面。
            </span>
          </div>
        )}

        {backendOnline === true && allReady && (
          <div style={{ marginTop: 12, fontSize: '.84rem', color: 'var(--green)', fontWeight: 500 }}>
            ✓ 全部 3 个模型文件已就绪，可以开始预测。
          </div>
        )}
      </div>

      <div className="nav-actions">
        <div style={{ flex: 1 }} />
        <button
          className="btn btn-primary"
          disabled={!(canContinueFile || canContinueGeo)}
          onClick={handleContinue}
        >
          下一步：配置参数 →
        </button>
      </div>
      {!(canContinueFile || canContinueGeo) && (
        <p style={{ fontSize: '.78rem', color: 'var(--mid)', textAlign: 'right', marginTop: 6 }}>
          {backendOnline === false
            ? '请先启动后端服务。'
            : mode === 'file' && !file
            ? '请先上传数据文件。'
            : '请下载模型文件后点击「刷新状态」。'}
        </p>
      )}
    </>
  )
}

// ─── Config step ──────────────────────────────────────────────────────────────

function ConfigStep({
  file,
  geoConfig,
  geoPreview,
  onRun,
  onBack,
  locale,
}: {
  file: File | null
  geoConfig: GeoRunConfig | null
  geoPreview: GeoPreviewResponse | null
  onRun: (nTop: number, prep?: PrepOptions) => void
  onBack: () => void
  locale: Locale
}) {
  const isZh = locale === 'zh'
  const [nTop, setNTop] = useState(50)
  const inputName = file ? file.name : geoConfig?.accession || 'Input'
  const inputMeta = file
    ? fmtBytes(file.size)
    : geoPreview
    ? `${geoPreview.n_samples} samples`
    : 'GEO series'
  const [preparation, setPreparation] = useState<PreparationSummary | null>(null)
  const [preparing, setPreparing] = useState(false)
  const [prepError, setPrepError] = useState<string | null>(null)
  const [groupColumn, setGroupColumn] = useState('')
  const [group1Value, setGroup1Value] = useState('')
  const [group2Value, setGroup2Value] = useState('')
  const [sampleIdColumn, setSampleIdColumn] = useState('')

  useEffect(() => {
    if (!file) {
      setPreparation(null)
      setPrepError(null)
      return
    }

    const prep: PrepOptions | undefined =
      groupColumn && group1Value && group2Value
        ? {
            group_column: groupColumn,
            group1_value: group1Value,
            group2_value: group2Value,
            sample_id_column: sampleIdColumn || undefined,
          }
        : sampleIdColumn
        ? { sample_id_column: sampleIdColumn }
        : undefined

    let cancelled = false
    setPreparing(true)
    setPrepError(null)
    prepareFileInput(file, prep)
      .then(summary => {
        if (!cancelled) setPreparation(summary)
      })
      .catch(err => {
        if (!cancelled) {
          setPreparation(null)
          setPrepError((err as Error).message)
        }
      })
      .finally(() => {
        if (!cancelled) setPreparing(false)
      })

    return () => {
      cancelled = true
    }
  }, [file, groupColumn, group1Value, group2Value, sampleIdColumn])

  useEffect(() => {
    if (!preparation || preparation.candidate_groups.length === 0) return
    if (groupColumn) return
    const best = preparation.candidate_groups[0]
    if (best.values.length === 2 && best.score >= 8) {
      setGroupColumn(best.column)
      setGroup1Value(best.values[0])
      setGroup2Value(best.values[1])
    }
  }, [preparation, groupColumn])

  useEffect(() => {
    if (!preparation?.candidate_sample_ids?.length) return
    if (sampleIdColumn) return
    setSampleIdColumn(preparation.used_sample_id_column || preparation.candidate_sample_ids[0].column)
  }, [preparation, sampleIdColumn])

  const selectedCandidate = preparation?.candidate_groups.find(item => item.column === groupColumn) || null
  const selectedSampleCandidate =
    preparation?.candidate_sample_ids?.find(item => item.column === sampleIdColumn) || null
  const prepReady = !!geoConfig || preparation?.status === 'ready'
  const currentPrepOptions =
    groupColumn && group1Value && group2Value
      ? {
          group_column: groupColumn,
          group1_value: group1Value,
          group2_value: group2Value,
          sample_id_column: sampleIdColumn || undefined,
        }
      : sampleIdColumn
      ? { sample_id_column: sampleIdColumn }
      : undefined

  return (
    <>
      <div className="card">
        <div className="card-title">⚙️ 预测参数设置</div>
        <div className="card-desc">确认参数后，点击「开始预测」运行集成模型。</div>

        <div style={{ marginBottom: 20 }}>
          <label className="field-label">当前文件</label>
          <div className="file-pill" style={{ marginTop: 4 }}>
            <span>{file ? '📄' : 'GEO'}</span>
            <span className="file-pill-name">{inputName}</span>
            <span className="file-pill-size">{inputMeta}</span>
          </div>
        </div>

        {file && (
          <>
            <div style={{ marginBottom: 18 }}>
              <label className="field-label">自动准备摘要</label>
              {preparing && (
                <div className="info-box" style={{ marginTop: 8 }}>
                  正在客观分析上传数据，并在可能的情况下自动整理成 DrugReflector 可用的签名格式。
                </div>
              )}
              {prepError && (
                <div className="error-box" style={{ marginTop: 8, marginBottom: 0 }}>
                  <div className="error-title">准备失败</div>
                  <div className="error-msg">{prepError}</div>
                </div>
              )}
              {!preparing && preparation && (
                <div className={preparation.status === 'ready' ? 'info-box' : 'warn-box'} style={{ marginTop: 8 }}>
                  <div><strong>识别模式：</strong><code>{preparation.mode}</code></div>
                  <div><strong>原始维度：</strong><code>{preparation.original_shape.join(' × ')}</code></div>
                  {preparation.prepared_shape && (
                    <div><strong>整理后维度：</strong><code>{preparation.prepared_shape.join(' × ')}</code></div>
                  )}
                  {preparation.notes.map(note => (
                    <div key={note}>{note}</div>
                  ))}
                  {preparation.top_genes.length > 0 && (
                    <div style={{ marginTop: 8 }}>
                      <strong>高权重基因预览：</strong>
                      <div className="prep-chip-row">
                        {preparation.top_genes.map(item => (
                          <span key={item.gene} className="prep-chip">
                            {item.gene} {item.score.toFixed(2)}
                          </span>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              )}
            </div>

            {preparation?.candidate_sample_ids && preparation.candidate_sample_ids.length > 0 && (
              <div style={{ marginBottom: 18 }}>
                <label className="field-label">{isZh ? '样本编号列' : 'Sample-id column'}</label>
                <select
                  className="sample-select"
                  value={sampleIdColumn}
                  onChange={e => setSampleIdColumn(e.target.value)}
                >
                  <option value="">{isZh ? '请选择样本编号列' : 'Select a sample-id column'}</option>
                  {preparation.candidate_sample_ids.map(item => (
                    <option key={item.column} value={item.column}>
                      {item.column} ({item.n_unique} {isZh ? '个样本' : 'samples'})
                    </option>
                  ))}
                </select>
                <p className="field-hint" style={{ marginTop: 6 }}>
                  {isZh
                    ? '单细胞数据会先按这个列做伪批量，再进入后续差异 signature 计算。'
                    : 'For single-cell data, cells are pseudobulked by this column before signature computation.'}
                  {selectedSampleCandidate ? ` ${isZh ? '当前候选评分' : 'Current heuristic score'}: ${selectedSampleCandidate.score}` : ''}
                </p>
              </div>
            )}

            {preparation?.status === 'needs_configuration' && preparation.candidate_groups.length > 0 && (
              <div className="prep-grid" style={{ marginBottom: 18 }}>
                <div>
                  <label className="field-label" htmlFor="group-column">分组列</label>
                  <select
                    id="group-column"
                    className="sample-select"
                    value={groupColumn}
                    onChange={e => {
                      const next = e.target.value
                      const candidate = preparation.candidate_groups.find(item => item.column === next)
                      setGroupColumn(next)
                      setGroup1Value(candidate?.values[0] || '')
                      setGroup2Value(candidate?.values[1] || '')
                    }}
                  >
                    <option value="">请选择元数据列</option>
                    {preparation.candidate_groups.map(item => (
                      <option key={item.column} value={item.column}>
                        {item.column} ({item.n_unique} groups)
                      </option>
                    ))}
                  </select>
                </div>
                <div>
                  <label className="field-label" htmlFor="group1-value">参考组</label>
                  <select
                    id="group1-value"
                    className="sample-select"
                    value={group1Value}
                    onChange={e => setGroup1Value(e.target.value)}
                    disabled={!selectedCandidate}
                  >
                    <option value="">请选择</option>
                    {(selectedCandidate?.values || []).map(value => (
                      <option key={value} value={value}>{value}</option>
                    ))}
                  </select>
                </div>
                <div>
                  <label className="field-label" htmlFor="group2-value">目标组</label>
                  <select
                    id="group2-value"
                    className="sample-select"
                    value={group2Value}
                    onChange={e => setGroup2Value(e.target.value)}
                    disabled={!selectedCandidate}
                  >
                    <option value="">请选择</option>
                    {(selectedCandidate?.values || []).map(value => (
                      <option key={value} value={value}>{value}</option>
                    ))}
                  </select>
                </div>
              </div>
            )}
          </>
        )}

        {geoConfig && geoPreview && (
          <div className="info-box" style={{ marginBottom: 18 }}>
            <div><strong>自动准备来源：</strong>GEO <code>{geoConfig.accession}</code></div>
            <div><strong>样本数：</strong>{geoPreview.n_samples}</div>
            {geoPreview.detected_grouping && (
              <div>
                <strong>识别到的分组：</strong><code>{geoPreview.detected_grouping.group_column}</code> :
                <code> {geoPreview.detected_grouping.group1_value}</code> vs
                <code> {geoPreview.detected_grouping.group2_value}</code>
              </div>
            )}
          </div>
        )}

        <div>
          <label className="field-label" htmlFor="ntop">
            返回候选化合物数量
            <Tip text="模型会对约 2,400 个化合物打分并排名，这里设置返回排名前 N 的化合物。建议设为 50–100。" />
          </label>
          <input
            id="ntop"
            type="number"
            min={10}
            max={2000}
            step={10}
            value={nTop}
            onChange={e => setNTop(Number(e.target.value))}
          />
          <p className="field-hint">
            模型将对约 2,400 个化合物打分，返回概率最高的 <strong>{nTop}</strong> 个。
          </p>
        </div>

        <div className="info-box" style={{ marginTop: 18 }}>
          💡 DrugReflector 使用 <strong>3 个神经网络</strong>（3 折交叉验证训练）的集成模型，
          对三个模型的预测结果取平均，以提升可靠性。
        </div>
      </div>

      <div className="nav-actions">
        <button className="btn btn-secondary" onClick={onBack}>← 重新上传</button>
        <div style={{ flex: 1 }} />
        <button
          className="btn btn-primary"
          disabled={!prepReady || preparing}
          onClick={() => onRun(nTop, currentPrepOptions)}
        >
          ⚡ 开始预测
        </button>
      </div>
      {!prepReady && !preparing && (
        <p style={{ fontSize: '.78rem', color: 'var(--mid)', textAlign: 'right', marginTop: 6 }}>
          当前数据还不能被客观地整理成可预测格式。请基于实际分组选择元数据列，或上传差异结果 / 签名文件。
        </p>
      )}
    </>
  )
}

// ─── Running step ─────────────────────────────────────────────────────────────

const RUNNING_MSGS = [
  '正在加载模型 checkpoint 文件…',
  '正在预处理基因名（HGNC 格式转换）…',
  '正在运行前向推理（3 个模型）…',
  '正在平均集成模型预测结果…',
  '正在整理输出格式…',
]

function RunningStep() {
  const [msgIdx, setMsgIdx] = useState(0)
  const timer = useRef<ReturnType<typeof setInterval> | null>(null)

  useEffect(() => {
    timer.current = setInterval(() => {
      setMsgIdx(i => Math.min(i + 1, RUNNING_MSGS.length - 1))
    }, 2400)
    return () => { if (timer.current) clearInterval(timer.current) }
  }, [])

  return (
    <div className="card">
      <div className="loading-wrap">
        <div className="spinner" />
        <div className="loading-title">DrugReflector 正在运行中</div>
        <div className="loading-sub">通常需要 10–60 秒，请耐心等待。</div>
        <ul className="loading-steps">
          {RUNNING_MSGS.map((msg, i) => (
            <li key={msg} className={i < msgIdx ? 'done' : i === msgIdx ? 'active' : ''}>
              <span>{i < msgIdx ? '✓' : i === msgIdx ? '⏳' : '○'}</span>
              {msg}
            </li>
          ))}
        </ul>
      </div>
    </div>
  )
}

// ─── Results step ─────────────────────────────────────────────────────────────

const BLUE_SCALE = ['#0071e3', '#1a8fff', '#34aadc', '#5bbce8', '#84cef0', '#b0d9f8']

function BarTip({ active, payload }: { active?: boolean; payload?: { payload: CompoundResult }[] }) {
  if (!active || !payload?.length) return null
  const d = payload[0].payload
  return (
    <div style={{
      background: '#fff', border: '1px solid var(--border)', borderRadius: 10,
      padding: '10px 14px', boxShadow: 'var(--shadow-md)', fontSize: '.82rem',
    }}>
      <div style={{ fontWeight: 600, marginBottom: 4 }}>{d.compound}</div>
      <div style={{ color: 'var(--mid)' }}>
        概率：<strong style={{ color: 'var(--blue)' }}>{d.prob.toFixed(4)}</strong>
      </div>
      <div style={{ color: 'var(--mid)' }}>排名：<strong>#{d.rank}</strong></div>
    </div>
  )
}

function ResultsStep({ response, onReset }: { response: PredictionResponse; onReset: () => void }) {
  const { data, meta } = response
  const [sample, setSample] = useState(data.samples[0])
  const [tab, setTab] = useState<'chart' | 'scatter' | 'table'>('chart')
  const [sortKey, setSortKey] = useState<SortKey>('rank')
  const [sortAsc, setSortAsc] = useState(true)
  const [nShow, setNShow] = useState(25)

  const raw = data.results[sample] ?? []
  const maxProb = [...raw].sort((a, b) => b.prob - a.prob)[0]?.prob ?? 1

  const sorted = [...raw].sort((a, b) => {
    const av = a[sortKey], bv = b[sortKey]
    const cmp = typeof av === 'string'
      ? av.localeCompare(bv as string)
      : (av as number) - (bv as number)
    return sortAsc ? cmp : -cmp
  })

  const topN = [...raw].sort((a, b) => b.prob - a.prob).slice(0, nShow)

  function toggleSort(k: SortKey) {
    if (k === sortKey) setSortAsc(p => !p)
    else { setSortKey(k); setSortAsc(k === 'rank') }
  }

  const Th = ({ k, label }: { k: SortKey; label: string }) => (
    <th onClick={() => toggleSort(k)}>
      {label} {sortKey === k ? (sortAsc ? '↑' : '↓') : ''}
    </th>
  )

  return (
    <>
      {/* Metrics */}
      <div className="metrics">
        <div className="metric">
          <div className="metric-val">{meta.n_samples}</div>
          <div className="metric-label">样本数</div>
        </div>
        <div className="metric">
          <div className="metric-val">{meta.n_compounds.toLocaleString()}</div>
          <div className="metric-label">候选化合物</div>
        </div>
        <div className="metric">
          <div className="metric-val">{meta.n_genes.toLocaleString()}</div>
          <div className="metric-label">输入基因数</div>
        </div>
      </div>

      {/* Sample selector */}
      {data.samples.length > 1 && (
        <div className="card" style={{ padding: '14px 22px' }}>
          <label className="field-label">当前查看的样本</label>
          <select className="sample-select" value={sample} onChange={e => setSample(e.target.value)}>
            {data.samples.map(s => <option key={s}>{s}</option>)}
          </select>
        </div>
      )}

      {/* Result card */}
      <div className="card">
        <div className="card-title">🎯 预测结果 — {sample}</div>
        <div className="card-desc">
          化合物按预测概率从高到低排列。概率越高，越可能逆转您的基因表达特征。
        </div>
        <div className="warn-box" style={{ marginTop: 0, marginBottom: 16 }}>
          ⚠️ <strong>重要提示：</strong>以下结果是模型给出的
          <strong>候选化合物优先级</strong>，不是实验验证结论。页面中的排名已按科研阅读习惯从
          <strong> 1 </strong>开始显示。
        </div>

        <div className="tabs">
          {([
            ['chart',   '📊 候选化合物'],
            ['scatter', '🔬 全局概览'],
            ['table',   '📋 完整列表'],
          ] as const).map(([t, label]) => (
            <button key={t} className={`tab ${tab === t ? 'active' : ''}`} onClick={() => setTab(t)}>
              {label}
            </button>
          ))}
        </div>

        {/* Bar chart */}
        {tab === 'chart' && (
          <>
            <div style={{ display: 'flex', alignItems: 'center', gap: 10, marginBottom: 14, flexWrap: 'wrap' }}>
              <span style={{ fontSize: '.8rem', color: 'var(--mid)' }}>显示前</span>
              {[15, 25, 50].map(n => (
                <button key={n}
                  className={`btn btn-sm ${nShow === n ? 'btn-primary' : 'btn-secondary'}`}
                  onClick={() => setNShow(n)}
                >
                  {n} 个
                </button>
              ))}
              <span style={{ fontSize: '.78rem', color: 'var(--mid)', marginLeft: 4 }}>
                颜色越深 = 排名越高
              </span>
            </div>
            <ResponsiveContainer width="100%" height={Math.max(300, topN.length * 30)}>
              <BarChart
                layout="vertical"
                data={[...topN].reverse()}
                margin={{ left: 10, right: 55, top: 0, bottom: 10 }}
              >
                <CartesianGrid strokeDasharray="3 3" stroke="#f0f0f5" horizontal={false} />
                <XAxis
                  type="number"
                  domain={[0, maxProb * 1.15]}
                  tickFormatter={v => v.toFixed(3)}
                  tick={{ fontSize: 11, fill: 'var(--mid)' }}
                  axisLine={false} tickLine={false}
                  label={{ value: '预测概率', position: 'insideBottom', offset: -4, fontSize: 11, fill: 'var(--mid)' }}
                />
                <YAxis
                  type="category" dataKey="compound" width={148}
                  tick={{ fontSize: 11, fill: 'var(--dark)' }}
                  axisLine={false} tickLine={false}
                />
                <Tooltip content={<BarTip />} cursor={{ fill: 'rgba(0,113,227,.05)' }} />
                <Bar dataKey="prob" radius={[0, 6, 6, 0]} maxBarSize={16}>
                  {[...topN].reverse().map((_, i) => {
                    const ci = Math.floor((i / Math.max(topN.length - 1, 1)) * (BLUE_SCALE.length - 1))
                    return <Cell key={i} fill={BLUE_SCALE[BLUE_SCALE.length - 1 - ci]} />
                  })}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </>
        )}

        {/* Scatter */}
        {tab === 'scatter' && (
          <>
            <div className="info-box" style={{ marginBottom: 14, marginTop: 0 }}>
              💡 蓝色点 = 排名前 10 的化合物，灰色点 = 其余化合物。
              左上角（排名低 + 概率高）的点是最佳候选。
            </div>
            <ResponsiveContainer width="100%" height={400}>
              <ScatterChart margin={{ left: 10, right: 20, top: 10, bottom: 30 }}>
                <CartesianGrid stroke="#f0f0f5" />
                <XAxis
                  type="number" dataKey="rank" name="排名"
                  label={{ value: '排名', position: 'insideBottom', offset: -16, fontSize: 12, fill: 'var(--mid)' }}
                  tick={{ fontSize: 11, fill: 'var(--mid)' }} axisLine={false} tickLine={false}
                />
                <YAxis
                  type="number" dataKey="prob" name="概率"
                  label={{ value: '预测概率', angle: -90, position: 'insideLeft', fontSize: 12, fill: 'var(--mid)' }}
                  tick={{ fontSize: 11, fill: 'var(--mid)' }} axisLine={false} tickLine={false}
                />
                <ZAxis range={[28, 28]} />
                <Tooltip
                  content={({ active, payload }) => {
                    if (!active || !payload?.length) return null
                    const d = payload[0].payload as CompoundResult
                    return (
                      <div style={{
                        background: '#fff', border: '1px solid var(--border)', borderRadius: 10,
                        padding: '10px 14px', boxShadow: 'var(--shadow-md)', fontSize: '.82rem',
                      }}>
                        <div style={{ fontWeight: 600, marginBottom: 4 }}>{d.compound}</div>
                        <div>概率：<strong style={{ color: 'var(--blue)' }}>{d.prob.toFixed(4)}</strong></div>
                        <div>排名：<strong>#{d.rank}</strong></div>
                      </div>
                    )
                  }}
                />
                <Scatter data={raw.filter(d => d.rank > 10)} fill="#d2d2d7" opacity={0.45} />
                <Scatter data={raw.filter(d => d.rank <= 10)} fill="#0071e3" opacity={0.9} />
              </ScatterChart>
            </ResponsiveContainer>
          </>
        )}

        {/* Table */}
        {tab === 'table' && (
          <>
            <div className="info-box" style={{ marginBottom: 14, marginTop: 0 }}>
              💡 点击列标题可排序。🥇 金色 = 前 3，🔵 蓝色 = 前 10。显示前 200 条。
            </div>
            <div className="table-wrap">
              <table>
                <thead>
                  <tr>
                    <Th k="rank" label="排名" />
                    <Th k="compound" label="化合物名称" />
                    <Th k="prob" label="预测概率" />
                    <Th k="logit" label="Logit 值" />
                  </tr>
                </thead>
                <tbody>
                  {sorted.slice(0, 200).map(row => {
                    const r = row.rank ?? 0
                    const bc = r <= 3 ? 'top3' : r <= 10 ? 'top10' : 'other'
                    return (
                      <tr key={row.compound}>
                        <td><span className={`rank-badge ${bc}`}>{r}</span></td>
                        <td style={{ fontWeight: 500 }}>{row.compound}</td>
                        <td>
                          <div className="prob-cell">
                            <span style={{ fontVariantNumeric: 'tabular-nums', fontSize: '.82rem', minWidth: 54 }}>
                              {row.prob.toFixed(4)}
                            </span>
                            <div className="prob-bar-track">
                              <div className="prob-bar-fill" style={{ width: `${(row.prob / maxProb) * 100}%` }} />
                            </div>
                          </div>
                        </td>
                        <td style={{ color: 'var(--mid)', fontVariantNumeric: 'tabular-nums', fontSize: '.82rem' }}>
                          {row.logit.toFixed(4)}
                        </td>
                      </tr>
                    )
                  })}
                </tbody>
              </table>
            </div>
          </>
        )}
      </div>

      {/* Actions */}
      <div className="dl-row">
        <button
          className="btn btn-primary"
          onClick={() => downloadCSV(sorted, `drugreflector_${sample}.csv`)}
        >
          ⬇ 下载 CSV 结果
        </button>
        <div style={{ flex: 1 }} />
        <button className="btn btn-secondary" onClick={onReset}>↺ 新建预测</button>
      </div>

      {/* Interpretation tips */}
      <details style={{ marginTop: 18 }}>
        <summary style={{ cursor: 'pointer', fontSize: '.86rem', color: 'var(--mid)', padding: '8px 0', userSelect: 'none' }}>
          ❓ 如何进一步使用这份结果？
        </summary>
        <div className="card" style={{ marginTop: 8 }}>
          <p style={{ fontSize: '.86rem', lineHeight: 1.8, color: 'var(--mid)' }}>
            <strong style={{ color: 'var(--dark)' }}>① 锁定候选化合物</strong><br />
            优先关注排名前 10–20 的化合物，对其进行湿实验验证（如细胞活力、基因表达回复实验）。<br /><br />
            <strong style={{ color: 'var(--dark)' }}>② 查询化合物背景信息</strong><br />
            将化合物名称在 <a href="https://www.drugbank.ca" style={{ color: 'var(--blue)' }}>DrugBank</a> 或{' '}
            <a href="https://www.ebi.ac.uk/chembl" style={{ color: 'var(--blue)' }}>ChEMBL</a> 中搜索，
            了解其已知的作用靶点、临床阶段和毒性信息。<br /><br />
            <strong style={{ color: 'var(--dark)' }}>③ 结合文献验证</strong><br />
            在 PubMed 中搜索「化合物名 + 您的疾病领域」，查看是否有已报道的体外或体内实验结果。
          </p>
        </div>
      </details>
    </>
  )
}

void ResultsStep

const EXPORT_BLUE_SCALE = ['#0b3d91', '#1358bf', '#1c78df', '#3d97ea', '#74b5f0', '#b8d5f5']

function FigureBarTip({ active, payload }: { active?: boolean; payload?: { payload: CompoundResult }[] }) {
  if (!active || !payload?.length) return null
  const item = payload[0].payload
  return (
    <div className="compound-tooltip">
      <div className="compound-tooltip-name">{getCompoundDisplayName(item)}</div>
      {getCompoundSecondaryLabel(item) && (
        <div className="compound-tooltip-id">{getCompoundSecondaryLabel(item)}</div>
      )}
      <div>Probability: <strong>{item.prob.toFixed(4)}</strong></div>
      <div>Rank: <strong>#{item.rank}</strong></div>
      {item.annotation?.target && <div>Target: <strong>{item.annotation.target}</strong></div>}
    </div>
  )
}

function EnhancedResultsStep({
  response,
  onReset,
  locale,
}: {
  response: PredictionResponse
  onReset: () => void
  locale: Locale
}) {
  const { data, meta } = response
  const [sample, setSample] = useState(data.samples[0])
  const isZh = locale === 'zh'
  const text = {
    publicationResults: isZh ? '科研结果视图' : 'Publication-style results',
    resultDescription: isZh
      ? '该面板按模型概率对候选化合物进行排序。这是基于表达 signature 的优先级筛选，不应解释为直接蛋白结合证据。'
      : 'This panel prioritizes compounds by model-assigned probability for the selected signature. It is a signature-based screening result and should not be interpreted as direct protein-binding evidence.',
    annotationPolicy: isZh
      ? '客观注释原则：名称、靶点、机制和结构仅在匹配到 Broad Repurposing Hub 官方注释及其关联 PubChem 记录时显示；没有可靠来源的字段会保留为空，不会猜测。'
      : 'Objective annotation policy: names, targets, mechanisms and structures are shown only when matched from official Broad Repurposing Hub annotations and linked PubChem records. Unmatched fields are left blank instead of being guessed.',
    inputQuality: isZh ? '输入质量' : 'Input quality',
    evidenceLimitations: isZh ? '证据与限制' : 'Evidence and limitations',
    annotationCoverage: isZh ? '注释覆盖度' : 'Annotation coverage',
    mechanismAggregation: isZh ? '机制聚合' : 'Mechanism aggregation',
    targetAggregation: isZh ? '靶点聚合' : 'Target aggregation',
    samples: isZh ? '样本数' : 'Samples',
    compoundsScreened: isZh ? '筛选化合物数' : 'Compounds screened',
    genesUsed: isZh ? '使用基因数' : 'Genes used',
    annotatedHits: isZh ? '客观注释命中' : 'Objectively annotated hits',
    selectSample: isZh ? '选择样本 / signature' : 'Select sample / signature',
    topCompounds: isZh ? '候选化合物' : 'Top compounds',
    globalView: isZh ? '全局视图' : 'Global view',
    annotatedTable: isZh ? '注释表格' : 'Annotated table',
    showTop: isZh ? '显示前' : 'Show top',
    directionColumn: isZh ? '方向性' : 'Direction',
    probability: isZh ? '预测概率' : 'Probability',
    logit: 'Logit',
    compound: isZh ? '化合物' : 'Compound',
    annotation: isZh ? '注释' : 'Annotation',
    selectedCompound: isZh ? '当前选中化合物' : 'Selected compound',
    newPrediction: isZh ? '重新开始预测' : 'New prediction',
    downloadSvg: 'Download SVG',
    downloadPng: 'Download PNG',
    downloadPdf: 'Download PDF',
    downloadCsv: 'Download CSV',
    downloadReport: isZh ? '下载中文报告' : 'Download report',
    noObjectiveEvidence: isZh ? '无客观证据' : 'No objective evidence',
    unavailable: isZh ? '不可用' : 'Unavailable',
    noSignedConnectivity: isZh ? '暂无外部签名连通性证据' : 'No signed external connectivity evidence',
    noObjectiveAnnotation: isZh ? '暂无客观注释' : 'No objective annotation',
    noObjectiveMoa: isZh ? '暂无客观 MOA 聚合信息。' : 'No objective MOA aggregation available.',
    noObjectiveTarget: isZh ? '暂无客观靶点聚合信息。' : 'No objective target aggregation available.',
    detailDirection: isZh ? '与输入 signature 的方向关系' : 'Direction vs input signature',
    directionSource: isZh ? '证据来源' : 'Source',
    directionScore: isZh ? '方向分数' : 'Score',
    queryGenes: isZh ? '方向查询基因' : 'Query genes for direction check',
    up: isZh ? '上调' : 'up',
    down: isZh ? '下调' : 'down',
    figureTitle: isZh
      ? `药物筛选结果摘要：${getDisplaySampleLabel(sample, locale)}`
      : `Drug screening summary: ${getDisplaySampleLabel(sample, locale)}`,
    figureKicker: isZh ? 'DrugReflector 结果图' : 'DrugReflector result figure',
    figureSubtitle: isZh
      ? `输入：${meta.filename} | 化合物：${meta.n_compounds.toLocaleString()} | 基因：${meta.n_genes.toLocaleString()}`
      : `Input: ${meta.filename} | Compounds: ${meta.n_compounds.toLocaleString()} | Genes: ${meta.n_genes.toLocaleString()}`,
  }
  const [tab, setTab] = useState<'chart' | 'scatter' | 'table'>('chart')
  const [sortKey, setSortKey] = useState<SortKey>('rank')
  const [sortAsc, setSortAsc] = useState(true)
  const [nShow, setNShow] = useState(25)
  const [isExporting, setIsExporting] = useState(false)
  const [selectedCompound, setSelectedCompound] = useState<CompoundResult | null>(null)
  const figureRef = useRef<HTMLDivElement | null>(null)

  const raw = data.results[sample] ?? []
  const topByProb = [...raw].sort((a, b) => b.prob - a.prob)
  const maxProb = topByProb[0]?.prob ?? 1
  const chartData = [...topByProb.slice(0, nShow)].reverse().map(item => ({
    ...item,
    label: getCompoundShortLabel(item, 30),
  }))
  const safeSample = sample.replace(/[^\w.-]+/g, '_')
  const annotatedCount = raw.filter(
    item => Boolean(item.annotation?.display_name || item.annotation?.target || item.annotation?.pubchem_cid),
  ).length
  const inputQuality = (meta.input_quality as {
    input_gene_count?: number
    model_gene_count?: number
    overlap_gene_count?: number
    missing_model_gene_count?: number
    overlap_ratio?: number
  } | undefined) ?? undefined
  const orthologMapping = (meta.ortholog_mapping as {
    source?: string
    input_genes?: number
    mapped_input_genes?: number
    unique_human_symbols?: number
    orthology_types?: string[]
  } | null | undefined) ?? null
  const topWindow = topByProb.slice(0, Math.min(20, topByProb.length))
  const topScatterRows = raw
    .filter(item => item.rank <= 8)
    .sort((a, b) => a.rank - b.rank)
    .map(item => ({
      ...item,
      scatterLabel: getCompoundShortLabel(item, 18),
    }))
  const scatterTopList = raw
    .filter(item => item.rank <= 10)
    .sort((a, b) => a.rank - b.rank)
  const moaSummary = summarizeTokens(topWindow.map(item => item.annotation?.moa))
  const targetSummary = summarizeTokens(topWindow.map(item => item.annotation?.target))
  const completenessSummary = topWindow.reduce(
    (acc, item) => {
      const tier = getAnnotationCompleteness(item).label
      acc[tier] = (acc[tier] ?? 0) + 1
      return acc
    },
    { High: 0, Medium: 0, Low: 0 } as Record<string, number>,
  )

  useEffect(() => {
    setSelectedCompound((data.results[sample] ?? [])[0] ?? null)
  }, [sample, data.results])

  const sorted = [...raw].sort((a, b) => {
    const av = a[sortKey]
    const bv = b[sortKey]
    const cmp = typeof av === 'string'
      ? av.localeCompare(bv as string)
      : (av as number) - (bv as number)
    return sortAsc ? cmp : -cmp
  })

  function toggleSort(k: SortKey) {
    if (k === sortKey) setSortAsc(prev => !prev)
    else {
      setSortKey(k)
      setSortAsc(k === 'rank')
    }
  }

  function selectCompound(value: unknown) {
    if (!value || typeof value !== 'object') return

    if ('payload' in value && value.payload && typeof value.payload === 'object') {
      setSelectedCompound(value.payload as CompoundResult)
      return
    }

    if ('compound' in value) {
      setSelectedCompound(value as CompoundResult)
    }
  }

  async function exportFigure(kind: 'png' | 'pdf' | 'svg') {
    if (!figureRef.current || !raw.length) return
    setIsExporting(true)
    try {
      if (kind === 'svg') {
        const dataUrl = await toSvg(figureRef.current, {
          cacheBust: true,
          backgroundColor: '#ffffff',
        })
        const link = document.createElement('a')
        link.href = dataUrl
        link.download = `drugreflector_${safeSample}_${tab}.svg`
        link.click()
        return
      }

      const dataUrl = await toPng(figureRef.current, {
        cacheBust: true,
        pixelRatio: 2.5,
        backgroundColor: '#ffffff',
      })

      if (kind === 'png') {
        const link = document.createElement('a')
        link.href = dataUrl
        link.download = `drugreflector_${safeSample}_${tab}.png`
        link.click()
        return
      }

      const image = new Image()
      image.src = dataUrl
      await new Promise<void>((resolve, reject) => {
        image.onload = () => resolve()
        image.onerror = () => reject(new Error('Could not render figure image for PDF export.'))
      })

      const orientation = image.width > image.height ? 'landscape' : 'portrait'
      const pdf = new jsPDF({ orientation, unit: 'mm', format: 'a4' })
      const pageWidth = pdf.internal.pageSize.getWidth()
      const pageHeight = pdf.internal.pageSize.getHeight()
      const margin = 12
      const scale = Math.min(
        (pageWidth - margin * 2) / image.width,
        (pageHeight - margin * 2) / image.height,
      )
      const width = image.width * scale
      const height = image.height * scale
      pdf.addImage(dataUrl, 'PNG', (pageWidth - width) / 2, (pageHeight - height) / 2, width, height)
      pdf.save(`drugreflector_${safeSample}_${tab}.pdf`)
    } finally {
      setIsExporting(false)
    }
  }

  const Th = ({ k, label }: { k: SortKey; label: string }) => (
    <th onClick={() => toggleSort(k)}>
      {label} {sortKey === k ? (sortAsc ? '↑' : '↓') : ''}
    </th>
  )

  const FigureScatterLabel = (props: any) => {
    const { cx, cy, fill, payload } = props
    if (!cx || !cy || !payload?.scatterLabel) return null
    const dy = payload.rank <= 3 ? -10 : payload.rank % 2 === 0 ? -8 : 15
    return (
      <g>
        <circle cx={cx} cy={cy} r={5} fill={fill} stroke="#ffffff" strokeWidth={1} />
        <text
          x={cx + 8}
          y={cy + dy}
          fontSize={11}
          fill="#23405f"
          fontWeight={600}
        >
          {payload.scatterLabel}
        </text>
      </g>
    )
  }

  const detail = selectedCompound ?? topByProb[0] ?? null
  const activeDirectionLabel = getDirectionDisplayLabel(detail?.direction, locale)
  const activeDirectionSummary = getDirectionSummary(detail?.direction, locale)
  const detailCompleteness = detail ? getAnnotationCompleteness(detail) : null

  return (
    <>
      <div className="metrics">
        <div className="metric">
          <div className="metric-val">{meta.n_samples}</div>
          <div className="metric-label">{text.samples}</div>
        </div>
        <div className="metric">
          <div className="metric-val">{meta.n_compounds.toLocaleString()}</div>
          <div className="metric-label">{text.compoundsScreened}</div>
        </div>
        <div className="metric">
          <div className="metric-val">{meta.n_genes.toLocaleString()}</div>
          <div className="metric-label">{text.genesUsed}</div>
        </div>
        <div className="metric">
          <div className="metric-val">{annotatedCount}</div>
          <div className="metric-label">{text.annotatedHits}</div>
        </div>
      </div>

      {data.samples.length > 1 && (
        <div className="card" style={{ padding: '14px 22px' }}>
          <label className="field-label">{text.selectSample}</label>
          <select className="sample-select" value={sample} onChange={e => setSample(e.target.value)}>
            {data.samples.map(item => <option key={item}>{getDisplaySampleLabel(item, locale)}</option>)}
          </select>
        </div>
      )}

      <div className="card">
        <div className="card-title">{text.publicationResults}</div>
        <div className="card-desc">{text.resultDescription}</div>
        <div className="warn-box" style={{ marginTop: 0, marginBottom: 16 }}>
          <strong>{locale === 'zh' ? '客观注释原则：' : 'Objective annotation policy: '}</strong>
          {text.annotationPolicy}
        </div>
        <div className="warn-box" style={{ marginTop: 0, marginBottom: 16 }}>
          <strong>{locale === 'zh' ? '方向性判读原则：' : 'Direction interpretation policy: '}</strong>
          {locale === 'zh'
            ? '这里只判断药物相对输入 signature 更像逆转、更像增强，还是暂无客观证据；这不直接等于某个目标基因已经被上调或下调，最终仍需实验验证。'
            : 'Direction is limited to reverse, mimic, or no objective evidence relative to the input signature; it does not directly mean a target gene has been up- or down-regulated, and final confirmation still requires experiments.'}
        </div>

        <div className="results-grid">
          <div className="results-panel">
            <div className="results-panel-title">{text.inputQuality}</div>
            <div className="results-panel-body">
              <div><strong>{locale === 'zh' ? '与模型基因重叠：' : 'Gene overlap with model:'}</strong> {inputQuality?.overlap_gene_count ?? text.unavailable} / {inputQuality?.model_gene_count ?? text.unavailable} ({formatRatio(inputQuality?.overlap_ratio)})</div>
              <div><strong>{locale === 'zh' ? '输入基因数：' : 'Input genes detected:'}</strong> {inputQuality?.input_gene_count ?? text.unavailable}</div>
              <div><strong>{locale === 'zh' ? '缺失模型基因：' : 'Missing model genes:'}</strong> {inputQuality?.missing_model_gene_count ?? text.unavailable}</div>
              {orthologMapping && (
                <>
                  <div><strong>{locale === 'zh' ? '跨物种映射：' : 'Cross-species mapping:'}</strong> {orthologMapping.source || text.unavailable}</div>
                  <div><strong>{locale === 'zh' ? '成功映射基因：' : 'Mapped genes:'}</strong> {orthologMapping.mapped_input_genes ?? text.unavailable} / {orthologMapping.input_genes ?? text.unavailable}</div>
                  <div><strong>{locale === 'zh' ? '唯一人类符号：' : 'Unique human symbols:'}</strong> {orthologMapping.unique_human_symbols ?? text.unavailable}</div>
                </>
              )}
            </div>
          </div>

          <div className="results-panel">
            <div className="results-panel-title">{text.evidenceLimitations}</div>
            <div className="results-panel-body">
              <div><strong>{locale === 'zh' ? '筛选类型：' : 'Screen type:'}</strong> {locale === 'zh' ? '基于 signature 的优先级筛选，不是直接结合证据。' : 'Signature-based prioritization, not direct binding evidence.'}</div>
              <div><strong>{locale === 'zh' ? '方向状态：' : 'Direction status:'}</strong> {activeDirectionLabel}</div>
              <div><strong>{locale === 'zh' ? '方向来源：' : 'Direction source:'}</strong> {detail?.direction?.source || text.unavailable}</div>
              <div><strong>{locale === 'zh' ? '当前限制：' : 'Current limitation:'}</strong> {detail?.direction?.reason || (locale === 'zh' ? '暂无额外签名连通性证据。' : 'No additional signed connectivity evidence is available.')}</div>
              <div><strong>{locale === 'zh' ? '方向性说明：' : 'Direction interpretation:'}</strong> {activeDirectionSummary}</div>
              <div><strong>{locale === 'zh' ? '结构说明：' : 'Structure caveat:'}</strong> {locale === 'zh' ? '导出的 SVG 对图表布局仍是矢量，但 PubChem 结构图本身是位图来源。' : 'Exported SVG remains vector for the chart layout, but PubChem structure images are raster source assets.'}</div>
            </div>
          </div>

          <div className="results-panel">
            <div className="results-panel-title">{text.annotationCoverage}</div>
            <div className="results-panel-body">
              <div><strong>{locale === 'zh' ? 'Top 20 高覆盖：' : 'Top 20 high coverage:'}</strong> {completenessSummary.High}</div>
              <div><strong>{locale === 'zh' ? 'Top 20 中覆盖：' : 'Top 20 medium coverage:'}</strong> {completenessSummary.Medium}</div>
              <div><strong>{locale === 'zh' ? 'Top 20 低覆盖：' : 'Top 20 low coverage:'}</strong> {completenessSummary.Low}</div>
              <div><strong>{locale === 'zh' ? '客观注释命中：' : 'Objectively annotated hits:'}</strong> {annotatedCount} / {raw.length}</div>
            </div>
          </div>
        </div>

        <div className="results-grid">
          <div className="results-panel">
            <div className="results-panel-title">{text.mechanismAggregation}</div>
            <div className="token-list">
              {moaSummary.length ? moaSummary.map(([token, count]) => (
                <span key={token} className="token-pill">{token} ({count})</span>
              )) : <span className="direction-note">{text.noObjectiveMoa}</span>}
            </div>
          </div>

          <div className="results-panel">
            <div className="results-panel-title">{text.targetAggregation}</div>
            <div className="token-list">
              {targetSummary.length ? targetSummary.map(([token, count]) => (
                <span key={token} className="token-pill">{token} ({count})</span>
              )) : <span className="direction-note">{text.noObjectiveTarget}</span>}
            </div>
          </div>
        </div>

        <div className="tabs">
          {([
            ['chart', text.topCompounds],
            ['scatter', text.globalView],
            ['table', text.annotatedTable],
          ] as const).map(([value, label]) => (
            <button
              key={value}
              className={`tab ${tab === value ? 'active' : ''}`}
              onClick={() => setTab(value)}
            >
              {label}
            </button>
          ))}
        </div>

        <div className="result-toolbar">
          <div className="result-toolbar-group">
            {tab === 'chart' && (
              <>
                <span className="result-toolbar-label">{text.showTop}</span>
                {[15, 25, 50].map(n => (
                  <button
                    key={n}
                    className={`btn btn-sm ${nShow === n ? 'btn-primary' : 'btn-secondary'}`}
                    onClick={() => setNShow(n)}
                  >
                    {n}
                  </button>
                ))}
              </>
            )}
            {tab === 'scatter' && (
              <span className="result-toolbar-label">
                {locale === 'zh' ? '蓝点表示预测概率排名前 10 的化合物。' : 'Blue points represent the top 10 compounds by predicted probability.'}
              </span>
            )}
            {tab === 'table' && (
              <span className="result-toolbar-label">
                {locale === 'zh' ? '点击行可查看化合物详情，点击列表头可排序。' : 'Click a row to inspect one compound in detail. Column headers sort the table.'}
              </span>
            )}
          </div>
          <div className="result-toolbar-group">
            <button
              className="btn btn-secondary"
              onClick={() => exportFigure('svg')}
              disabled={isExporting || !raw.length}
            >
              {locale === 'zh' ? '下载 SVG' : text.downloadSvg}
            </button>
            <button
              className="btn btn-secondary"
              onClick={() => exportFigure('png')}
              disabled={isExporting || !raw.length}
            >
              {locale === 'zh' ? '下载 PNG' : text.downloadPng}
            </button>
            <button
              className="btn btn-secondary"
              onClick={async () => {
                try {
                  setIsExporting(true)
                  await downloadWordReport(response, sample, locale)
                } finally {
                  setIsExporting(false)
                }
              }}
              disabled={isExporting || !raw.length}
            >
              {text.downloadReport}
            </button>
            <button
              className="btn btn-primary"
              onClick={() => downloadCSV(sorted, `drugreflector_${safeSample}.csv`)}
            >
              {locale === 'zh' ? '下载 CSV' : text.downloadCsv}
            </button>
          </div>
        </div>

        <div ref={figureRef} className="figure-sheet">
          <div className="figure-sheet-header">
            <div>
              <div className="figure-sheet-kicker">{text.figureKicker}</div>
              <div className="figure-sheet-title">{text.figureTitle}</div>
              <div className="figure-sheet-subtitle">{text.figureSubtitle}</div>
            </div>
          </div>

          {tab === 'chart' && (
            <>
              <ResponsiveContainer width="100%" height={Math.max(360, chartData.length * 31)}>
                <BarChart layout="vertical" data={chartData} margin={{ left: 8, right: 30, top: 10, bottom: 18 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#eef1f6" horizontal={false} />
                  <XAxis
                    type="number"
                    domain={[0, maxProb * 1.15]}
                    tickFormatter={value => value.toFixed(3)}
                    tick={{ fontSize: 11, fill: '#596273' }}
                    axisLine={false}
                    tickLine={false}
                    label={{ value: locale === 'zh' ? '预测概率' : 'Predicted probability', position: 'insideBottom', offset: -6, fontSize: 11, fill: '#596273' }}
                  />
                  <YAxis
                    type="category"
                    dataKey="label"
                    width={220}
                    tick={{ fontSize: 11, fill: '#1d1d1f' }}
                    axisLine={false}
                    tickLine={false}
                  />
                  <Tooltip content={<FigureBarTip />} cursor={{ fill: 'rgba(11,61,145,.05)' }} />
                  <Bar dataKey="prob" radius={[0, 7, 7, 0]} maxBarSize={16} onClick={selectCompound}>
                    {chartData.map((_, index) => {
                      const colorIndex = Math.floor((index / Math.max(chartData.length - 1, 1)) * (EXPORT_BLUE_SCALE.length - 1))
                      return <Cell key={index} fill={EXPORT_BLUE_SCALE[EXPORT_BLUE_SCALE.length - 1 - colorIndex]} />
                    })}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
              <div className="figure-sheet-caption">
                {locale === 'zh'
                  ? '化合物标签优先显示有客观来源的注释名称；如果没有，就保留原始 BRD 编号。'
                  : 'Compound labels prefer objective annotated names when available; otherwise the original BRD identifier is preserved.'}
              </div>
            </>
          )}

          {tab === 'scatter' && (
            <>
              <ResponsiveContainer width="100%" height={430}>
                <ScatterChart margin={{ left: 10, right: 24, top: 14, bottom: 32 }}>
                  <CartesianGrid stroke="#eef1f6" />
                  <XAxis
                    type="number"
                    dataKey="rank"
                    name={locale === 'zh' ? '排名' : 'Rank'}
                    label={{ value: locale === 'zh' ? '排名' : 'Rank', position: 'insideBottom', offset: -18, fontSize: 12, fill: '#596273' }}
                    tick={{ fontSize: 11, fill: '#596273' }}
                    axisLine={false}
                    tickLine={false}
                  />
                  <YAxis
                    type="number"
                    dataKey="prob"
                    name={locale === 'zh' ? '预测概率' : 'Predicted probability'}
                    label={{ value: locale === 'zh' ? '预测概率' : 'Predicted probability', angle: -90, position: 'insideLeft', fontSize: 12, fill: '#596273' }}
                    tick={{ fontSize: 11, fill: '#596273' }}
                    axisLine={false}
                    tickLine={false}
                  />
                  <ZAxis range={[28, 28]} />
                  <Tooltip content={<FigureBarTip />} />
                  <Scatter data={raw.filter(item => item.rank > 10)} fill="#ccd5e1" opacity={0.5} onClick={selectCompound} />
                  <Scatter data={raw.filter(item => item.rank <= 10)} fill="#1358bf" opacity={0.95} onClick={selectCompound} />
                  <Scatter data={topScatterRows} fill="#1358bf" shape={<FigureScatterLabel />} onClick={selectCompound} />
                </ScatterChart>
              </ResponsiveContainer>
              <div className="scatter-insight-grid">
                {scatterTopList.map(item => (
                  <div key={item.compound} className="scatter-insight-card">
                    <div className="scatter-insight-rank">#{item.rank}</div>
                    <div className="scatter-insight-name">{getCompoundDisplayName(item)}</div>
                    <div className="scatter-insight-meta">
                      {locale === 'zh' ? '预测概率' : 'Probability'} {item.prob.toFixed(4)}
                    </div>
                  </div>
                ))}
              </div>
              <div className="figure-sheet-caption">
                {locale === 'zh'
                  ? '点击柱子或散点即可更新下方详情卡。越靠左上通常表示排名更靠前且模型概率更高。'
                  : 'Leading compounds are labeled directly in the scatter plot, and the Top 10 list below keeps names and probabilities readable after export.'}
              </div>
            </>
          )}

          {tab === 'table' && (
            <>
              <div className="figure-table-wrap">
                <table>
                  <thead>
                    <tr>
                      <Th k="rank" label={locale === 'zh' ? '排名' : 'Rank'} />
                      <Th k="compound" label={text.compound} />
                      <th>{text.directionColumn}</th>
                      <Th k="prob" label={text.probability} />
                      <Th k="logit" label={text.logit} />
                      <th>{text.annotation}</th>
                    </tr>
                  </thead>
                  <tbody>
                    {sorted.slice(0, 200).map(row => {
                      const rank = row.rank ?? 0
                      const badgeClass = rank <= 3 ? 'top3' : rank <= 10 ? 'top10' : 'other'
                      const completeness = getAnnotationCompleteness(row)
                      return (
                        <tr key={row.compound} className="clickable-row" onClick={() => setSelectedCompound(row)}>
                          <td><span className={`rank-badge ${badgeClass}`}>{rank}</span></td>
                          <td>
                            <div className="compound-name-cell">
                              {getCompoundDisplayName(row)}
                              <span className={`annotation-badge ${completeness.className}`}>{completeness.label}</span>
                            </div>
                            {getCompoundSecondaryLabel(row) && (
                              <div className="compound-id-cell">{getCompoundSecondaryLabel(row)}</div>
                            )}
                          </td>
                          <td className="compound-meta-cell">
                            <div className={`direction-pill ${getDirectionClass(row.direction)}`}>
                              {getDirectionDisplayLabel(row.direction, locale)}
                            </div>
                            <div className="direction-note">
                              {row.direction?.source || text.noSignedConnectivity}
                            </div>
                          </td>
                          <td>
                            <div className="prob-cell">
                              <span style={{ fontVariantNumeric: 'tabular-nums', fontSize: '.82rem', minWidth: 54 }}>
                                {row.prob.toFixed(4)}
                              </span>
                              <div className="prob-bar-track">
                                <div className="prob-bar-fill" style={{ width: `${(row.prob / maxProb) * 100}%` }} />
                              </div>
                            </div>
                          </td>
                          <td style={{ color: 'var(--mid)', fontVariantNumeric: 'tabular-nums', fontSize: '.82rem' }}>
                            {row.logit.toFixed(4)}
                          </td>
                          <td className="compound-meta-cell">
                            {row.annotation?.target || row.annotation?.moa
                              ? `${row.annotation?.target || (locale === 'zh' ? '靶点不可用' : 'Target unavailable')} | ${row.annotation?.moa || (locale === 'zh' ? '机制不可用' : 'MOA unavailable')}`
                              : text.noObjectiveAnnotation}
                          </td>
                        </tr>
                      )
                    })}
                  </tbody>
                </table>
              </div>
              <div className="figure-sheet-caption">
                {locale === 'zh'
                  ? '如果 BRD 编号没有匹配到 Broad Repurposing Hub 官方记录，对应表格字段会保留为空。'
                  : 'Table fields remain blank when no official Broad Repurposing Hub match is available for a BRD identifier.'}
              </div>
            </>
          )}
        </div>
      </div>

      {detail && (
        <div className="card compound-detail-card">
          <div className="card-title">{text.selectedCompound}</div>
          <div className="compound-detail-grid">
            <div className="compound-detail-main">
              <div className="compound-detail-name">
                {getCompoundDisplayName(detail)}
                {detailCompleteness && (
                  <span className={`annotation-badge ${detailCompleteness.className}`}>{detailCompleteness.label}</span>
                )}
              </div>
              <div className="compound-detail-id">{detail.compound}</div>
              <div className="compound-detail-stats">
                <span>{locale === 'zh' ? `排名 #${detail.rank}` : `Rank #${detail.rank}`}</span>
                <span>{locale === 'zh' ? `概率 ${detail.prob.toFixed(4)}` : `Probability ${detail.prob.toFixed(4)}`}</span>
                <span>{text.logit} {detail.logit.toFixed(4)}</span>
              </div>
              <div className="compound-direction-block">
                <div className={`direction-pill ${getDirectionClass(detail.direction)}`}>
                  {text.detailDirection}: {getDirectionDisplayLabel(detail.direction, locale)}
                </div>
                <div className="direction-note">
                  {text.directionScore}: {detail.direction?.score ?? text.unavailable} | {text.directionSource}: {detail.direction?.source || text.unavailable}
                </div>
                <div className="direction-note">
                  {detail.direction?.reason || (locale === 'zh' ? '该化合物暂无客观签名连通性证据。' : 'No objective signed connectivity evidence is available for this compound.')}
                </div>
                <div className="direction-note">
                  {getDirectionSummary(detail.direction, locale)}
                </div>
                <div className="direction-note">
                  {text.queryGenes}: {text.up} {detail.direction?.query_up_genes ?? 0}, {text.down} {detail.direction?.query_down_genes ?? 0}
                </div>
              </div>
              <div className="compound-detail-meta">
                <div><strong>Annotated name:</strong> {detail.annotation?.display_name || 'Unavailable from objective source'}</div>
                <div><strong>Alternate label:</strong> {detail.annotation?.vendor_name || 'Unavailable'}</div>
                <div><strong>Trade / brand name:</strong> {detail.annotation?.trade_name || 'No reliable structured trade-name source available'}</div>
                <div><strong>Chemical name:</strong> {detail.annotation?.chemical_name || 'Unavailable'}</div>
                <div><strong>PubChem title:</strong> {detail.annotation?.pubchem_title || 'Unavailable'}</div>
                <div><strong>Molecular formula:</strong> {detail.annotation?.molecular_formula || 'Unavailable'}</div>
                <div><strong>Molecular weight:</strong> {detail.annotation?.molecular_weight ? detail.annotation.molecular_weight.toFixed(4) : 'Unavailable'}</div>
                <div><strong>Target:</strong> {detail.annotation?.target || 'Unavailable'}</div>
                <div><strong>Mechanism:</strong> {detail.annotation?.moa || 'Unavailable'}</div>
                <div><strong>Clinical phase:</strong> {detail.annotation?.clinical_phase || 'Unavailable'}</div>
                <div><strong>Disease area:</strong> {detail.annotation?.disease_area || 'Unavailable'}</div>
                <div><strong>Indication:</strong> {detail.annotation?.indication || 'Unavailable'}</div>
                <div><strong>SMILES:</strong> <span className="compound-mono">{detail.annotation?.smiles || 'Unavailable'}</span></div>
                <div><strong>Canonical SMILES:</strong> <span className="compound-mono">{detail.annotation?.canonical_smiles || 'Unavailable'}</span></div>
                <div><strong>Isomeric SMILES:</strong> <span className="compound-mono">{detail.annotation?.isomeric_smiles || 'Unavailable'}</span></div>
                <div><strong>InChIKey:</strong> <span className="compound-mono">{detail.annotation?.inchikey || 'Unavailable'}</span></div>
                <div><strong>Synonyms:</strong> {detail.annotation?.synonyms?.length ? detail.annotation.synonyms.join(' | ') : 'Unavailable'}</div>
                <div><strong>Source:</strong> {detail.annotation?.source || 'No objective annotation match'}</div>
                {detail.annotation?.pubchem_url && (
                  <div>
                    <strong>PubChem:</strong>{' '}
                    <a href={detail.annotation.pubchem_url} target="_blank" rel="noreferrer">
                      Open compound record
                    </a>
                  </div>
                )}
              </div>
            </div>
            <div className="compound-detail-structure">
              {detail.annotation?.structure_image ? (
                <>
                  <img
                    src={detail.annotation.structure_image}
                    alt={`Structure of ${getCompoundDisplayName(detail)}`}
                    className="compound-structure-image"
                    loading="lazy"
                    crossOrigin="anonymous"
                  />
                  <div className="compound-structure-caption">
                    {locale === 'zh'
                      ? '仅当存在客观的 PubChem CID 交叉引用时，才显示该 2D 结构图。'
                      : '2D structure rendered from PubChem only when an objective CID cross-reference is available.'}
                  </div>
                </>
              ) : (
                <div className="compound-structure-empty">
                  {locale === 'zh' ? '该 BRD 编号暂无客观结构图可用。' : 'No objective structure image is available for this BRD identifier.'}
                </div>
              )}
            </div>
          </div>
        </div>
      )}

      <div className="dl-row">
        <button className="btn btn-secondary" onClick={onReset}>{text.newPrediction}</button>
      </div>
    </>
  )
}

// ─── Error box ────────────────────────────────────────────────────────────────

function ErrorBox({ message, onDismiss }: { message: string; onDismiss: () => void }) {
  return (
    <div className="error-box">
      <div className="error-title">⚠️ 预测失败</div>
      <div className="error-msg">{message}</div>
      <details className="error-detail">
        <summary>常见解决方法</summary>
        <pre>{`• 确认 checkpoints/ 目录下的 3 个模型文件完整
• 基因名请使用 HGNC 格式（如 TP53、EGFR、CDKN1A）
• CSV 格式：首列为样本名，其余列为基因
• 如有 GPU 相关报错，可以换用 CPU`}</pre>
      </details>
      <button className="btn btn-sm btn-secondary" onClick={onDismiss} style={{ marginTop: 12 }}>
        返回修改设置
      </button>
    </div>
  )
}

// ─── App ─────────────────────────────────────────────────────────────────────

export default function App() {
  const [locale, setLocale] = useState<Locale>(() => {
    if (typeof window === 'undefined') return 'zh'
    return window.localStorage.getItem('drugreflector_locale') === 'en' ? 'en' : 'zh'
  })
  const [page, setPage] = useState<Page>('tool')
  const [step, setStep] = useState<Step>('upload')
  const [file, setFile] = useState<File | null>(null)
  const [geoConfig, setGeoConfig] = useState<GeoRunConfig | null>(null)
  const [geoPreviewSelection, setGeoPreviewSelection] = useState<GeoPreviewResponse | null>(null)
  const [checkpoints, setCheckpoints] = useState<CheckpointStatus | null>(null)
  const [backendOnline, setBackendOnline] = useState<boolean | null>(null)
  const [results, setResults] = useState<PredictionResponse | null>(null)
  const [error, setError] = useState<string | null>(null)

  const refreshCheckpoints = useCallback(() => {
    setBackendOnline(null)
    fetchCheckpoints()
      .then(data => { setCheckpoints(data); setBackendOnline(true) })
      .catch(() => { setCheckpoints(null); setBackendOnline(false) })
  }, [])

  useEffect(() => { refreshCheckpoints() }, [])
  useEffect(() => {
    if (typeof window !== 'undefined') {
      window.localStorage.setItem('drugreflector_locale', locale)
    }
  }, [locale])

  async function handleRun(nTop: number, prep?: PrepOptions) {
    if (!file && !geoConfig) return
    setStep('running')
    setError(null)
    try {
      const res = geoConfig
        ? await runGeoPrediction(geoConfig, nTop)
        : await runPrediction(file as File, nTop, prep)
      setResults(res)
      setStep('results')
    } catch (e) {
      setError((e as Error).message)
      setStep('config')
    }
  }

  function handleReset() {
    setStep('upload')
    setFile(null)
    setGeoConfig(null)
    setGeoPreviewSelection(null)
    setResults(null)
    setError(null)
    fetchCheckpoints().then(setCheckpoints).catch(() => null)
  }

  return (
    <div className="app">
      {/* Nav */}
      <nav className="nav">
        <div className="container nav-inner">
          <button className="nav-logo" onClick={() => setPage('tool')} style={{ background: 'none', border: 'none' }}>
            🧬 DrugReflector
          </button>
          <div className="nav-spacer" />
          <div className="lang-switch" role="group" aria-label="Language switch">
            <button
              className={`lang-chip ${locale === 'zh' ? 'active' : ''}`}
              onClick={() => setLocale('zh')}
            >
              中文
            </button>
            <button
              className={`lang-chip ${locale === 'en' ? 'active' : ''}`}
              onClick={() => setLocale('en')}
            >
              EN
            </button>
          </div>
          <button
            className={`nav-tab ${page === 'tutorial' ? 'active' : ''}`}
            onClick={() => setPage('tutorial')}
          >
            {locale === 'zh' ? '使用教程' : 'Tutorial'}
          </button>
          <button
            className={`nav-tab ${page === 'tool' ? 'primary' : ''}`}
            onClick={() => setPage('tool')}
          >
            {locale === 'zh' ? '开始使用 →' : 'Start →'}
          </button>
        </div>
      </nav>

      <main className="container">
        {/* ── Tutorial page ── */}
        {page === 'tutorial' && (
          <div style={{ paddingTop: 32 }}>
            <TutorialPage onStart={() => setPage('tool')} />
          </div>
        )}

        {/* ── Tool page ── */}
        {page === 'tool' && (
          <>
            <div className="hero">
              <span className="hero-icon">🔬</span>
              <h1 className="hero-title">{locale === 'zh' ? '虚拟药物筛选' : 'Virtual Drug Screening'}</h1>
              <p className="hero-sub">
                {locale === 'zh'
                  ? '上传基因表达特征，预测哪些化合物最有可能逆转您的细胞状态变化，无需编程基础。'
                  : 'Upload a gene-expression signature to prioritize compounds that may reverse the observed cellular state, with no coding required.'}
              </p>
            </div>

            <StepIndicator current={step} locale={locale} />

            {error && step === 'config' && (
              <ErrorBox message={error} onDismiss={() => setError(null)} />
            )}

            {step === 'upload' && (
              <UploadStep
                checkpoints={checkpoints}
                backendOnline={backendOnline}
                locale={locale}
                onReady={f => {
                  setFile(f)
                  setGeoConfig(null)
                  setGeoPreviewSelection(null)
                  setStep('config')
                }}
                onGeoReady={(config, preview) => {
                  setFile(null)
                  setGeoConfig(config)
                  setGeoPreviewSelection(preview)
                  setStep('config')
                }}
                onTutorial={() => setPage('tutorial')}
                onRefreshCheckpoints={refreshCheckpoints}
              />
            )}
            {step === 'config' && !error && (file || geoConfig) && (
              <ConfigStep
                file={file}
                geoConfig={geoConfig}
                geoPreview={geoPreviewSelection}
                onRun={handleRun}
                onBack={() => setStep('upload')}
                locale={locale}
              />
            )}
            {step === 'running' && <RunningStep />}
            {step === 'results' && results && (
              <EnhancedResultsStep response={results} onReset={handleReset} locale={locale} />
            )}
          </>
        )}
      </main>

      <footer className="footer">
        DrugReflector · Cellarity Inc. ·{' '}
        <a href="https://doi.org/10.5281/zenodo.16912444">模型文件（Zenodo）</a>
        {' · '}
        <button
          style={{ background: 'none', border: 'none', color: 'var(--mid)', cursor: 'pointer', fontSize: 'inherit', fontFamily: 'inherit' }}
          onClick={() => setPage('tutorial')}
        >
          使用教程
        </button>
      </footer>
    </div>
  )
}
