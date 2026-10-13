




from pathlib import Path
import argparse, csv, json, hashlib, math

base=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=base/'output')
args=parser.parse_args()
out=args.output.resolve()
def table(p):
    with p.open(newline='',encoding='utf-8-sig') as f: return list(csv.DictReader(f))
def near(a,b): return math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=1e-12)
checks=[]
def check(name, condition):
    assert condition, name
    checks.append({'check':name,'status':'PASS'})

meta=table(base/'input/sample_metadata.csv')
wide=table(base/'input/subtype_abundance.csv')
samples=table(out/'sample_summaries.csv')
check('72 unique sample columns and 2646 unique subtypes',len(meta)==72 and len(wide)==2646 and
      len({r['Sample'] for r in meta})==72 and len({r['Subtype'] for r in wide})==2646)
check('each input abundance is explicit finite and nonnegative',all(
      math.isfinite(float(row[m['Sample']])) and float(row[m['Sample']])>=0 for row in wide for m in meta))
for sample in samples:
    values=[float(row[sample['Sample']]) for row in wide]
    check(sample['Sample']+' richness and total reconcile independently',
          int(sample['Richness'])==sum(v>0 for v in values) and near(sample['Total_abundance'],math.fsum(values)))
external=table(base/'input/external_benchmark.csv')
merged=table(out/'panel_d_global/merged_sample_with_meta_cleaned.csv')
check('35 benchmark input rows, 33 external plus 18 assigned demo totals',len(external)==35 and len(merged)==51)
external_original={r['Sample_ID']:r['Copies_per_cell'] for r in external if r['COUNTRY']!='This Study'}
external_output={r['Sample_ID']:r['Copies_per_cell'] for r in merged if r['COUNTRY']!='This Study'}
check('all 33 external values preserved in analysis',set(external_original)==set(external_output) and all(
      near(external_original[k],external_output[k]) for k in external_original))
check('both original study benchmark examples replaced, not mixed with demo totals',
      not ({r['Sample_ID'] for r in external if r['COUNTRY']=='This Study'} & {r['Sample_ID'] for r in merged}))
study_output={r['Sample_ID']:r['Copies_per_cell'] for r in merged if r['COUNTRY']=='This Study'}
study_totals={r['Sample']:r['Total_abundance'] for r in samples if r['Setting']=='WWTP'}
check('18 canonical assigned WWTP totals reconcile',set(study_output)==set(study_totals) and
      all(near(study_output[k],study_totals[k]) for k in study_output))
for panel in ['panel_b_NB','panel_c_LMM']:
    pairs=table(out/panel/'six_setting_contrasts_BH.csv')
    check(panel+' exactly six finite contrasts and BH q >= p',len(pairs)==6 and all(
       math.isfinite(float(r['P'])) and float(r['q_BH'])+1e-14>=float(r['P']) for r in pairs))
g=table(out/'panels_gh/tables/g_marker_detection_fraction_by_setting.csv')
h=table(out/'panels_gh/tables/h_monthly_mean_abundance_by_setting.csv')
check('g complete 9 by 4 detection fractions',len(g)==36 and all(
      int(r['n_samples'])==18 and near(r['detection_fraction'],int(r['positive_n'])/18) for r in g))
check('h complete 8 by 4 by 3 monthly observations',len(h)==96 and all(int(r['n_samples'])==6 for r in h))
coverage=table(out/'panel_coverage.csv')
check('all required Figure 2 panels ran',all(r['Status']=='RUN' for r in coverage) and
      {'Fig2a','Fig2b','Fig2c','Fig2b/c plots','Fig2d','Fig2e','Fig2f','Fig2g/h'} <= {r['Panel'] for r in coverage})
check('figure PDF and PNG exports are nonempty',len(list(out.rglob('*.pdf')))>=12 and
      len(list(out.rglob('*.png')))>=10 and all(p.stat().st_size>1000 for p in out.rglob('*.pdf')))
checksums=table(base/'input_checksums.csv')
check('input_sha256',all(
      hashlib.sha256((base/'input'/row['path']).read_bytes()).hexdigest()==row['sha256']
      for row in checksums))
result={'status':'PASS','checks':checks,'input_cells':72*2646,'benchmark_analysis_rows':51,
        'pdf_files':len(list(out.rglob('*.pdf'))),'png_files':len(list(out.rglob('*.png')))}
(out/'validation_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='checks'}))
