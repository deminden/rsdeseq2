"""Generate tests/run_equivalence.rs: every public DeseqBuilder::fit* method
against its DeseqBuilder::run request (design note 4.2, step 1).

Usage: python3 tools/validation/gen_run_equivalence.py > tests/run_equivalence.rs && cargo fmt
"""
import sys
import glob, re, os
HERE = os.path.dirname(os.path.abspath(__file__))
BUILDER = os.path.join(HERE, "..", "..", "src", "core", "implementation", "builder")
HEAD = os.path.join(HERE, "run_equivalence_head.rs")
rows = []
for path in sorted(glob.glob(os.path.join(BUILDER, "**", "*.rs"), recursive=True)):
    text = open(path).read()
    for m in re.finditer(r"pub fn (fit[a-z0-9_]*)\s*(<[^(]*>)?\((.*?)\)\s*->\s*(.*?)\{", text, re.S):
        if m.group(1) == "fit_type":
            continue
        args = " ".join(m.group(3).split()).replace("&self,", "").strip().rstrip(",")
        rows.append((m.group(1), args, " ".join(m.group(4).split())))
out=[]
def arg_expr(a, formula):
    name,ty=[x.strip() for x in a.split(':',1)]
    if name=='counts': return '&f.counts'
    if name in('design','full_design'): return '&f.design'
    if name=='reduced_design': return '&f.reduced'
    if name=='coefficient': return '2'
    if name=='coefficient_name': return '"condition_B_vs_A"'
    if name in ('formula','full_formula'): return 'FULL'
    if name=='reduced_formula': return 'REDUCED'
    if name=='dispersions': return '&f.dispersions'
    if name=='replacement_options': return '&f.cooks'
    if name=='model_frame': return '&f.frame'
    if name=='sample_levels': return 'LEVELS'
    if name=='contrast':
        if 'ContrastSpec' in ty: return '&f.spec'
        if 'FactorLevel' in ty: return 'f.factor_level()'
        if 'ResultsContrast' in ty: return 'RC'
        return '&f.numeric'
    raise Exception(a)
kinds={'Result<(DeseqFit, DeseqResults), DeseqError>':'fr','Result<DeseqFit, DeseqError>':'fo',
 'Result<CooksReplacementWaldOutput, DeseqError>':'rw','Result<CooksReplacementLrtOutput, DeseqError>':'rl',
 'Result<CooksReplacementTestOutput, DeseqError>':'rt'}
n=0
for name,args,ret in rows:
    if 'fast_vst' in name: continue
    kind=kinds[ret]
    arglist=[a for a in args.split(', ') if a] if args else []
    # split respecting generics: args contain no commas inside except Option<&[S]> fine
    formula='formula' in name
    stage=any(s in name for s in ['size_factors','gene_wise','_trend_','map_dispersions'])
    lrt_named='lrt' in name
    leaf = ('glm_mu' in name) or ('linear_mu' in name) or ('fixed_dispersion' in name)
    has_levels=any(a.startswith('sample_levels') for a in arglist)
    rc_variants=[('character','f.character()'),('numeric','f.numeric_request()')] if any('ResultsContrast' in a for a in arglist) else [(None,None)]
    lv_variants=[('levels','Some(f.levels.as_slice())'),('nolevels','None::<&[String]>')] if has_levels else [(None,None)]
    builders=['wald']
    if not lrt_named and not leaf and not stage and 'with_results' not in name.replace('with_test_results','') :
        builders=['wald','lrt']
    if 'with_test_results' in name: builders=['wald','lrt']
    if name in ('fit','fit_formula','fit_formula_name','fit_contrast','fit_contrast_spec','fit_factor_level_contrast'): builders=['wald','lrt']
    if not lrt_named and not leaf and not stage and name.startswith('fit_with_results') and 'with_test' not in name and 'cooks' not in name:
        builders=['wald','lrt']
    for bname in builders:
      for rcn,rce in rc_variants:
        for lvn,lve in lv_variants:
            # analysis
            a=['Analysis::formula(&f.counts, FULL)' if formula else ('Analysis::without_design(&f.counts)' if name=='fit_size_factors_and_base_means' else 'Analysis::new(&f.counts, &f.design)')]
            if lrt_named or (bname=='lrt' and not leaf):
                if formula and lrt_named: a.append('.lrt_formula(REDUCED)')
                else: a.append('.lrt(&f.reduced)')
            if 'fixed_dispersion' in name: a.append('.fixed_dispersion(&f.dispersions)')
            if 'glm_mu' in name: a.append('.means(MeanModel::GlmMu)')
            if 'linear_mu' in name: a.append('.means(MeanModel::LinearMu)')
            if stage and 'glm_mu' not in name and 'linear_mu' not in name: pass
            if 'parametric' in name: a.append('.fit_type(FitType::Parametric)')
            if '_mean_dispersion_trend' in name: a.append('.fit_type(FitType::Mean)')
            if '_local_dispersion_trend' in name: a.append('.fit_type(FitType::Local)')
            if stage:
                st='SizeFactors' if 'size_factors' in name else 'GeneWise' if 'gene_wise' in name else 'Trend' if '_trend_' in name else 'Map'
                if name!='fit_size_factors_and_base_means': a.append(f'.stop_after(Stage::{st})')
            else:
                if any(a_.startswith('model_frame') for a_ in arglist):
                    a.append('.result(ResultColumn::Results { contrast: RC, sample_levels: None, model_frame: Some(&f.frame) })'.replace('RC','&'+rce.replace('f.','f.') if False else 'RCREF'))
                elif rce is not None:
                    a.append('.result(ResultColumn::Results { contrast: RCREF, sample_levels: LVREF, model_frame: None })')
                elif any('ContrastSpec' in a_ for a_ in arglist): a.append('.result(ResultColumn::Spec(&f.spec))')
                elif any('FactorLevel' in a_ for a_ in arglist): a.append('.result(ResultColumn::FactorLevel(f.factor_level()))')
                elif any(a_.startswith('contrast') for a_ in arglist): a.append('.result(ResultColumn::Numeric(&f.numeric))')
                elif any(a_.startswith('coefficient_name') for a_ in arglist): a.append('.result(ResultColumn::Name("condition_B_vs_A"))')
                elif any(a_.startswith('coefficient:') for a_ in arglist): a.append('.result(ResultColumn::Coefficient(2))')
            if 'cooks' in name: a.append('.cooks_replacement(&f.cooks)')
            analysis=''.join(a)
            call_args=[arg_expr(x,formula) for x in arglist]
            call=', '.join(call_args)
            tag=f'{name}/{bname}' + (f'/{rcn}' if rcn else '') + (f'/{lvn}' if lvn else '')
            lv_run = 'Some(f.levels.as_slice())' if lvn=='levels' else 'None'
            body=f'''    {{
        let b = f.builder("{bname}");
{'        let rc = '+rce+';' if rce else ''}
        let method = {kind}(b.{name}({call}));
        let run = run_{kind}(b.run(&{analysis}));
        check(&mut report, "{tag}", method, run);
    }}
'''
            body=body.replace('FULL','f.full_formula').replace('REDUCED','f.reduced_formula').replace('RCREF','&rc').replace('RC','&rc').replace('LVREF',lv_run).replace('LEVELS',lve if lve else 'None::<&[String]>')
            if 'FactorLevel' in body or 'factor_level' in body: pass
            out.append(body); n+=1
print(open(HEAD).read(), end='')
print('// @generated below by tools/validation/gen_run_equivalence.py from the builder')
print('// signatures: one block per method and builder/request variant.')
print(f'// {n} comparisons.')
print('fn compare_all(f: &Fixture) -> Vec<String> {\n    let mut report = Vec::new();')
print(''.join(out))
print('    report\n}')
