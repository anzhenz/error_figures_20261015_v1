#!/usr/bin/env python3
"""Build the offline image index, verified statistics, report examples and registry."""
import sys
sys.dont_write_bytecode=True
import argparse
from collections import Counter
import csv
import hashlib
import html
import io
import json
from pathlib import Path
import shutil
from PIL import Image
import plot_error_cases as plot

EXAMPLES=[
 ("01_QRS_like_high_confidence_FP","MITDB_105_FP_ap00455","QRS样局部突变被高置信放行",
  "错误AP附近可见先负后正的快速局部变化，MEAN_TEST也有响应；模型p_real=0.9141。最近专家参考位置在AP后199.5ms，超过评价容差，当前候选仍被保留为FP。图像提示局部形态具有迷惑性，但不能仅凭本图确定噪声的生理来源。"),
 ("02_visible_V_beat_FN","MITDB_203_FN_event01942_ap01984","窗口内V类真实心搏被模型删除",
  "专家V类心搏位于AP后16ms，已处于波形输入窗内。局部ECG及MEAN_TEST存在多个峰谷，模型p_real=0.2832，关联候选被删除。该例说明异常形态心搏即使参考位置可见，也可能被模型判负。"),
 ("03_duplicate_trigger_FP","MITDB_228_FP_ap00329","同一专家心搏附近的重复候选",
  "专家V类参考位置在当前AP前86.5ms；同一专家事件已由另一候选匹配。当前候选仍以p_real=0.7227保留，形成额外FP。该例属于匹配关系可直接核实的重复候选。"),
 ("04_reference_after_deadline_FN","MITDB_105_FN_event01870_ap01905","专家参考位置超过波形截止的FN",
  "专家N类参考位置位于AP后119.5ms，晚于波形输入的+100ms截止；模型p_real=0.1436，事件最终漏检。这说明参考位置尚未进入当前波形窗，但不能据此断言此前完全没有可辨识线索。"),
 ("05_RR_false_anchor_FN_1717","MITDB_105_FN_event01705_ap01717","保留FP成为RR时间锚点，继而删除真实心搏",
  "当前N类心搏先通过模型，随后被Rule 2删除。上一最终KEEP是FP，算法RR仅201ms，而专家相邻心搏RR为628ms；错误时间锚点使真实心搏看起来明显提前。"),
 ("06_initial_forced_keep_FP","MITDB_228_FP_ap00010","前30个AP直接放行造成的FP",
  "当前候选属于前30个原始AP，未进入模型，p_real=1是测试程序设定的直接放行值。它最终为FP，最近专家位置在AP后185.5ms。该例不能列为模型高置信误判。"),
 ("07_RR_false_anchor_FN_2497","MITDB_105_FN_event02447_ap02497","第二例RR时间锚点污染导致的FN",
  "另一N类心搏同样先通过模型再被Rule 2删除。上一最终KEEP仍是FP，算法RR为266.5ms，而专家相邻RR为617ms；两例共同提示应分析保留FP如何影响后续RR状态。"),
]


def csv_bytes(rows):
    s=io.StringIO(newline="")
    w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return s.getvalue().encode("utf-8-sig")


def finish(input_path,output):
    plot.setup_font(output)
    z,cases,tables=plot.load_pack(input_path)
    results=[]
    for c in cases:
        cat,text=plot.classify(c,tables)
        rel=Path("images")/c["case_type"]/c["record"]/(c["case_id"]+".png")
        image=output/rel
        with Image.open(image) as im:
            im.verify()
        # Regenerate the two RR explanations with the verified false-anchor evidence.
        if cat=="FN_RR":
            plot.make_figure(c,z,tables,image)
        results.append({**c,"error_category":cat,"error_category_zh":plot.CATEGORY_NAMES[cat],
                        "explanation":text,"image":rel.as_posix(),
                        "rr_postprocessor_evaluated":"是" if c["v707_keep"]=="1" else "否：模型已删除"})
    plot.write_atomic(output/"case_descriptions.csv",csv_bytes(results))
    counts=Counter(r["error_category"] for r in results)
    distribution=[]
    for cat in ("FP_EARLY","FP_LATE","FP_DUP","FN_WINDOW","FN_MODEL","FN_RR"):
        total=97 if cat.startswith("FP") else 80
        distribution.append({"error_category":cat,"类别":plot.CATEGORY_NAMES[cat],"数量":counts[cat],
                             "分母":total,"组内占比":f"{counts[cat]/total*100:.2f}%",
                             "全部177例占比":f"{counts[cat]/177*100:.2f}%"})
    plot.write_atomic(output/"error_category_counts.csv",csv_bytes(distribution))
    record_counts=plot.rows(z.read("error_counts_by_record.csv"))
    plot.write_atomic(output/"error_counts_by_record.csv",csv_bytes(record_counts))
    plot.write_atomic(output/"source_summary.json",z.read("source_summary.json"))
    plot.write_atomic(output/"plot_standard.json",z.read("plot_standard.json"))
    summary=json.loads(z.read("source_summary.json"))
    high_conf=sum(r["case_type"]=="FP" and r["model_evaluated"]=="1" and float(r["p_real"])>=.9 for r in results)
    stats={"cases":177,"final_fp":97,"added_fn":80,"model_stage_added_fn":78,"rr_added_fn":2,
           "forced_keep_fp":sum(r["case_type"]=="FP" and r["model_evaluated"]=="0" for r in results),
           "model_output_fp_p_real_ge_0_90":high_conf,
           "fn_expert_position_after_100ms":35,"fn_expert_position_inside_waveform":45,
           "fn_with_rr_input_index_after_100ms":sum(r["case_type"]=="FN" and r["rr_input_exceeds_100ms"]=="1" for r in results),
           "categories":distribution,"source_zip_sha256":plot.digest_bytes(input_path.read_bytes()),
           "classification_basis":"Verified matching, reference timing and deletion stage; not all morphology causes are established",
           "all_unique_case_pngs_verified":True,"model_or_gpu_inference_performed":False}
    plot.write_atomic(output/"analysis_summary.json",json.dumps(stats,ensure_ascii=False,indent=2).encode())
    example_dir=output/"report_examples";example_dir.mkdir(exist_ok=True)
    by_id={r["case_id"]:r for r in results}
    example_rows=[]
    for name,case_id,title,caption in EXAMPLES:
        c=by_id[case_id]
        image=example_dir/(name+".png")
        plot.make_figure(c,z,tables,image,svg=True)
        example_rows.append({"figure":name,"case_id":case_id,"title":title,"caption":caption,
                             "png":image.relative_to(output).as_posix(),
                             "svg":image.with_suffix('.svg').relative_to(output).as_posix()})
    plot.write_atomic(output/"representative_cases.csv",csv_bytes(example_rows))
    captions="# 报告代表案例\n\n前三张对应：QRS样局部变化FP、窗口内异常心搏FN、重复候选FP。其余为补充。\n\n"
    for r in example_rows:
        captions+=f"## {r['title']}\n\n案例：`{r['case_id']}`\n\n{r['caption']}\n\n![{r['title']}]({r['png']})\n\n"
    plot.write_atomic(output/"report_captions.md",captions.encode())
    # Self-contained, local index; all image/data links remain inside this folder.
    cards=[]
    for r in results:
        text=html.escape(r["explanation"])
        cards.append(f'''<article data-kind="{r['case_type']}" data-record="{r['record']}" data-category="{r['error_category']}">
<h3>{html.escape(r['case_id'])}</h3><div class="tag">{html.escape(r['error_category_zh'])}</div>
<a href="{r['image']}" target="_blank"><img loading="lazy" src="{r['image']}" alt="{html.escape(r['case_id'])}"></a>
<p>{text}</p></article>''')
    options=''.join(f'<option value="{x}">{x}</option>' for x in sorted({r['record'] for r in results}))
    cat_options=''.join(f'<option value="{c}">{html.escape(plot.CATEGORY_NAMES[c])}</option>' for c in counts)
    table=''.join(f"<tr><td>{html.escape(r['类别'])}</td><td>{r['数量']}</td><td>{r['组内占比']}</td></tr>" for r in distribution)
    page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>ECG 全部177例错误波形图</title><style>
*{box-sizing:border-box}body{margin:0;background:#f3f6f8;color:#183244;font:15px/1.6 system-ui,"Microsoft YaHei",sans-serif}
header,main{max-width:1500px;margin:auto;padding:24px}header{background:white;border-bottom:1px solid #d9e1e7}h1{margin:0;font-size:28px}
.meta{color:#627787}.controls{display:flex;gap:10px;flex-wrap:wrap;margin:18px 0}input,select{padding:9px;border:1px solid #c4d0da;border-radius:5px;background:white}
input{min-width:280px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(440px,1fr));gap:20px}article{background:white;border:1px solid #dae3e9;border-radius:8px;padding:16px}
h3{font-size:16px;margin:0 0 5px}.tag{font-size:13px;color:#526a7d}img{width:100%;height:auto;display:block}p{margin:10px 0 0;font-size:14px}
details{background:white;padding:12px 16px;margin-bottom:18px}table{border-collapse:collapse}td,th{padding:6px 15px;border-bottom:1px solid #e0e7ec;text-align:left}
a{color:#1e62a1}article[hidden]{display:none}@media(max-width:500px){.grid{grid-template-columns:1fr}header,main{padding:15px}}
</style><header><h1>ECG 全部177例错误波形图</h1><div class="meta">当前冻结权重 · TEST27 · 97个最终FP / 80个新增FN（含2个RR新增FN）</div>
<p>红色五角星：当前错误候选；绿色倒三角与点线：原始AP；蓝色上三角与虚线：专家标注。浅色区域：−150～+100ms波形输入窗。点击图片查看原图。</p>
<p><a href="README.md">阅读口径说明</a> · <a href="case_descriptions.csv">全部逐例说明CSV</a> · <a href="report_captions.md">7个报告代表案例</a> · <a href="error_category_counts.csv">分类统计CSV</a></p></header><main>
<details><summary>可核实分类与比例（按FP/FN各自为分母）</summary><table><tr><th>类别</th><th>数量</th><th>组内占比</th></tr>TABLE</table>
<p>这些分类描述匹配关系、专家位置与判决阶段，并不把所有FP直接诊断为伪迹。专家蓝标不要求在正向波峰顶点。</p></details>
<div class="controls"><select id="kind"><option value="">全部 FP / FN</option><option>FP</option><option>FN</option></select><select id="record"><option value="">全部记录</option>RECORDS</select>
<select id="category"><option value="">全部类别</option>CATEGORIES</select><input id="query" placeholder="搜索AP编号、专家类型、规则或说明"><span id="count"></span></div>
<section class="grid">CARDS</section></main><script>
const articles=[...document.querySelectorAll('article')];function filter(){const k=document.getElementById('kind').value,r=document.getElementById('record').value,c=document.getElementById('category').value,q=document.getElementById('query').value.toLowerCase();let n=0;
for(const a of articles){const show=(!k||a.dataset.kind===k)&&(!r||a.dataset.record===r)&&(!c||a.dataset.category===c)&&(!q||a.textContent.toLowerCase().includes(q));a.hidden=!show;if(show)n++;}document.getElementById('count').textContent=`显示 ${n} / 177 例`;}
for(const id of ['kind','record','category','query'])document.getElementById(id).addEventListener('input',filter);filter();</script></html>'''
    page=page.replace('TABLE',table).replace('RECORDS',options).replace('CATEGORIES',cat_options).replace('CARDS','\n'.join(cards))
    plot.write_atomic(output/"index.html",page.encode())
    table_md="\n".join(f"| {r['类别']} | {r['数量']} | {r['组内占比']} |" for r in distribution)
    weight=summary["input_sha256_before"]["runs/V7.0.7/checkpoints/best_model_v7_0_7_best.pth"]
    readme=f'''# ECG当前冻结权重：完整错误图集

双击 `index.html` 可筛选查看全部177例，或直接打开 `images/`。每例均有四个视图和图内简要说明；完整字段与说明在 `case_descriptions.csv`。

## 数量与当前评价结果

- 最终FP：97个；相对V11_B新增FN：80个，恢复事件0个，因此净新增也是80。
- 模型阶段新增FN 78个；RR阶段再新增FN 2个。
- 最终 TP=138410、FP=97、FN=428；F1={summary['totals']['final']['f1']*100:.4f}%。
- RR共删除25个候选，FP由120降到97，同时FN由426增至428。
- 旧报告99 FP / 75净新增FN对应历史模型；据原作者说明权重曾覆盖。本图集使用当前已冻结并登记哈希的模型，不能称为旧权重精确复现。
- 当前权重SHA256：`{weight}`。

## 分类数量

以下为可复核的匹配、时序、判决阶段分类；形态与生理来源不由这些条件直接确定。

| 类别 | 数量 | FP组或FN组内占比 |
|---|---:|---:|
{table_md}

80个新增FN中，35个专家参考位置在波形+100ms截止之后，45个在波形窗内；窗内45个包括43个模型阶段漏检和2个RR新增漏检。
35个窗外参考位置不等于35个必然不可判样本：波形早期线索与RR分支可能仍有信息。

## 图例与读取方式

红色五角星=当前错误候选AP；绿色倒三角及绿色点线=全部原始动态阈值AP；蓝色上三角及蓝色虚线=专家标注位置。
浅色区域仅表示波形输入窗−150～+100ms，实际右端最后一个采样点为+99.5ms。
四个视图统一为局部ECG±0.5s、局部MEAN_TEST±0.5s、ECG节律上下文±3s、ECG历史背景−20～0s。
20秒历史仅用于人工复核，模型实际使用前30个原始AP的窗口和RR序列；历史与节律长图用分箱最小/最大值保留极值以改善渲染，局部两图保留全部采样点。
专家参考位置不一定落在正向峰顶。FP的专家类型字段属于“最近专家参考”，不能把FP本身标为N类或V类心搏。

## 两点需要特别区分

1. 2个FP来自前30个AP直接放行，p_real=1不是模型输出；其余真实模型输出中，p_real≥0.90的FP有{high_conf}个。0.90仅用于描述置信度，不是新判决阈值。
2. 被模型删除的候选不会进入RR规则评估；原始CSV的空RR字段及history_stable=0为初始化/未评估值，不能据此判断当时节律不稳定。

两例RR新增FN的前一个最终KEEP均为FP：AP#1717的算法RR为201ms、专家相邻RR为628ms；AP#2497分别为266.5ms与617ms。
这条状态链表明保留FP成为RR时间锚点，使后续真实心搏看起来提前。两例当前专家参考点均在波形窗内。

## 实时性与本轮范围

TEST27共有81个模型样本的RR输入引用了AP后超过100ms的位置；本图集80个新增FN中有28个带此标记，97个FP中为0。
本轮是当前冻结模型的离线错误画像，不能作为严格+100ms实时性能证明；上游滤波及特征生成的可用时延也尚未核实。
没有在TEST27上搜索阈值、改变权重或重新训练；本地绘图没有重新推理。

## 文件

- `images/FP/`：97个FP原图；`images/FN/`：80个新增FN原图。
- `report_examples/`：7个代表案例的PNG副本与SVG；前三个对应报告常用三类示例。
- `report_captions.md`：代表案例的简要描述。
- `case_descriptions.csv`：177例完整说明、概率、专家位置、RR状态、规则与文件路径。
- `error_category_counts.csv` / `error_counts_by_record.csv`：分类及记录统计。
- `plot_error_cases.py` / `finalize_error_gallery.py`：绘图与图集代码；`assets/plot_font.ttf`为已验证可显示中文的字体。
- `file_registry.json`：所有交付文件的SHA256登记（登记文件自身不递归计算）。

复现：安装numpy、matplotlib、pillow；首次缺少字体时还需pymupdf。保留assets字体可免去字体提取。
运行 `python -B plot_error_cases.py 原始波形包.zip 新输出目录` 会生成177张PNG，随后用同目录 `finalize_error_gallery.py` 完成索引、统计及代表图。
'''
    plot.write_atomic(output/"README.md",readme.encode())
    here=Path(__file__).resolve().parent
    for name in ("plot_error_cases.py","finalize_error_gallery.py"):
        if (here/name).resolve()!=(output/name).resolve():
            shutil.copyfile(here/name,output/name)
    # Final independent validation covers every case and the seven duplicate report figures.
    for p in output.rglob("*.png"):
        with Image.open(p) as im: im.verify()
    files=[]
    for p in sorted(output.rglob("*")):
        if p.is_file() and p.name!="file_registry.json":
            files.append({"path":p.relative_to(output).as_posix(),"bytes":p.stat().st_size,
                          "sha256":plot.digest_bytes(p.read_bytes())})
    plot.write_atomic(output/"file_registry.json",json.dumps({"status":"COMPLETE_177_CASES",
        "source_zip_sha256":stats["source_zip_sha256"],"files":files,
        "registry_self":{"path":"file_registry.json","sha256":None,"note":"Self-hash omitted to avoid recursion"}},ensure_ascii=False,indent=2).encode())
    z.close()
    print("COMPLETE",stats,flush=True)
    return stats


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("input",type=Path);p.add_argument("output",type=Path)
    args=p.parse_args();finish(args.input,args.output)
