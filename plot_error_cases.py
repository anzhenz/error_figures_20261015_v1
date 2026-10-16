#!/usr/bin/env python3
"""Draw all frozen ECG errors from the exported wavepack, CPU only.

Usage: python -B plot_error_cases.py INPUT.zip OUTPUT_DIR [--sample]
All categories are verified matching/timing/stage descriptions, not inferred
medical diagnoses. FP nearest-reference beat types never label the FP itself.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import csv
import hashlib
import html
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import textwrap
import zipfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from PIL import Image

WAVE = "#455c6b"
GREEN = "#188452"
BLUE = "#2563cf"
RED = "#d22132"
SHADE = "#dcecf5"
CATEGORY_NAMES = {
    "FP_EARLY": "FP｜最近专家标注前的额外候选",
    "FP_LATE": "FP｜最近专家标注后的额外候选",
    "FP_DUP": "FP｜同一专家心搏附近的重复候选",
    "FN_WINDOW": "FN｜专家参考位置超出波形输入窗",
    "FN_MODEL": "FN｜模型阶段漏检，专家位置在窗口内",
    "FN_RR": "FN｜RR Rule 2删除真实心搏候选",
}


def digest_bytes(data):
    return hashlib.sha256(data).hexdigest()


def rows(data):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))


def load_pack(path):
    z = zipfile.ZipFile(path)
    if z.testzip() is not None:
        raise ValueError("ZIP CRC check failed")
    registry = json.loads(z.read("file_registry.json"))
    for entry in registry["files"]:
        name = entry["path"]
        if name.startswith("../"):
            continue
        if digest_bytes(z.read(name)) != entry["sha256"]:
            raise ValueError("Input hash mismatch: " + name)
    cases = rows(z.read("error_cases.csv"))
    if len(cases) != 177 or sum(c["case_type"]=="FP" for c in cases)!=97:
        raise ValueError("Unexpected case inventory")
    tables = {}
    for c in cases:
        key = (c["db_name"],c["record"])
        if key not in tables:
            prefix = f"tables/{key[0]}/{key[1]}/"
            ca, ev = rows(z.read(prefix+"all_ap_details.csv")), rows(z.read(prefix+"expert_events.csv"))
            tables[key] = (ca, ev)
    return z, cases, tables


def setup_font(output):
    asset = output/"assets/plot_font.ttf"
    asset.parent.mkdir(parents=True,exist_ok=True)
    if not asset.exists():
        import fitz
        asset.write_bytes(fitz.Font("cjk").buffer)
    fm.fontManager.addfont(str(asset))
    name = fm.FontProperties(fname=str(asset)).get_name()
    plt.rcParams.update({"font.family":[name,"DejaVu Sans"],"axes.unicode_minus":False,
                         "font.size":9,"axes.labelsize":8.6,"xtick.labelsize":8.4,
                         "ytick.labelsize":8,"axes.linewidth":0.65,
                         "savefig.facecolor":"white","svg.fonttype":"path"})


def classify(c, tables):
    ca, ev = tables[(c["db_name"],c["record"])]
    offset = float(c["expert_offset_ms"])
    if c["case_type"]=="FP":
        if abs(offset)<=150:
            e = ev[int(c["expert_event_id"])]
            partner = e["final_matched_ap_index"]
            if partner=="" or int(partner)==int(c["ap_index"]):
                raise ValueError("FP within tolerance lacks a competing final TP")
            category = "FP_DUP"
            explanation = (f"候选距最近专家标注{offset:+.1f}ms，仍在±150ms匹配容差内；"
                           f"该专家事件已由AP#{partner}匹配，当前候选形成额外FP。")
        else:
            category = "FP_EARLY" if offset>150 else "FP_LATE"
            direction = "前" if offset>0 else "后"
            explanation = (f"候选位于最近专家标注{direction}{abs(offset):.1f}ms，超出±150ms评价容差；"
                           "最终保留为FP。局部形态需结合波形观察，不能仅凭距离判为伪迹。")
    elif c["stage_route"]=="RR_ADDED_FN":
        category = "FN_RR"
        expert_id=int(c["expert_event_id"])
        expert_rr=(int(ev[expert_id]["expert_position"])-int(ev[expert_id-1]["expert_position"]))/2.0 if expert_id>0 else None
        explanation = (f"专家{c['expert_symbol']}类事件在模型阶段检出，随后关联候选被Rule 2删除。"
                       f"RR={float(c['rr_ms']):.1f}ms、参考RR={float(c['rr_ref_ms']):.1f}ms，"
                       f"比值{float(c['rr_ratio']):.3f}；专家位置相对AP为{offset:+.1f}ms，处于波形输入窗内。")
        if c["previous_final_keep_label"]=="FP" and expert_rr is not None:
            explanation+=(f"上一最终KEEP为FP（AP#{c['previous_final_keep_ap_index']}），"
                          f"而相邻专家心搏间隔为{expert_rr:.1f}ms；保留FP成为时间锚点，使当前RR显得过短。")
    elif offset>100 or offset < -150:
        category = "FN_WINDOW"
        explanation = (f"专家{c['expert_symbol']}类事件最终漏检；其参考位置在AP后{offset:.1f}ms，"
                       "超出波形输入的+100ms截止。该事实说明波形观察受限，不能单独证明模型无法识别。")
    else:
        category = "FN_MODEL"
        explanation = (f"专家{c['expert_symbol']}类事件在V11_B阶段检出，模型阶段起漏检；"
                       f"参考位置相对AP为{offset:+.1f}ms，在−150～+100ms波形窗内。"
                       "关联候选的低模型概率表明窗口内真实心搏未被放行。")
    if c["model_evaluated"]=="0":
        explanation += "该AP属于前30个直接放行候选，p_real=1并非模型输出。"
    if c["rr_input_exceeds_100ms"]=="1":
        explanation += "本样本RR输入引用了+100ms之外的位置，不能作为严格实时性能证据。"
    return category, explanation


def minmax_decimate(x,y,max_points=3600):
    if len(x)<=max_points:
        return x,y
    width=int(np.ceil(len(x)/(max_points//2)))
    keep=[]
    for start in range(0,len(y),width):
        end=min(len(y),start+width)
        lo=start+int(np.argmin(y[start:end])); hi=start+int(np.argmax(y[start:end]))
        keep.extend(sorted(set((lo,hi))))
    return x[keep],y[keep]


def write_atomic(path, data):
    fd, temp_name=tempfile.mkstemp(prefix=path.name+".",suffix=".tmp",dir=path.parent)
    try:
        view=memoryview(data)
        while view:
            written=os.write(fd,view)
            if written<=0:
                raise OSError("Incomplete artifact write")
            view=view[written:]
        os.fsync(fd)
        os.close(fd);fd=-1
        os.replace(temp_name,path)
        if digest_bytes(path.read_bytes()) != digest_bytes(data):
            raise OSError("Artifact verification failed: "+str(path))
    finally:
        if fd>=0:
            os.close(fd)
        if os.path.exists(temp_name):
            os.unlink(temp_name)


def make_figure(c,z,tables,out,svg=False):
    ca, ev = tables[(c["db_name"],c["record"])]
    ap0=int(c["ap_position"])
    aps=np.asarray([int(r["ap_position"]) for r in ca])
    gt=np.asarray([int(e["expert_position"]) for e in ev])
    with np.load(io.BytesIO(z.read(c["waveform_file"])),allow_pickle=False) as w:
        ecg=w["ecg"].copy(); es=int(w["ecg_start"])
        mean=w["mean"].copy(); ms=int(w["mean_start"])
    category, explanation=classify(c,tables)
    fig=plt.figure(figsize=(11.2,9.1))
    grid=fig.add_gridspec(4,1,left=.071,right=.98,top=.815,bottom=.183,
                        hspace=.63,height_ratios=[1.12,1.0,1.03,.92])
    axes=[fig.add_subplot(grid[i]) for i in range(4)]
    specs=[(ecg,es,-.5,.5,"a  ECG_FILT：局部波形 ±0.5 s"),
           (mean,ms,-.5,.5,"b  MEAN_TEST：局部特征 ±0.5 s"),
           (ecg,es,-3,3,"c  ECG_FILT：节律上下文 ±3 s（离线复核）"),
           (ecg,es,-20,0,"d  ECG_FILT：前20秒历史背景（离线复核）")]
    for panel,(ax,(signal,start,left,right,title)) in enumerate(zip(axes,specs)):
        sample_left=max(start,int(round(ap0+left*2000)))
        sample_right=min(start+len(signal)-1,int(round(ap0+right*2000)))
        if sample_right<sample_left:
            raise ValueError("Missing panel data: "+c["case_id"])
        pos=np.arange(sample_left,sample_right+1,dtype=np.int64)
        xx=(pos-ap0)/2000
        yy=signal[pos-start]
        xp,yp=minmax_decimate(xx,yy,3600 if panel>=2 else 8000)
        ax.plot(xp,yp,color=WAVE,lw=.77 if panel<3 else .58,zorder=2)
        lo,hi=float(yy.min()),float(yy.max())
        span=max(hi-lo,1.0)
        ax.set_ylim(lo-.14*span,hi+.20*span)
        ax.set_xlim(left,right)
        if panel<3:
            ax.axvspan(-.15,.1,color=SHADE,alpha=.82,zorder=0)
            ax.axvline(.1,color="#7f91a1",lw=.8,ls="-.",zorder=1)
        # Original APs are shown regardless of the model/RR keep decision.
        aa=aps[(aps>=sample_left)&(aps<=sample_right)]
        gg=gt[(gt>=sample_left)&(gt<=sample_right)]
        for p in aa:
            ax.axvline((p-ap0)/2000,color=GREEN,ls=":",lw=.65,alpha=.40 if panel>=2 else .65,zorder=1)
        for p in gg:
            ax.axvline((p-ap0)/2000,color=BLUE,ls="--",lw=.68,alpha=.36 if panel>=2 else .72,zorder=1)
        marker_size=3.0 if panel>=2 else 5.0
        ax.plot((aa-ap0)/2000,signal[aa-start],ls="none",marker="v",ms=marker_size,
                color=GREEN,mec="white",mew=.22,zorder=4,clip_on=True)
        ax.plot((gg-ap0)/2000,signal[gg-start],ls="none",marker="^",ms=marker_size+.6,
                color=BLUE,mec="white",mew=.22,zorder=5,clip_on=True)
        if start<=ap0<start+len(signal):
            ax.plot(0,float(signal[ap0-start]),marker="*",ms=12 if panel<2 else 10,
                    mfc=RED,mec="white",mew=.5,ls="none",zorder=8,clip_on=False)
        ax.set_title(title,loc="left",fontsize=9.4,pad=4)
        ax.set_ylabel("幅值 / AD")
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
        ax.grid(axis="y",color="#e4e8eb",lw=.5,zorder=0)
        ax.tick_params(length=3,pad=2)
        for edge in ("top","right"):
            ax.spines[edge].set_visible(False)
        for edge in ("left","bottom"):
            ax.spines[edge].set_color("#a7b0b7")
        ax.set_xlabel("相对当前AP的时间 / s",labelpad=2)
        if panel<2:
            ax.set_xticks([-.5,-.25,0,.1,.25,.5])
        if panel==0:
            ax.text(.1,.98,"+100ms截止",transform=ax.get_xaxis_transform(),fontsize=7.8,
                    color="#61778a",ha="left",va="top")
            target=int(c["expert_position"])
            if sample_left<=target<=sample_right:
                ax.annotate(c["expert_symbol"],((target-ap0)/2000,float(signal[target-start])),
                            xytext=(3,8),textcoords="offset points",fontsize=8,color=BLUE,
                            bbox={"facecolor":"white","alpha":.78,"edgecolor":"none","pad":.5})
    label=CATEGORY_NAMES[category]
    fig.text(.071,.976,f"{c['db_name'].upper()} / {c['record']}   ·   AP #{int(c['ap_index']):05d}   ·   {label}",
             fontsize=12.0,fontweight="bold",ha="left",va="top",color="#152b3d")
    probability=(f"p_real={float(c['p_real']):.4f}  |  p_noise={float(c['p_noise']):.4f}"
                 if c["model_evaluated"]=="1" else "前30个AP直接放行：p_real=1（非模型输出）")
    expert_role="漏检专家事件" if c["case_type"]=="FN" else "最近专家参考"
    fig.text(.071,.943,f"{probability}   |   {expert_role} {c['expert_symbol']}：{float(c['expert_offset_ms']):+.1f} ms   |   AP={ap0}",
             fontsize=9.1,ha="left",va="top",color="#31485a")
    reason="Rule 2删除" if c["delete_rule"]=="RULE2_COOLDOWN_NOISE" else "Rule 1删除" if c["delete_rule"] else "模型删除" if c["v707_keep"]=="0" else "最终保留"
    def number(k,places=1):
        try:
            value=float(c[k]); return f"{value:.{places}f}" if np.isfinite(value) else "未建立"
        except (ValueError,TypeError):
            return "未建立"
    if c["v707_keep"]=="0":
        rrline="模型已删除；RR后处理未评估此候选，表中空RR字段不代表节律不稳定"
    else:
        rrline=(f"RR={number('rr_ms')} ms / 参考={number('rr_ref_ms')} ms / 比值={number('rr_ratio',3)}"
                f"   |   稳定={c['history_stable']}，保护计数={c['cooldown_before']}   |   {reason}")
    rrline+=f"   |   模型RR输入位置超过+100ms：{'是' if c['rr_input_exceeds_100ms']=='1' else '否'}"
    fig.text(.071,.916,rrline,fontsize=8.6,ha="left",va="top",color="#526677")
    handles=[Line2D([],[],color=RED,marker="*",ms=10,ls="none",label="当前错误候选"),
             Line2D([],[],color=GREEN,marker="v",ms=5,ls=":",lw=.8,label="动态阈值候选AP"),
             Line2D([],[],color=BLUE,marker="^",ms=5,ls="--",lw=.8,label="专家标注"),
             Patch(facecolor=SHADE,label="波形输入窗 −150～+100ms")]
    fig.legend(handles=handles,loc="upper left",bbox_to_anchor=(.063,.887),ncol=4,
               frameon=False,fontsize=8.9,columnspacing=2.0,handlelength=2.0)
    # A short factual caption is included in every image, with full text in CSV.
    if category=="FP_DUP":
        short="同一专家心搏已有其他候选匹配，当前候选仍被保留，形成重复FP。"
    elif category.startswith("FP"):
        side="前" if float(c["expert_offset_ms"])>0 else "后"
        short=f"最终保留的候选位于最近专家标注{side}{abs(float(c['expert_offset_ms'])):.1f}ms，超出±150ms匹配容差。"
    elif category=="FN_RR":
        j=int(c["expert_event_id"])
        gt_rr=(int(ev[j]["expert_position"])-int(ev[j-1]["expert_position"]))/2.0
        short=(f"上一KEEP为{c['previous_final_keep_label']}，算法RR仅{float(c['rr_ms']):.1f}ms；"
               f"专家相邻RR为{gt_rr:.1f}ms。Rule 2删除了真实心搏候选。")
    elif category=="FN_WINDOW":
        short=f"专家参考位置在AP后{float(c['expert_offset_ms']):.1f}ms，超过波形+100ms截止；事件最终漏检。"
    else:
        short=f"专家{c['expert_symbol']}类心搏参考位置在波形输入窗内，但模型阶段起漏检。"
    if c["model_evaluated"]=="0":
        short+="该点为起始直接放行，不能解释为模型高置信。"
    short="\n".join(textwrap.wrap("案例说明："+short,width=68,break_long_words=False))
    fig.text(.071,.102,short,fontsize=9.5,ha="left",va="top",color="#243c4e",linespacing=1.25)
    fig.text(.071,.063,"说明：窗口外及历史波形仅供人工复核；模型历史为前30个原始AP，20秒视图仅作背景。",
             fontsize=8.1,ha="left",va="top",color="#657685")
    fig.text(.071,.034,c["case_id"]+"   |   当前冻结权重：离线错误分析，严格实时性审计未通过",
             fontsize=7.8,ha="left",va="top",color="#7a8792")
    out.parent.mkdir(parents=True,exist_ok=True)
    png_buffer=io.BytesIO()
    fig.savefig(png_buffer,dpi=180,format="png")
    payload=png_buffer.getvalue()
    with Image.open(io.BytesIO(payload)) as picture:
        picture.verify()
    write_atomic(out,payload)
    if svg:
        svg_buffer=io.BytesIO()
        fig.savefig(svg_buffer,format="svg")
        write_atomic(out.with_suffix(".svg"),svg_buffer.getvalue())
    plt.close(fig)
    return category, explanation


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("input",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--sample",action="store_true")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    setup_font(args.output)
    z,cases,tables=load_pack(args.input)
    chosen=cases
    if args.sample:
        wanted={"MITDB_105_FP_ap00455","MITDB_228_FP_ap00329",
                "MITDB_203_FN_event01942_ap01984","MITDB_105_FN_event01705_ap01717",
                "MITDB_105_FN_event02447_ap02497","MITDB_105_FN_event01870_ap01905",
                "MITDB_228_FP_ap00010"}
        chosen=[c for c in cases if c["case_id"] in wanted]
    results=[]
    for i,c in enumerate(chosen,1):
        path=args.output/"images"/c["case_type"]/c["record"]/(c["case_id"]+".png")
        cat,text=make_figure(c,z,tables,path)
        results.append({**c,"error_category":cat,"error_category_zh":CATEGORY_NAMES[cat],
                        "explanation":text,"image":path.relative_to(args.output).as_posix()})
        if i%15==0 or i==len(chosen):
            print(f"Rendered {i}/{len(chosen)}",flush=True)
    with (args.output/"case_descriptions.csv").open("w",encoding="utf-8-sig",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    z.close()
    if not args.sample:
        from finalize_error_gallery import finish
        finish(args.input,args.output)
    print("DONE",args.output,flush=True)


if __name__=="__main__":
    main()
