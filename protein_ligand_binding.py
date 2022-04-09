""" 
MIT License

Copyright (c) 2021-2022 Wen Jiang

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

def import_with_auto_install(packages, scope=locals()):
    if isinstance(packages, str): packages=[packages]
    for package in packages:
        if package.find(":")!=-1:
            package_import_name, package_pip_name = package.split(":")
        else:
            package_import_name, package_pip_name = package, package
        try:
            scope[package_import_name] = __import__(package_import_name)
        except ImportError:
            import subprocess
            subprocess.call(f'pip install {package_pip_name}', shell=True)
            scope[package_import_name] =  __import__(package_import_name)
required_packages = "streamlit numpy pandas bokeh".split()
import_with_auto_install(required_packages)

import streamlit as st
import numpy as np
from bokeh.plotting import figure
from bokeh.models import LegendItem

def main():
    title = "Protein-Ligand Binding"
    st.set_page_config(page_title=title, layout="wide")

    with st.sidebar:    # sidebar at the left of the screen
        st.title(title)
        st.write(r'For a reversible binding reaction:$\\P + L\  {\rightleftharpoons}\  P{\cdot}L$')
        st.write(r"the binding affinity $K_D$ is defined as  $K_D = \frac{[P]_{free}[L]_{free}}{[P {\cdot}L]}$")
        st.write(r'then, the fraction of protein in ligand-bound state is:')
        st.latex(r'f = \frac{[P{\cdot}L]}{[P{\cdot}L]+[P]_{free}} = \frac{1}{1 + \frac{K_D}{[L]_{free}}}')
        st.latex(r'f = 0.5 \  when \  [L]_{free} = K_D')
        st.write('Note that $[L]_{free}$ means free, unbound ligand concentration that we cannot control directly. Instead, we can only control the total concentration of the protein ($P_T=[P{\cdot}L]+[P]_{free}$) and ligand ($L_T=[P{\cdot}L]+[L]_{free}$). So it is useful to express $f$ in terms of $P_T$ and $L_T$:')
        st.latex(r'f = \frac{1}{1 + \frac{K_D}{[L]_{free}}} = \frac{1}{1 + \frac{K_D}{[L_T]-[P_T]f}}')
        st.latex(r'P_Tf^2-(P_T+L_T+K_D)f+L_T=0')
        st.latex(r'f=\frac{P_T+L_T+K_D-\sqrt{(P_T+L_T+K_D)^2-4P_TL_T}}{2P_T}')
        st.write(r'So the fraction of protein with ligand bound is a function of total protein concentration, total ligand concentration, and $K_D$, instead of just the molar ratio of protein and ligand ($P_T:L_T$)')
        
        # make radio display horizontal
        st.markdown('<style>div.row-widget.stRadio > div{flex-direction:row;}</style>', unsafe_allow_html=True)

        hide_streamlit_style = """
        <style>
        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        </style>
        """
        st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

    # right-side main panel
    col1, col2 = st.columns((1, 4.5))
    with col1:
        mode = st.radio(label="Plotting Mode", options=["[P]", "[P]:[L]", "Kd"], index=0, help="Choose a plotting mode: Mode [P] will plot ligand concentration [L] vs fraction curves at different protein concentrations [P]; Mode [P]:[L] will plot protein concentration [P] vs fraction curves at different molar ratios of protein and ligand [P]:[L]; Mode Kd will plot ligand concentration [L] vs fraction curves at different binding affinities Kd")
        m = st.number_input('Protein mass per binding site (kDa)', value=100.0, min_value=1.0, step=1.0, format="%g", help="Same as the total protein mass for monomer protein. If the protein is an oligomer with each subunit having a binding site, the input value here should be the total protein mass divided by the number of subunits")
        if mode == "[P]":
            n = int(st.number_input('# of [P]', value=3, min_value=1, step=1, help="Number of protein concentrations to plot"))
            ps = [0] * n
            for i in range(n):
                c = 0.01 * np.power(10, i)
                ps[i] = st.number_input(f'[P] (mg/ml) - Curve {i+1}', value=c, min_value=0.0, step=0.1, format="%g", help=f"Total protein concentration for curve {i+1}")
            ps = np.sort(np.unique([ p for p in ps if p>0 ]))
            kd = st.number_input('Kd (nM)', value=1000.0, min_value=0.0, step=1.0, format="%g", help="Binding affinity")
            kd_molar = kd * 1e-9    # nM -> M
            p_molar_min = ps.min()/m * 1e3 # mg/ml -> μM
            p_molar_max = ps.max()/m * 1e3 # mg/ml -> μM
            lmin = st.number_input(r'[L] (μM) - min', value=p_molar_min*1e-2, min_value=0.0, step=p_molar_min*1e-2, format="%g", help="Minimal total ligand concentration to plot")
            lmax = st.number_input('[L] (μM) - max', value=p_molar_max*1e2, min_value=lmin*1e4, step=p_molar_max, format="%g", help="Maximal total ligand concentration to plot")
        elif mode == "[P]:[L]":
            n = int(st.number_input('# of [P]:[L]', value=3, min_value=1, step=1, help="Number of protein:ligand molar ratios to plot"))
            ratios = [1.0] * n
            for i in range(n):
                ratio = 1.0 * np.power(10, i)
                ratios[i] = st.number_input(f'[P]:[L]=1:? - Curve {i+1}', value=ratio, min_value=0.0, step=1.0, format="%g", help=f"Molar ratio of protein and ligand for curve {i+1}")
            ratios = np.sort(np.unique([ r for r in ratios if r>0 ]))
            kd = st.number_input('Kd (nM)', value=1000.0, min_value=0.0, step=1.0, format="%g", help="Binding affinity")
            kd_molar = kd * 1e-9    # nM -> M
            pmin = st.number_input('[P] (mg/ml) - min', value=0.001, min_value=0.0, step=0.01, format="%g", help="Minimal total protein concentration to plot")
            pmax = st.number_input('[P] (mg/ml) - max', value=10.0, min_value=pmin*10.0, step=0.1, format="%g", help="Maximal total protein concentration to plot")
        else:   # "Kd"
            ratio = st.number_input(f'[P]:[L]=1:?', value=1.0, min_value=0.0, step=1.0, format="%g", help="Molar ratio of total protein and total ligand")
            n = int(st.number_input('# of Kd', value=3, min_value=1, step=1, help="Number of binding affinities to plot"))
            kd_list = [1e-6, 1e-9, 1e-3, 1e-12, 1e-7, 1e-5, 1e-8, 1e-4, 1e-15, 1e-13, 1e-14]
            kds = [kd_list[i % len(kd_list)] * 1e9 for i in range(n)]
            for i, kd in enumerate(sorted(kds)):
                kds[i] = st.number_input(f'Kd (nM) - Curve {i+1}', value=kd, min_value=0.0, step=kd*1e-1, format="%g", help=f"Binding affinity for curve {i+1}")
            kds = np.sort(np.unique([ kd for kd in kds if kd>0 ]))
            kds_molar = kds * 1e-9    # nM -> M
            pmin = st.number_input('[P] (mg/ml) - min', value=0.001, min_value=0.0, step=0.01, format="%g", help="Minimal total protein concentration to plot")
            pmax = st.number_input('[P] (mg/ml) - max', value=10.0, min_value=pmin*10.0, step=0.1, format="%g", help="Maximal total protein concentration to plot")
        n_data_points = int(st.number_input('# of data points', value=100, min_value=10, step=10, help="Number of data points per curve"))
        xlog = st.checkbox('X-axis in log scale', value=True)
        

    with col2:
        y_label = "fraction of protein in ligand-bound state"
        x_axis_type = "log" if xlog else "linear"
        tools = 'box_zoom,crosshair,hover,pan,reset,save,wheel_zoom'
        line_dashes = 'solid dashed dotted dotdash dashdot'.split()
        legends = []
        raw_data = []
        if mode == "[P]":
            if xlog:
                l = np.logspace(np.log10(lmin), np.log10(lmax), n_data_points+1)
            else:
                l = np.linspace(lmin, lmax, n_data_points+1)
            l_molar = l * 1e-6  # μM -> M

            x = l
            x_label = r"$$[L]_T\ (μM)$$"

            hover_tips = [("[P]:[L]", "1:@ratio"), ("[P]", "@p_txt"), ("[L]", "$x μM"), ("Kd", "@kd nM"), ("fraction", "$y")]
            fig = figure(title="", x_axis_type=x_axis_type, x_axis_label=x_label, y_axis_label=y_label, tools=tools, tooltips=hover_tips)

            for pi, p in enumerate(ps):
                p_molar = p/(m*1e3) # mg/ml -> M
                p_txt = f"{p:.3g} mg/ml / {p_molar*1e6:.3g} μM"
                ratio = l_molar/p_molar
                f = binding_fraction(kd=kd_molar, protein_concentration=p_molar, ligand_concentration=l_molar)
                source = dict(x=x, y=f, p_txt=[p_txt]*len(x), ratio=ratio, kd=[kd]*len(x))
                line_dash = line_dashes[pi%len(line_dashes)]
                line_width = 2
                line = fig.line(x='x', y='y', source=source, line_dash=line_dash, line_width=line_width)
                label = f"[P] = {p_txt}"
                legends.append(LegendItem(label=label, renderers=[line]))
                raw_data.append((label, f))        
            fig.x_range.start = l[0]
            fig.x_range.end = l[-1]
        elif mode == "[P]:[L]":
            if xlog:
                p = np.logspace(np.log10(pmin), np.log10(pmax), n_data_points+1)
            else:
                p = np.linspace(min, pmax, n_data_points+1)
            p_molar = p / (m*1e3)  # mg/ml -> M

            x = p
            x_label = r"$$[P]_T\ (mg/ml)$$"

            hover_tips = [("[P]:[L]", "@ratio_txt"), ("[P]", "$x mg/ml / @x_mc μM"), ("[L]", "@l μM"), ("Kd", "@kd nM"), ("fraction", "$y")]
            fig = figure(title="", x_axis_type=x_axis_type, x_axis_label=x_label, y_axis_label=y_label, tools=tools, tooltips=hover_tips)

            for ri, r in enumerate(ratios):
                ratio_txt = f"1:{r:g}"
                l_molar = p_molar * r
                f = binding_fraction(kd=kd_molar, protein_concentration=p_molar, ligand_concentration=l_molar)
                source = dict(x=x, y=f, x_mc=p_molar*1e6, ratio_txt=[ratio_txt]*len(x), l=l_molar*1e6, kd=[kd]*len(x))
                line_dash = line_dashes[ri%len(line_dashes)]
                line_width = 2
                line = fig.line(x='x', y='y', source=source, line_dash=line_dash, line_width=line_width)
                label = f"[P]:[L] = {ratio_txt}"
                legends.append(LegendItem(label=label, renderers=[line]))
                raw_data.append((label, f))
            fig.x_range.start = p[0]
            fig.x_range.end = p[-1]        
        else:   # "Kd"
            if xlog:
                p = np.logspace(np.log10(pmin), np.log10(pmax), n_data_points+1)
            else:
                p = np.linspace(min, pmax, n_data_points+1)
            p_molar = p / (m*1e3)  # mg/ml -> M

            x = p
            x_label = r"$$[P]_T\ (mg/ml)$$"

            hover_tips = [("[P]:[L]", "@ratio_txt"), ("[P]", "$x mg/ml / @x_mc μM"), ("[L]", "@l μM"), ("Kd", "@kd_txt"), ("fraction", "$y")]
            fig = figure(title="", x_axis_type=x_axis_type, x_axis_label=x_label, y_axis_label=y_label, tools=tools, tooltips=hover_tips)

            for ki, kd_molar in enumerate(kds_molar):
                kd_txt = f"{kd_molar*1e9:.3g} nM"
                ratio_txt = f"1:{ratio:g}"
                l_molar = p_molar * ratio
                f = binding_fraction(kd=kd_molar, protein_concentration=p_molar, ligand_concentration=l_molar)
                source = dict(x=x, y=f, x_mc=p_molar*1e6, ratio_txt=[ratio_txt]*len(x), l=l_molar*1e6, kd_txt=[kd_txt]*len(x))
                line_dash = line_dashes[ki%len(line_dashes)]
                line_width = 2
                line = fig.line(x='x', y='y', source=source, line_dash=line_dash, line_width=line_width)
                label = f"Kd = {kd_txt}"
                legends.append(LegendItem(label=label, renderers=[line]))
                raw_data.append((label, f))
            fig.x_range.start = p[0]
            fig.x_range.end = p[-1]        
        fig.y_range.start = 0
        fig.y_range.end = 1
        fig.yaxis[0].ticker.desired_num_ticks = 10
        from bokeh.models import CustomJS
        from bokeh.events import MouseEnter
        title_js = CustomJS(args=dict(title=title), code="""
            document.title=title
        """)
        fig.js_on_event(MouseEnter, title_js)
        if len(legends):
            from bokeh.models import Legend
            legend = Legend(items=legends)
            fig.add_layout(legend)
            fig.legend[0].location = "top_left"
            fig.legend.click_policy= "hide"
            from bokeh.models import CustomJS
            from bokeh.events import DoubleTap
            toggle_legend_js = CustomJS(args=dict(leg=fig.legend[0]), code="""
                if (leg.visible) {
                    leg.visible = false
                    }
                else {
                    leg.visible = true
                }
            """)
            fig.js_on_event(DoubleTap, toggle_legend_js)
        st.text("") # workaround for a layout bug in streamlit 
        st.bokeh_chart(fig, use_container_width=True)

        show_data = st.checkbox('show raw data', value=False)
        if show_data:
            import pandas as pd
            data = np.zeros((len(x), 1+len(raw_data)))
            data[:, 0] = x
            columns = [x_label]
            for i, (label, f) in enumerate(raw_data):
                data[:, 1+i] = f
                columns.append(('fraction-'+label))
            columns = [col.replace(" ", "").rjust(32) for col in columns]
            df = pd.DataFrame(data, columns=columns)
            st.dataframe(df, width=None)
            st.markdown(get_table_download_link(df), unsafe_allow_html=True)

        st.markdown("*Developed by the [Jiang Lab@Purdue University](https://jiang.bio.purdue.edu). Report problems to Wen Jiang (jiang12 at purdue.edu)*")

def binding_fraction(kd, protein_concentration, ligand_concentration):
    # unit of all concentrations: M
    t = protein_concentration + ligand_concentration + kd
    f = (t - np.sqrt(t*t - 4*protein_concentration*ligand_concentration))/(2*protein_concentration)
    return f

def get_table_download_link(df):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=False)
    import base64
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}" download="raw_data.csv">Download the raw data</a>'
    return href

@st.cache(persist=True, show_spinner=False)
def setup_anonymous_usage_tracking():
    try:
        import pathlib, stat
        index_file = pathlib.Path(st.__file__).parent / "static/index.html"
        index_file.chmod(stat.S_IRUSR|stat.S_IWUSR|stat.S_IRGRP|stat.S_IROTH)
        txt = index_file.read_text()
        if txt.find("gtag/js?")==-1:
            txt = txt.replace("<head>", '''<head><script async src="https://www.googletagmanager.com/gtag/js?id=G-6ZB67M09Q9"></script><script>window.dataLayer = window.dataLayer || [];function gtag(){dataLayer.push(arguments);}gtag('js', new Date());gtag('config', 'G-6ZB67M09Q9');</script>''')
            index_file.write_text(txt)
    except:
        pass

if __name__ == "__main__":
    setup_anonymous_usage_tracking()
    main()
