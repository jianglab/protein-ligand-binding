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

from secrets import randbelow
import streamlit as st
import numpy as np
from bokeh.plotting import figure
from bokeh.models import LegendItem 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from bokeh.models import CustomJS
from bokeh.events import MouseEnter
from bokeh.models import Legend
from bokeh.models import CustomJS
from bokeh.events import DoubleTap




def main():
    title = "Multi-Protein Multi-Ligand Binding"
    st.set_page_config(page_title=title, layout="wide")

    hosted, host = is_hosted(return_host=True)
    if hosted and host in ['heroku']:
        st.error(f"This app hosted on Heroku will be unavailable starting November 28, 2022 [when Heroku discontinues free hosting service](https://blog.heroku.com/next-chapter). Please switch to [the same app hosted elsewhere](https://jianglab-protein-ligand-binding-streamlit-app-7qd88x.streamlitapp.com)")

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
        
        hide_streamlit_style = """
        <style>
        #MainMenu {visibility: hidden;}
        footer {visibility: hidden;}
        </style>
        """
        st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

    # right-side main panel
    col1, col2 = st.columns((2, 4))


    with col1:
        
        num_p = st.number_input('Number of Proteins', value=2, min_value=2, step=1)
        num_l = st.number_input('Number of Ligands', value=2, min_value=2, step=1)
        num_p1 = num_p+1
        num_l1 = num_l+1
        mode = st.radio(label="Ligand", options=list(reversed(range(int(num_l1))))[:-1], index=1, horizontal=True, help="Choose which ligand to plot.")

        column_p = ['Protein %d' % (i+1) for i in range(int(num_p))]
        df_p = pd.DataFrame(np.tile([0.6], (1, int(num_p))), ['Prot. Conc. (μM)'], column_p)
        edited_df_p = st.experimental_data_editor(df_p) 

        column_l = ['Ligand %d' % (i+1) for i in range(int(num_l))]
        df_l = pd.DataFrame(np.tile([0.6], (1, int(num_l))), ['Lig. Conc. (μM)'], column_l)
        edited_df_l = st.experimental_data_editor(df_l)

        column_kd = ['Protein %d Kd' % (i+1) for i in range(int(num_p))]
        df_kd = pd.DataFrame(np.tile([0.6], (1, int(num_p))), [f'Ligand {mode}'], column_kd)
        edited_df_kd = st.experimental_data_editor(df_kd) 

        xlog = st.checkbox('X-axis in log scale', value=True)


        



    with col2:
        y_label = "fraction of protein in ligand-bound state"
        x_axis_type = "log" if xlog else "linear"
        tools = 'box_zoom,crosshair,hover,pan,reset,save,wheel_zoom'
        line_dashes = 'solid dashed dotted dotdash dashdot '.split()
        legends = []
        raw_data = []
        n_data_points = 100
        lmax = 1000
        lmin = 0.0001
        pmin = 0.001
        i = 0
        m = 100
        l_molar = 1000
        n = int(num_p)
        ratios = [1.0] * n


        if mode == 1:
            pmax = edited_df_l.loc['Lig. Conc. (μM)', 'Ligand 1']
            if xlog:
                l = np.logspace(np.log10(pmin), np.log10(pmax), n_data_points+1)
            else:
                l = np.linspace(pmin, pmax, n_data_points+1)
            l_molar = l / (m*1e3)  # mg/ml -> M

            x = l
            x_label = r"$$[L]_T\ (μM)$$"
            y = int(edited_df_l.loc['Lig. Conc. (μM)', 'Ligand 1'])

            hover_tips = [("[P]:[L]", "@ratio_txt"), ("[P]", "$x mg/ml / @x_mc μM"), ("[L]", "@l μM"), ("Kd", "@kd nM"), ("fraction", "$y")]
            fig = figure(title="", x_axis_type=x_axis_type, x_axis_label=x_label, y_axis_label=y_label, tools=tools, tooltips=hover_tips)

            for ri, r in enumerate(ratios):  #plot
                i=i+1
                kd = edited_df_kd.loc['Ligand 1', f'Protein {i} Kd']
                kd_molar = kd 
                kd_b_molar = edited_df_kd.loc['Ligand 1', f'Protein {(i%2)+1} Kd'] 
                p_mol = edited_df_p.loc['Prot. Conc. (μM)', f'Protein {i}'] 
                p_mol_b = edited_df_p.loc['Prot. Conc. (μM)', f'Protein {(i%2)+1}']
                l_mol = edited_df_l.loc['Lig. Conc. (μM)','Ligand 1'] 
                ratio_txt = f"{edited_df_p.loc['Prot. Conc. (μM)', f'Protein {i}']}:{edited_df_l.loc['Lig. Conc. (μM)','Ligand 1']}"
                p_molar = l_molar * r
                if int(num_p) == 2:
                    total = (percent_bound(kd=kd_molar, kd_b=kd_b_molar, prot_conc=p_mol, prot_conc_b=p_mol_b, lig_conc=l))/(edited_df_kd.loc['Ligand 1', f'Protein {i} Kd'])

                elif int(num_p) == 3: 
                        total
                else:
                    return
                source = dict(x=x, y=total, x_mc=p_molar*1e6, ratio_txt=[ratio_txt]*len(x), l=l_molar*1e6, kd=[kd]*len(x))
                line_dash = line_dashes[ri%len(line_dashes)]
                fig.segment(l_mol, p_mol, l_mol, p_mol, color="lightgrey", line_width=2)
                line = fig.line(x='x', y='y', source=source, line_dash=line_dash, line_width=2)
                label = f"[P{i}]:[L{mode}] = {ratio_txt}"
                legends.append(LegendItem(label=label, renderers=[line]))
                raw_data.append((label, total))
            fig.x_range.start = l[0]
            fig.x_range.end = l[-1]  
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
        title_js = CustomJS(args=dict(title=title), code="""
            document.title=title
        """)
        fig.js_on_event(MouseEnter, title_js)
        if len(legends):
            legend = Legend(items=legends)
            fig.add_layout(legend)
            fig.legend[0].location = "top_left"
            fig.legend.click_policy= "hide"
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


def percent_bound(kd, kd_b, prot_conc, prot_conc_b, lig_conc):
    #total = (prot_conc*lig_conc)/(kd*(1+(prot_conc_b/kd_b))+prot_conc)
    total = (prot_conc*lig_conc)/((kd*(1+(prot_conc_b/kd_b)))+lig_conc)
    return total

def get_username():
    from getpass import getuser
    return getuser()

def get_hostname():
    import socket
    fqdn = socket.getfqdn()
    return fqdn

def is_hosted(return_host=False):
    hosted = False
    host = ""
    fqdn = get_hostname()
    if fqdn.find("heroku")!=-1:
        hosted = True
        host = "heroku"
    username = get_username()
    if username.find("appuser")!=-1:
        hosted = True
        host = "streamlit"
    if not host:
        host = "localhost"
    if return_host:
        return hosted, host
    else:
        return hosted

def qr_code(url=None, size = 8):
    import_with_auto_install(["qrcode"])
    import qrcode
    if url is None: # ad hoc way before streamlit can return the url
        _, host = is_hosted(return_host=True)
        if len(host)<1: return None
        if host == "streamlit":
            url = "https://share.streamlit.io/wjiang/ctfsimulation/master/"
        elif host == "heroku":
            url = "https://ctfsimulation.herokuapp.com/"
        else:
            url = f"http://{host}:8501/"
        import urllib
        params = st.experimental_get_query_params()
        d = {k:params[k][0] for k in params}
        url += "?" + urllib.parse.urlencode(d)
    if not url: return None
    img = qrcode.make(url)  # qrcode.image.pil.PilImage
    data = np.array(img.convert("RGBA"))
    return data

if __name__ == "__main__":
    main()