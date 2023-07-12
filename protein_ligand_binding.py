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
from ssl import PROTOCOL_TLS_CLIENT
import streamlit as st
import numpy as np
import scipy
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
from sympy import roots, cos, nsolve, Symbol, nroots, solve, poly
from sympy.abc import H




def main():
    title = "Single-Ligand Multi-Protein Binding"
    st.set_page_config(page_title=title, layout="wide")

    hosted, host = is_hosted(return_host=True)
    if hosted and host in ['heroku']:
        st.error(f"This app hosted on Heroku will be unavailable starting November 28, 2022 [when Heroku discontinues free hosting service](https://blog.heroku.com/next-chapter). Please switch to [the same app hosted elsewhere](https://jianglab-protein-ligand-binding-streamlit-app-7qd88x.streamlitapp.com)")

    with st.sidebar:    # sidebar at the left of the screen
        st.title(title)
        st.write(r'For a competitive binding reaction:$\\P_{A} + L\  {\rightleftharpoons}\  P_{A}{\cdot}L$ $\\P_{B} + L\  {\rightleftharpoons}\  P_{B}{\cdot}L$')
        st.write(r'the binding affinities, $K_A$ and $K_B$, are defined as  $K_D = \frac{[P]_{free}[L]_{free}}{[P {\cdot}L]}$.')
        st.write(r'Then the fraction of protein, $P_{A}$, in a competitive ligand-bound state is:')
        st.latex(r'[P_A{\cdot}L] = \frac{[L_{free}]{\cdot}[P_{A_{0}}]}{K_{A}{\cdot}(1+\frac{[P_{B_{0}}]}{K_B})+[L]_{free}}')
        st.write(r'Note that $[L]_{free}$ means free, unbound ligand concentration that we cannot control directly. Instead, we can only control the total concentration of the protein ($P_{A_{0}}=[P_{A}{\cdot}L]+[P]_{A_{free}},$ and $P_{B_{0}}=[P_{B}{\cdot}L]+[P]_{B_{free}}$) and ligand ($L_0=[P_{A}{\cdot}L]+[P_{B}{\cdot}L]+[L]_{free}$). So to generalize this to $n$ number of proteins we see:')
        st.latex(r'P_{A_{0}}=[P_{A}{\cdot}L]+[P]_{A_{free}}')
        st.latex(r'L_0=[P_{A}{\cdot}L]+[P_{B}{\cdot}L]+[P_{C}{\cdot}L]+...+[P_{n}{\cdot}L]+[L]_{free}')
        st.write(r'So we can numerically solve for $L_{free}$. Plugging the estimated value of $L_{free}$ back into our fraction bound equation, we can approximate a solution.')
        
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

        num_p1 = num_p+1

        column_l = ['Ligand']
        df_l = pd.DataFrame(np.tile([1000], (1, 1)), ['Lig. Conc. (μM)'], column_l)
        edited_df_l = st.data_editor(df_l)

        column_p = ['Protein %d' % (i+1) for i in range(int(num_p))]
        df_p = pd.DataFrame(np.tile([300], (1, int(num_p))), ['Prot. Conc. (μM)'], column_p)
        edited_df_p = st.data_editor(df_p) 


        column_kd = ['Protein %d Kd (μM)' % (i+1) for i in range(int(num_p))]
        df_kd = pd.DataFrame(np.tile([500], (1, int(num_p))), [f'Ligand'], column_kd)
        edited_df_kd = st.data_editor(df_kd) 

        xlog = st.checkbox('X-axis in log scale', value=True)


        



    with col2:
        y_label = "Fraction of protein in ligand-bound state"
        x_axis_type = "log" if xlog else "linear"
        tools = 'box_zoom,crosshair,hover,pan,reset,save,wheel_zoom'
        line_dashes = 'solid dashed dotted dotdash dashdot '.split()
        legends = []
        raw_data = []
        i = 0
        ratios = [1.0] * int(num_p)


          
        lmax = edited_df_l.loc['Lig. Conc. (μM)','Ligand']
        if xlog:
            l = np.logspace(np.log10(1e-3), np.log10(lmax), 101)
        else:
            l = np.linspace(1e-3, lmax, 101)
        l_molar = l / (1e3)  # mg/ml -> M

        x = l
        x_label = r"$$[L]_T\ (μM)$$"

        hover_tips = [("[P]:[L]", "@ratio_txt"), ("[P]", "$x mg/ml / @x_mc μM"), ("[L]", "@l μM"), ("Kd", "@kd nM"), ("fraction", "$y")]
        fig = figure(title="", x_axis_type=x_axis_type, x_axis_label=x_label, y_axis_label=y_label, tools=tools, tooltips=hover_tips)

            
        for ri, r in enumerate(ratios):  #plot
            frac = H
            i=i+1
            l_mol = edited_df_l.loc['Lig. Conc. (μM)','Ligand']
            p_mol = edited_df_p.loc['Prot. Conc. (μM)', f'Protein {i}'] 
            kd_molar = edited_df_kd.loc['Ligand', f'Protein {i} Kd (μM)'] 
            kd_list_molar = pd.concat([edited_df_kd.iloc[0][0:i-1],edited_df_kd.iloc[0][i:]])
            p_mol_list = pd.concat([edited_df_p.iloc[0][0:i-1],edited_df_p.iloc[0][i:]])
            ratio_txt = f"{edited_df_p.loc['Prot. Conc. (μM)', f'Protein {i}']}:{edited_df_l.loc['Lig. Conc. (μM)','Ligand']}"
            frac += percent_bound(kd=kd_molar, kd_list=kd_list_molar, prot_conc=p_mol, prot_conc_list=p_mol_list, H=H)
            sol_l=[]
            for j in l:
                tmp=frac-j
                sol_l.append(nsolve(tmp, H, 1))
            sol_l=np.array(sol_l,dtype=np.float64)
            total = (percent_bound(kd=kd_molar, kd_list=kd_list_molar, prot_conc=p_mol, prot_conc_list=p_mol_list, H = sol_l*1000))/(edited_df_p.loc['Prot. Conc. (μM)', f'Protein {i}'])
            source = dict(x=x, y=total, x_mc=l_molar*r*1e6, ratio_txt=[ratio_txt]*len(x), l=l_molar*1e6, kd=[kd_molar]*len(x))
            line_dash = line_dashes[ri%len(line_dashes)]
            line = fig.line(x='x', y='y', source=source, line_dash=line_dash, line_width=2)
            label = f"Protein {i}"
            legends.append(LegendItem(label=label, renderers=[line]))
            raw_data.append((label, total))
        fig.x_range.start = l[0]
        fig.x_range.end = l[-1]
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
 

def percent_bound(kd, kd_list, prot_conc, prot_conc_list, H):
    adj_kd = kd
    for i in range(len(kd_list)):
        adj_kd = adj_kd*(1+prot_conc_list[i]/kd_list[i])
    total = (prot_conc*H)/(adj_kd+H)
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