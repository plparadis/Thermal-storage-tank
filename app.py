# Author PLParadis
import streamlit as st
from matplotlib import pyplot as plt


TextColor = "#404041"
plt.rcParams['text.color'] = TextColor
plt.rcParams['axes.labelcolor'] = TextColor
plt.rcParams['xtick.color'] = TextColor
plt.rcParams['ytick.color'] = TextColor
fontname = 'Helvetica LT Std'

st.set_page_config(page_title="Thermal Storage Tank",
                    page_icon="static/favicon.ico",
                    layout="centered",
                    initial_sidebar_state="auto")

col1, col2 = st.beta_columns((3,1))
col1.title("Thermal Storage Tank")
with col2:
    col2.image("static/logo.png",use_column_width=True,output_format='PNG')

with st.beta_expander("Tool description",expanded=True):
    st.markdown("A simple simulation tool to simulate stratified thermal storage tank as seen below.")
    st.image("static/Schematic.JPG",use_column_width=True,output_format='JPG')
    
    st.markdown("A heat balance on a slice of the tank gives")
    st.write(r'''
        $$
        \underbrace{(rho*c_p)_{SHW} * 
        \frac{\partial T_{SHW}}{\partial t}}_{\text{Transient term}}
        =
        \underbrace{\frac{\partial }{\partial y} \lbrack (k_{SHW}+\Delta k) \frac{\partial T_{SHW}}{\partial y} \rbrack}_{\text{1-D conduction term}} -
        \underbrace{\frac{(\dot{m}c_p)_{SHW}}{A_{SHW}}\frac{\partial T_{SHW}}{\partial y}}_{\text{Advection term}} +
        $$
        ''')
    st.write(r'''
        $$
        \underbrace{\frac{dq_{DCW} - dq_{loss}}{d\forall_{SHW}}}_{\text{Source terms}}     
        $$
        ''')




