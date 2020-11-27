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

col1, col2 = st.beta_columns((2,1))
col1.title("Thermal Storage Tank")
with col2:
    col2.image("static/logo.png",use_column_width=True,output_format='PNG')

with st.beta_expander("Tool description"):
    st.markdown("A simple simulation tool to simulate stratified thermal storage tank.")

